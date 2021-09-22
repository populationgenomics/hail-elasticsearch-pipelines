"""
Functions to find the pipeline inputs and communicate with the SM server
"""

import logging
import sys
from dataclasses import dataclass
from typing import List, Dict, Optional, Set, Collection

import pandas as pd
from hailtop.batch import Batch
from hailtop.batch.job import Job
from sample_metadata import (
    AnalysisApi,
    SequenceApi,
    SampleApi,
    exceptions,
    AnalysisUpdateModel,
    AnalysisModel,
)

import utils


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


sapi = SampleApi()
aapi = AnalysisApi()
seqapi = SequenceApi()


@dataclass
class Analysis:
    """
    Represents the analysis DB entry
    """

    id: str
    type: str
    status: str
    sample_ids: Set[str]
    output: Optional[str]


def _parse_analysis(data: Dict) -> Optional[Analysis]:
    if not data:
        return None
    if 'id' not in data:
        logger.error(f'Analysis data doesn\'t have id: {data}')
        return None
    if 'type' not in data:
        logger.error(f'Analysis data doesn\'t have type: {data}')
        return None
    if 'status' not in data:
        logger.error(f'Analysis data doesn\'t have status: {data}')
        return None
    a = Analysis(
        id=data['id'],
        type=data['type'],
        status=data['status'],
        sample_ids=set(data.get('sample_ids', [])),
        output=data.get('output', None),
    )
    return a


def process_existing_analysis(
    sample_ids: Collection[str],
    completed_analysis: Optional[Analysis],
    analysis_type: str,
    analysis_project: str,
    analysis_sample_ids: Collection[str],
    expected_output_fpath: str,
    skip_stage: bool,
) -> Optional[str]:
    """
    Checks whether existing analysis exists, and output matches the expected output
    file. Invalidates bad analysis by setting status=failure, and submits a
    status=completed analysis if the expected output already exists.

    Returns the path to the output if it can be reused, otherwise None.

    :param sample_ids: sample IDs to pull the analysis for
    :param completed_analysis: existing completed analysis of this type for these samples
    :param analysis_type: cram, gvcf, joint_calling
    :param analysis_project: analysis project name (e.g. seqr)
    :param analysis_sample_ids: sample IDs that analysis refers to
    :param expected_output_fpath: where the pipeline expects the analysis output file
        to sit on the bucket (will invalidate the analysis if it doesn't match)
    :param skip_stage: if not skip_stage and analysis output is not as expected,
        we invalidate the analysis and set its status to failure
    :return: path to the output if it can be reused, otherwise None
    """
    label = f'type={analysis_type}'
    if len(analysis_sample_ids) > 1:
        label += f' for {", ".join(analysis_sample_ids)}'

    found_output_fpath = None
    if not completed_analysis:
        logger.warning(
            f'Not found completed analysis {label} for '
            f'{f"sample {sample_ids}" if len(sample_ids) == 1 else f"{len(sample_ids)} samples" }'
        )
    elif not completed_analysis.output:
        logger.error(
            f'Found a completed analysis {label}, '
            f'but the "output" field does not exist or empty'
        )
    else:
        found_output_fpath = str(completed_analysis.output)
        if found_output_fpath != expected_output_fpath:
            logger.error(
                f'Found a completed analysis {label}, but the "output" path '
                f'{found_output_fpath} does not match the expected path '
                f'{expected_output_fpath}'
            )
            found_output_fpath = None
        elif not utils.file_exists(found_output_fpath):
            logger.error(
                f'Found a completed analysis {label}, '
                f'but the "output" file {found_output_fpath} does not exist'
            )
            found_output_fpath = None

    # skipping stage
    if skip_stage:
        if found_output_fpath:
            logger.info(f'Skipping stage, picking existing {found_output_fpath}')
            return found_output_fpath
        else:
            logger.info(
                f'Skipping stage, and expected {expected_output_fpath} not found, '
                f'so skipping {label}'
            )
            return None

    # completed and good exists, can reuse
    if found_output_fpath:
        logger.info(
            f'Completed analysis {label} exists, '
            f'reusing the result {found_output_fpath}'
        )
        return found_output_fpath

    # can't reuse, need to invalidate
    if completed_analysis:
        logger.warning(
            f'Invalidating the analysis {label} by setting the status to "failure", '
            f'and resubmitting the analysis.'
        )
        aapi.update_analysis_status(
            completed_analysis.id, AnalysisUpdateModel(status='failed')
        )

    # can reuse, need to create a completed one?
    if utils.file_exists(expected_output_fpath):
        logger.info(
            f'Output file {expected_output_fpath} already exists, so creating '
            f'an analysis {label} with status=completed'
        )
        am = AnalysisModel(
            type=analysis_type,
            output=expected_output_fpath,
            status='completed',
            sample_ids=analysis_sample_ids,
        )
        aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        return expected_output_fpath

    # proceeding with the standard pipeline (creating status=queued, submitting jobs)
    else:
        logger.info(
            f'Expected output file {expected_output_fpath} does not exist, '
            f'so queueing analysis {label}'
        )
        return None


def find_joint_calling_analysis(
    analysis_project: str,
    sample_ids: Collection[str],
) -> Optional[Analysis]:
    """
    Query the DB to find the last completed joint-calling analysis for the samples
    """
    data = aapi.get_latest_complete_analysis_for_type(
        project=analysis_project,
        analysis_type='joint-calling',
    )
    a = _parse_analysis(data)
    if not a:
        return None
    assert a.type == 'joint-calling', data
    assert a.status == 'completed', data
    if a.sample_ids != set(sample_ids):
        return None
    return a


def find_analyses_by_sid(
    sample_ids: Collection[str],
    analysis_project: str,
    analysis_type: str,
) -> Dict[str, Analysis]:
    """
    Query the DB to find the last completed analysis for the type and samples,
    one Analysis object per sample. Assumes the analysis is defined for a single
    sample (e.g. cram, gvcf)
    """
    analysis_per_sid: Dict[str, Analysis] = dict()
    datas = aapi.get_latest_analysis_for_samples_and_type(
        project=analysis_project,
        analysis_type=analysis_type,
        request_body=sample_ids,
    )
    for data in datas:
        a = _parse_analysis(data)
        if not a:
            continue
        assert a.type == analysis_type, data
        assert a.status == 'completed', data
        assert len(a.sample_ids) == 1, data
        analysis_per_sid[list(a.sample_ids)[0]] = a
    return analysis_per_sid


def make_sm_in_progress_job(
    b: Batch, analyais_type: str, analysis_project: str, analysis_id: str
) -> Job:
    """
    Creates a job that updates the sample metadata server entry analysis status
    to in-progress
    """
    return make_sm_update_status_job(
        b, analyais_type, 'in-progress', analysis_project, analysis_id
    )


def make_sm_completed_job(
    b: Batch, analyais_type: str, sm_db_name: str, analysis_id: str
) -> Job:
    """
    Creates a job that updates the sample metadata server entry analysis status
    to completed
    """
    return make_sm_update_status_job(
        b, analyais_type, 'completed', sm_db_name, analysis_id
    )


def make_sm_update_status_job(
    b: Batch, analysis_type: str, status: str, sm_db_name: str, analysis_id: str
) -> Job:
    """
    Creates a job that updates the sample metadata server entry analysis status.
    """
    assert status in ['in-progress', 'failed', 'completed', 'queued']
    j = b.new_job(f'SM: update {analysis_type} to {status}')
    j.image(utils.SM_IMAGE)
    j.command(
        f"""
set -o pipefail
set -ex

export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS
export SM_DEV_DB_PROJECT={sm_db_name}
export SM_ENVIRONMENT=PRODUCTION

cat <<EOT >> update.py
from sample_metadata.api import AnalysisApi
from sample_metadata import AnalysisUpdateModel
aapi = AnalysisApi()
aapi.update_analysis_status(
    analysis_id='{analysis_id}',
    analysis_update_model=AnalysisUpdateModel(status='{status}'),
)
EOT
python update.py
    """
    )
    return j


def replace_paths_to_test(s: Dict) -> Optional[Dict]:
    """
    Replace paths of all files in -main namespace to -test namespsace,
    and return None if files in -test are not found.
    :param s:
    :return:
    """

    def fix(fpath):
        fpath = fpath.replace('-main-upload/', '-test-upload/')
        if not utils.file_exists(fpath):
            return None
        return fpath

    try:
        reads_type = s['meta']['reads_type']
        if reads_type in ('bam', 'cram'):
            fpath = s['meta']['reads'][0]['location']
            fpath = fix(fpath)
            if not fpath:
                return None
            s['meta']['reads'][0]['location'] = fpath

            fpath = s['meta']['reads'][0]['secondaryFiles'][0]['location']
            fpath = fix(fpath)
            if not fpath:
                return None
            s['meta']['reads'][0]['secondaryFiles'][0]['location'] = fpath

        elif reads_type == 'fastq':
            for li in range(len(s['meta']['reads'])):
                for rj in range(len(s['meta']['reads'][li])):
                    fpath = s['meta']['reads'][li][rj]['location']
                    fpath = fix(fpath)
                    if not fpath:
                        return None
                    s['meta']['reads'][li][rj]['location'] = fpath

        logger.info(f'Found test sample {s["id"]}')
        return s
    except Exception:  # pylint: disable=broad-except
        return None


default_entry = {
    's': None,
    'external_id': None,
    'project': None,
    'continental_pop': '-',
    'subpop': '-',
    'gvcf': '-',
    'cram': '-',
    'crai': '-',
    'batch': '-',
    'operation': 'add',
    'flowcell_lane': '-',
    'library_id': '-',
    'platform': '-',
    'centre': '-',
    'r_contamination': None,
    'r_chimera': None,
    'r_duplication': None,
    'median_insert_size': None,
}


def find_inputs_from_db(
    input_projects: List[str],
    is_test: bool = False,
    skip_samples: Optional[Collection[str]] = None,
) -> pd.DataFrame:
    """
    Determine input samples and pull input files and metadata from
    the CPG sample-metadata server database.
    """
    inputs = []

    for proj in input_projects:
        logger.info(f'Processing project {proj}')
        samples = sapi.get_samples(
            body_get_samples_by_criteria_api_v1_sample_post={
                'project_ids': [proj],
                'active': True,
            }
        )
        logger.info(f'Found {len(samples)} samples')
        if not samples:
            logger.info(f'No samples to process, skipping project {proj}')
            continue

        if skip_samples:
            logger.info('Checking which samples need to skip')
            not_skipped_sids = []
            for s in samples:
                if skip_samples and s['id'] in skip_samples:
                    logger.info(f'Skiping sample: {s["id"]}')
                    continue
                not_skipped_sids.append(s['id'])
            logger.info(f'Excluding skipped samples: {len(not_skipped_sids)}')
            samples = [s for s in samples if s in not_skipped_sids]
            if not samples:
                logger.info(f'No samples to process, skipping project {proj}')
                continue

        logger.info('Checking GVCF analyses for samples')
        gvcf_analysis_per_sid = find_analyses_by_sid(
            sample_ids=[s['id'] for s in samples],
            analysis_type='gvcf',
            analysis_project=proj,
        )
        gvcf_by_sid = dict()
        sids_without_gvcf = []
        for s in samples:
            a = gvcf_analysis_per_sid.get(s['id'])
            if not a:
                sids_without_gvcf.append(s['id'])
                continue
            gvcf_path = a.output
            if not gvcf_path:
                logger.error(
                    f'"output" is not defined for the latest gvcf analysis, '
                    f'skipping sample {s["id"]}'
                )
                sids_without_gvcf.append(s['id'])
                continue
            gvcf_by_sid[s['id']] = gvcf_path

        if sids_without_gvcf:
            logger.warning(
                f'No gvcf found for {len(sids_without_gvcf)}/{len(samples)} samples: '
                f'{", ".join(sids_without_gvcf)}'
            )
            samples = [s for s in samples if s['id'] in gvcf_by_sid]
            if not samples:
                logger.info(f'No samples to process, skipping project {proj}')
                continue

        logger.info('Checking sequencing info for samples')
        try:
            seq_infos: List[Dict] = seqapi.get_sequences_by_sample_ids(
                request_body=[s['id'] for s in samples]
            )
        except exceptions.ApiException:
            logger.critical(f'Not for all samples sequencing data was found')
            raise

        seq_meta_by_sid: Dict = dict()
        sids_without_meta = []
        for si in seq_infos:
            if 'meta' not in si:
                sids_without_meta.append(si['sample_id'])
            else:
                seq_meta_by_sid[si['sample_id']] = si['meta']
        if sids_without_meta:
            logger.error(
                f'Found {len(sids_without_meta)} samples without "meta" in '
                f'sequencing info: {", ".join(sids_without_meta)}'
            )
        samples = [s for s in samples if s['id'] in seq_meta_by_sid]
        if not samples:
            logger.info(f'No samples to process, skipping project {proj}')
            continue

        logger.info(
            'Checking GVCFs for samples and collecting pipeline input data frame'
        )
        for s in samples:
            sample_id = s['id']
            external_id = s['external_id']
            seq_meta = seq_meta_by_sid[sample_id]
            gvcf_path = gvcf_by_sid[s['id']]

            # TODO: reenable once we support raw data and crams
            # if is_test:
            #     s = replace_paths_to_test(s)
            # if s:
            #     samples_by_project[proj].append(s)

            if is_test:
                if '/batch1/' not in gvcf_path:
                    continue
                gvcf_path = gvcf_path.replace(
                    f'gs://cpg-{proj}-main',
                    f'gs://cpg-{proj}-test',
                )
                gvcf_path = gvcf_path.replace(s['id'], s['external_id'])
                if not utils.file_exists(gvcf_path):
                    continue
                logger.info(f'Using {gvcf_path} for a test run')

            if not gvcf_path.endswith('.g.vcf.gz'):
                logger.warning(
                    f'GVCF analysis for sample ID {sample_id} "output" field '
                    f'is not a GVCF'
                )
                continue
            if not utils.file_exists(gvcf_path):
                logger.warning(
                    f'GVCF analysis for sample ID {sample_id} "output" file '
                    f'does not exist: {gvcf_path}'
                )
                continue
            if not utils.file_exists(gvcf_path + '.tbi'):
                logger.warning(
                    f'GVCF analysis for sample ID {sample_id} "output" field '
                    f'does not have a corresponding tbi index: {gvcf_path}.tbi'
                )
                continue
            entry = default_entry.copy()
            entry.update(
                {
                    's': sample_id,
                    'external_id': external_id,
                    'project': proj,
                    'gvcf': gvcf_path,
                    'batch': seq_meta.get('batch', '-'),
                    'flowcell_lane': seq_meta.get('sample.flowcell_lane', '-'),
                    'library_id': seq_meta.get('sample.library_id', '-'),
                    'platform': seq_meta.get('sample.platform', '-'),
                    'centre': seq_meta.get('sample.centre', '-'),
                    'r_contamination': seq_meta.get('raw_data.FREEMIX'),
                    'r_chimera': seq_meta.get('raw_data.PCT_CHIMERAS'),
                    'r_duplication': seq_meta.get('raw_data.PERCENT_DUPLICATION'),
                    'median_insert_size': seq_meta.get('raw_data.MEDIAN_INSERT_SIZE'),
                }
            )
            inputs.append(entry)

    if not inputs:
        logger.error('No found any projects with samples good for processing')
        sys.exit(1)

    df = pd.DataFrame(inputs).set_index('s', drop=False)
    return df


def add_validation_samples(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add NA12878 GVCFs and syndip BAM into the dataframe.
    """
    if 'syndip' not in df.s:
        entry = default_entry.copy()
        entry.update(
            {
                's': 'syndip',
                'external_id': 'syndip',
                'project': 'syndip',
                'cram': 'gs://cpg-reference/validation/syndip/raw/CHM1_CHM13_2.bam',
                'crai': 'gs://cpg-reference/validation/syndip/raw/CHM1_CHM13_2.bam.bai',
            }
        )
        # Can only append a dict if ignore_index=True. So then need to set index back.
        df = df.append(entry, ignore_index=True).set_index('s', drop=False)

    giab_samples = ['NA12878', 'NA12891', 'NA12892']
    for sn in giab_samples:
        if sn not in df.s:
            cram = f'gs://cpg-reference/validation/giab/cram/{sn}.cram'
            entry = default_entry.copy()
            entry.update(
                {
                    's': sn,
                    'external_id': sn,
                    'project': 'giab',
                    'cram': cram,
                    'crai': cram + '.crai',
                }
            )
            # Can only append a dict if ignore_index=True. So then need to set index back.
            df = df.append(entry, ignore_index=True).set_index('s', drop=False)
    return df


@dataclass
class AlignmentInput:
    """
    Sort of a union type for possible alignment inputs
    """

    bam_or_cram_path: Optional[str] = None
    index_path: Optional[str] = None
    fqs1: Optional[List[str]] = None
    fqs2: Optional[List[str]] = None
