"""
Functions to find the pipeline inputs and communicate with the SM server
"""

import logging
import traceback
from dataclasses import dataclass
from textwrap import dedent
from typing import List, Dict, Optional, Set, Collection

from hailtop.batch import Batch
from hailtop.batch.job import Job
from sample_metadata import (
    AnalysisApi,
    SequenceApi,
    SampleApi,
    AnalysisUpdateModel,
    AnalysisModel,
    exceptions,
)

import utils


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@dataclass
class Analysis:
    """
    Represents the analysis SampleMetadata database entry
    """

    id: str
    type: str
    status: str
    sample_ids: Set[str]
    output: Optional[str]


class SMDB:
    """
    Singleton class abstracting the communication with
    the SampleMetadata database
    """

    sapi = SampleApi()
    aapi = AnalysisApi()
    seqapi = SequenceApi()
    do_update_analyses: bool = False

    @classmethod
    def get_samples_by_project(
        cls,
        projects: List[str],
        namespace: str,
        skip_samples: Optional[List[str]] = None,
    ) -> Dict[str, List[Dict]]:
        """
        Returns a dictionary of samples per input projects
        """
        samples_by_project: Dict[str, List[Dict]] = dict()
        for proj_name in projects:
            logger.info(f'Finding samples for project {proj_name}')
            input_proj = proj_name
            if namespace != 'main':
                input_proj += '-test'
            samples = cls.sapi.get_samples(
                body_get_samples_by_criteria_api_v1_sample_post={
                    'project_ids': [input_proj],
                    'active': True,
                }
            )
            samples_by_project[proj_name] = []
            for s in samples:
                if skip_samples and s['id'] in skip_samples:
                    logger.info(f'Skipping sample: {s["id"]}')
                    continue
                samples_by_project[proj_name].append(s)
        return samples_by_project

    @classmethod
    def find_seq_info_by_sid(cls, sample_ids) -> Dict[List, Dict]:
        """
        Return a dict of "Sequence" entries by sample ID
        """
        seq_infos: List[Dict] = cls.seqapi.get_sequences_by_sample_ids(sample_ids)
        seq_info_by_sid = dict(zip(sample_ids, seq_infos))
        return seq_info_by_sid

    @classmethod
    def update_analysis(
        cls,
        analysis: Analysis,
        status: Optional[str] = None,
        output: Optional[str] = None,
    ):
        """
        Update "status" of an Analysis entry
        """
        if not cls.do_update_analyses:
            return
        try:
            cls.aapi.update_analysis_status(
                analysis.id,
                AnalysisUpdateModel(
                    status=status or analysis.status,
                    output=output or analysis.output,
                ),
            )
        except exceptions.ApiException:
            traceback.print_exc()
        else:
            if status:
                analysis.status = status
            if output:
                analysis.output = output

    @classmethod
    def find_joint_calling_analysis(
        cls,
        analysis_project: str,
        sample_ids: Collection[str],
    ) -> Optional[Analysis]:
        """
        Query the DB to find the last completed joint-calling analysis for the samples
        """
        data = cls.aapi.get_latest_complete_analysis_for_type(
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

    @classmethod
    def find_analyses_by_sid(
        cls,
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
        try:
            logger.info(f'Querying analysis entries for project {analysis_project}')
            datas = cls.aapi.get_latest_analysis_for_samples_and_type(
                project=analysis_project,
                analysis_type=analysis_type,
                request_body=sample_ids,
            )
        except exceptions.ApiException:
            return dict()

        for data in datas:
            a = _parse_analysis(data)
            if not a:
                continue
            if a.status == 'completed':
                assert a.type == analysis_type, data
                assert len(a.sample_ids) == 1, data
                analysis_per_sid[list(a.sample_ids)[0]] = a
        return analysis_per_sid

    @classmethod
    def make_sm_in_progress_job(cls, *args, **kwargs) -> Job:
        """
        Creates a job that updates the sample metadata server entry analysis status
        to "in-progress"
        """
        kwargs['status'] = 'in-progress'
        return cls.make_sm_update_status_job(*args, **kwargs)

    @classmethod
    def make_sm_completed_job(cls, *args, **kwargs) -> Job:
        """
        Creates a job that updates the sample metadata server entry analysis status
        to "completed"
        """
        kwargs['status'] = 'completed'
        return cls.make_sm_update_status_job(*args, **kwargs)

    @classmethod
    def make_sm_update_status_job(
        cls,
        b: Batch,
        project: str,
        analysis_id: str,
        analysis_type: str,
        status: str,
        sample_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Job:
        """
        Creates a job that updates the sample metadata server entry analysis status.
        """
        assert status in ['in-progress', 'failed', 'completed', 'queued']
        job_name = ''
        if project_name and sample_name:
            job_name += f'{project_name}/{sample_name}: '
        job_name += f'Update SM: {analysis_type} to {status}'

        if not cls.do_update_analyses:
            return b.new_job(f'{job_name} [skip]')

        j = b.new_job(job_name)
        j.image(utils.SM_IMAGE)
        j.command(
            dedent(
                f"""
        set -o pipefail
        set -ex
        
        export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
        gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS
        export SM_DEV_DB_PROJECT={project}
        export SM_ENVIRONMENT=PRODUCTION
        
        cat <<EOT >> update.py
        from sample_metadata.api import AnalysisApi
        from sample_metadata import AnalysisUpdateModel
        from sample_metadata import exceptions
        import traceback
        aapi = AnalysisApi()
        try:
            aapi.update_analysis_status(
                analysis_id='{analysis_id}',
                analysis_update_model=AnalysisUpdateModel(status='{status}'),
            )
        except exceptions.ApiException:
            traceback.print_exc()
        EOT
        python update.py
        """
            )
        )
        return j

    @classmethod
    def create_analysis(
        cls,
        project: str,
        type_: str,
        output: str,
        status: str,
        sample_ids: Collection[str],
    ) -> Optional[int]:
        """
        Tries to create an Analysis entry, returns its id if successfuly
        """
        if not cls.do_update_analyses:
            return None

        am = AnalysisModel(
            type=type_,
            output=output,
            status=status,
            sample_ids=sample_ids,
        )
        aid = cls.aapi.create_new_analysis(project=project, analysis_model=am)
        logger.info(f'Created analysis of type={type_} status={status} with ID: {aid}')
        return aid

    @classmethod
    def process_existing_analysis(
        cls,
        sample_ids: Collection[str],
        completed_analysis: Optional[Analysis],
        analysis_type: str,
        analysis_project: str,
        analysis_sample_ids: Collection[str],
        expected_output_fpath: str,
        skip_stage: bool,
        check_existence: bool,
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
        :param check_existence: check if the files are on the bucket
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
                logger.info(
                    f'Updating analysis {completed_analysis.id}: {completed_analysis}'
                )
                SMDB.update_analysis(completed_analysis, output=expected_output_fpath)
                found_output_fpath = None
            elif check_existence and not utils.file_exists(found_output_fpath):
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
            elif utils.file_exists(expected_output_fpath):
                logger.info(f'Skipping stage, picking existing {expected_output_fpath}')
                return expected_output_fpath
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
            SMDB.update_analysis(completed_analysis, status='failed')

        # can reuse, need to create a completed one?
        if utils.file_exists(expected_output_fpath):
            logger.info(
                f'Output file {expected_output_fpath} already exists, so creating '
                f'an analysis {label} with status=completed'
            )
            cls.create_analysis(
                project=analysis_project,
                type_=analysis_type,
                output=expected_output_fpath,
                status='completed',
                sample_ids=analysis_sample_ids,
            )
            return expected_output_fpath

        # proceeding with the standard pipeline (creating status=queued, submitting jobs)
        else:
            logger.info(
                f'Expected output file {expected_output_fpath} does not exist, '
                f'so queueing analysis {label}'
            )
            return None


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
