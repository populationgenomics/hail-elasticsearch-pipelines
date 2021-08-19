#!/usr/bin/env python3

"""
Driver for loading data into SEQR for the CPG. See the README for more information.

- 2021/04/16 Michael Franklin and Vlad Savelyev
"""

import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from os.path import join, dirname, abspath, splitext, basename
from typing import Optional, List, Tuple, Set, Dict
from dataclasses import dataclass
import pandas as pd
import click
import hailtop.batch as hb
from analysis_runner import dataproc
from hailtop.batch.job import Job
from sample_metadata import (
    SampleApi,
    AnalysisApi,
    AnalysisModel,
    AnalysisUpdateModel,
    AnalysisType,
    AnalysisStatus,
)

from find_inputs import sm_verify_reads_data, AlignmentInput
from vqsr import make_vqsr_jobs
import utils


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


sapi = SampleApi()
aapi = AnalysisApi()


@click.command()
@click.option(
    '-n',
    '--namespace',
    'output_namespace',
    type=click.Choice(['main', 'test', 'tmp']),
    help='The bucket namespace to write the results to',
)
@click.option(
    '--analysis-project',
    'analysis_project',
    default='seqr',
    help='SM project name to write the intermediate/joint-calling analysis entries to',
)
@click.option(
    '--input-project',
    'input_projects',
    multiple=True,
    required=True,
    help='Only read samples that belong to the project(s). Can be set multiple times.',
)
@click.option('--output-version', 'output_version', type=str, required=True)
@click.option(
    '--output-projects',
    'output_projects',
    multiple=True,
    help='Only create ES indicies for the project(s). Can be set multiple times. '
    'Defaults to --input-projects. The name of the ES index will be suffixed '
    'with the dataset version (set by --version)',
)
@click.option(
    '--start-from-stage',
    'start_from_stage',
    type=click.Choice(['cram', 'gvcf', 'joint_calling', 'annotate', 'load_to_es']),
)
@click.option(
    '--skip-sample',
    '-S',
    'skip_samples',
    multiple=True,
    help='Don\'t process specified samples. Can be set multiple times.',
)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
@click.option(
    '--overwrite/--reuse',
    'overwrite',
    is_flag=True,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option(
    '--make-checkpoints',
    'make_checkpoints',
    is_flag=True,
    help='Create checkpoints for intermediate Hail data',
)
@click.option(
    '--skip-ped-checks',
    'skip_ped_checks',
    is_flag=True,
    help='Skip checking provided sex and pedigree against the inferred one',
)
@click.option('--vep-block-size', 'vep_block_size')
def main(
    output_namespace: str,
    analysis_project: str,
    input_projects: List[str],
    output_version: str,
    output_projects: Optional[List[str]],
    start_from_stage: str,
    skip_samples: List[str],
    keep_scratch: bool,
    overwrite: bool,
    dry_run: bool,
    make_checkpoints: bool,  # pylint: disable=unused-argument
    skip_ped_checks: bool,  # pylint: disable=unused-argument
    vep_block_size: Optional[int] = None,  # pylint: disable=unused-argument
):  # pylint: disable=missing-function-docstring
    billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'

    if output_namespace in ['test', 'tmp']:
        tmp_bucket_suffix = 'test-tmp'
    else:
        tmp_bucket_suffix = 'main-tmp'

    if output_namespace in ['test', 'main']:
        output_suffix = output_namespace
        web_bucket_suffix = f'{output_namespace}-web'
    else:
        output_suffix = 'test-tmp'
        web_bucket_suffix = 'test-tmp'
    out_bucket = (
        f'gs://cpg-seqr-{output_suffix}/seqr_loader/{analysis_project}/{output_version}'
    )
    tmp_bucket = f'gs://cpg-seqr-{tmp_bucket_suffix}/seqr_loader/{analysis_project}/{output_version}'
    web_bucket = f'gs://cpg-seqr-{web_bucket_suffix}/seqr_loader/{analysis_project}/{output_version}'

    assert input_projects
    if output_projects:
        if not all(op in input_projects for op in output_projects):
            logger.critical(
                'All output projects must be contained within '
                'the specified input projects'
            )

    hail_bucket = os.environ.get('HAIL_BUCKET')
    if not hail_bucket or keep_scratch:
        # Scratch files are large, so we want to use the temporary bucket to put them in
        hail_bucket = f'{tmp_bucket}/hail'
    logger.info(
        f'Starting hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
    )
    b = hb.Batch(
        f'Seqr loading. '
        f'Project: {analysis_project}, '
        f'input projects: {input_projects}, '
        f'dataset version: {output_version}, '
        f'namespace: "{output_namespace}"',
        backend=backend,
    )

    local_tmp_dir = tempfile.mkdtemp()

    b = _add_jobs(
        b=b,
        tmp_bucket=tmp_bucket,
        out_bucket=out_bucket,
        web_bucket=web_bucket,
        local_tmp_dir=local_tmp_dir,
        output_version=output_version,
        overwrite=overwrite,
        prod=output_namespace == 'main',
        input_projects=input_projects,
        output_projects=output_projects or input_projects,
        vep_block_size=vep_block_size,
        analysis_project=analysis_project,
        start_from_stage=start_from_stage,
        skip_samples=skip_samples,
    )
    if b:
        b.run(dry_run=dry_run, delete_scratch_on_exit=not keep_scratch, wait=False)
    shutil.rmtree(local_tmp_dir)


@dataclass
class Analysis:
    """
    Represents the analysis DB entry
    """

    id: str
    type: AnalysisType
    status: AnalysisStatus
    output: Optional[str]
    sample_ids: List[str]

    @staticmethod
    def from_db(**kwargs):
        """
        Convert from db keys, mainly converting id to id_
        """
        analysis_type = kwargs.pop('type', None)
        status = kwargs.pop('status', None)
        sample_ids = kwargs['sample_ids']
        output = kwargs.pop('output', [])
        return Analysis(
            id=kwargs.pop('id'),
            type=AnalysisType(analysis_type),
            status=AnalysisStatus(status),
            sample_ids=list(set(sorted(sample_ids))),
            output=output,
        )


def _get_latest_complete_analysis(
    analysis_project: str,
) -> Dict[Tuple[str, Tuple], Analysis]:
    """
    Returns a dictionary that maps a tuple (analysis type, sample ids) to the
    lastest complete analysis record (represented by a AnalysisModel object)
    """
    latest_by_type_and_sids = dict()
    for a_type in ['cram', 'gvcf', 'joint-calling']:
        for a_data in aapi.get_latest_complete_analyses_by_type(
            project=analysis_project,
            analysis_type=a_type,
        ):
            a: Analysis = Analysis.from_db(**a_data)
            latest_by_type_and_sids[(a_type, tuple(a.sample_ids))] = a
    return latest_by_type_and_sids


def _add_jobs(
    b: hb.Batch,
    out_bucket,
    tmp_bucket,
    web_bucket,  # pylint: disable=unused-argument
    local_tmp_dir,
    output_version: str,
    overwrite: bool,
    prod: bool,
    input_projects: List[str],
    output_projects: List[str],
    vep_block_size,
    analysis_project: str,
    start_from_stage: Optional[str],  # pylint: disable=unused-argument
    skip_samples: List[str],
) -> Optional[hb.Batch]:
    genomicsdb_bucket = f'{out_bucket}/genomicsdbs'
    # pylint: disable=unused-variable
    fingerprints_bucket = f'{out_bucket}/fingerprints'

    latest_by_type_and_sids = _get_latest_complete_analysis(analysis_project)
    reference, bwa_reference, noalt_regions = utils.get_refs(b)
    intervals_j = _add_split_intervals_job(
        b=b,
        interval_list=utils.UNPADDED_INTERVALS,
        scatter_count=utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
        ref_fasta=utils.REF_FASTA,
    )

    gvcf_jobs = []
    gvcf_by_sid: Dict[str, str] = dict()

    samples_by_project: Dict[str, List[Dict]] = dict()
    for proj in input_projects:
        samples = sapi.get_samples(
            body_get_samples_by_criteria_api_v1_sample_post={
                'project_ids': [proj],
                'active': True,
            }
        )
        samples_by_project[proj] = []
        for s in samples:
            if s['id'] in skip_samples:
                logger.info(f'Skiping sample: {s["id"]}')
            else:
                samples_by_project[proj].append(s)

    # after dropping samples with incorrect metadata, missing inputs, etc
    good_samples = []
    for proj, samples in samples_by_project.items():
        for s in samples:
            logger.info(f'Project {proj}. Processing sample {s["id"]}')
            cram_analysis = latest_by_type_and_sids.get(('cram', (s['id'],)))

            # Requested to skip this stage, but an existing analysis is not found
            if start_from_stage is not None and start_from_stage not in ['cram']:
                cram_job = None
                if not cram_analysis:
                    logger.warning(
                        f'No completed CRAM analysis found for sample {s["id"]}, '
                        f'and start_from_stage is {start_from_stage}, so'
                        f'skipping this sample'
                    )
                    continue
                cram_fpath = cram_analysis.output
                if not cram_fpath or not utils.file_exists(cram_fpath):
                    logger.error(
                        f'Sample {s["id"]} has a completed CRAM analysis, '
                        f'but the output file {cram_fpath} does not exist, so '
                        f'skipping this sample'
                    )
                    continue
                cram_fpath = str(cram_fpath)
            else:
                alignment_input = sm_verify_reads_data(
                    s['meta'].get('reads'), s['meta'].get('reads_type')
                )
                if not alignment_input:
                    continue
                cram_job, cram_fpath = _make_realign_jobs(
                    b=b,
                    sample_name=s['id'],
                    project_name=proj,
                    alignment_input=alignment_input,
                    reference=bwa_reference,
                    out_bucket=out_bucket,
                    overwrite=overwrite,
                    completed_analysis=cram_analysis,
                    analysis_project=analysis_project,
                )

            gvcf_analysis = latest_by_type_and_sids.get(('gvcf', (s['id'],)))
            if start_from_stage is not None and start_from_stage not in [
                'cram',
                'gvcf',
            ]:
                if not gvcf_analysis:
                    logger.warning(
                        f'No completed GVCF analysis found for sample {s["id"]}, '
                        f'and start_from_stage is {start_from_stage}, so'
                        f'skipping this sample'
                    )
                    continue
                gvcf_fpath = gvcf_analysis.output
                if not gvcf_fpath or not utils.file_exists(gvcf_fpath):
                    logger.error(
                        f'Sample {s["id"]} has a completed GVCF analysis, '
                        f'but the output file {gvcf_fpath} does not exist, so '
                        f'skipping this sample'
                    )
                    continue
                gvcf_fpath = str(gvcf_fpath)
            else:
                gvcf_job, gvcf_fpath = _make_produce_gvcf_jobs(
                    b=b,
                    sample_name=s['id'],
                    project_name=proj,
                    cram_path=cram_fpath,
                    intervals_j=intervals_j,
                    reference=reference,
                    noalt_regions=noalt_regions,
                    out_bucket=out_bucket,
                    tmp_bucket=tmp_bucket,
                    overwrite=overwrite,
                    depends_on=[cram_job] if cram_job else [],
                    analysis_project=analysis_project,
                    completed_analysis=gvcf_analysis,
                )
                gvcf_jobs.append(gvcf_job)

            gvcf_by_sid[s['id']] = gvcf_fpath
            good_samples.append(s)

    if not good_samples:
        logger.info('No samples left to joint-call')
        return None

    # Is there a complete joint-calling analysis for the requested set of samples?
    jc_analysis = latest_by_type_and_sids.get(
        ('joint-calling', tuple(set(s['id'] for s in good_samples)))
    )
    if start_from_stage is None and start_from_stage not in [
        'cram',
        'gvcf',
        'joint_calling',
    ]:
        if not jc_analysis:
            logger.warning(
                f'No joint-caling analysis found, '
                f'and start_from_stage is {start_from_stage}, so '
                f'stopping here'
            )
            return None

    jc_job, jc_vcf_path = _make_joint_genotype_jobs(
        b=b,
        genomicsdb_bucket=genomicsdb_bucket,
        samples=good_samples,
        gvcf_by_sid=gvcf_by_sid,
        reference=reference,
        dbsnp=utils.DBSNP_VCF,
        out_bucket=out_bucket,
        tmp_bucket=out_bucket,
        local_tmp_dir=local_tmp_dir,
        overwrite=overwrite,
        analysis_project=analysis_project,
        completed_analysis=jc_analysis,
        depends_on=gvcf_jobs,
    )

    for project in output_projects:
        annotated_mt_path = join(out_bucket, 'annotation', f'{project}.mt')
        if utils.can_reuse(annotated_mt_path, overwrite):
            annotate_job = b.new_job(f'{project}: annotate [reuse]')
        else:
            annotate_job = dataproc.hail_dataproc_job(
                b,
                f'batch_seqr_loader/scripts/make_annotated_mt.py '
                f'--source-path {jc_vcf_path} '
                f'--dest-mt-path {annotated_mt_path} '
                f'--bucket {join(tmp_bucket, "annotation", project)} '
                '--disable-validation '
                '--make-checkpoints '
                + (f'--vep-block-size ' if vep_block_size else ''),
                max_age='16h',
                packages=utils.DATAPROC_PACKAGES,
                num_secondary_workers=utils.NUMBER_OF_DATAPROC_WORKERS,
                job_name=f'Annotate {project}',
                vep='GRCh38',
                depends_on=[jc_job],
            )

        dataproc.hail_dataproc_job(
            b,
            f'batch_seqr_loader/scripts/load_to_es.py '
            f'--mt-path {annotated_mt_path} '
            f'--es-index {project}-{output_version} '
            f'--es-index-min-num-shards 1 '
            f'--genome-version GRCh38 '
            f'{"--prod" if prod else ""}',
            max_age='16h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            job_name=f'{project}: add to the ES index',
            depends_on=[annotate_job],
            scopes=['cloud-platform'],
        )
    return b


def _pedigree_checks(
    b: hb.Batch,
    samples_df: pd.DataFrame,
    reference: hb.ResourceGroup,
    sites: hb.ResourceFile,
    ped_file: hb.ResourceFile,
    overwrite: bool,  # pylint: disable=unused-argument
    fingerprints_bucket: str,
    web_bucket: str,
    web_url: str,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[Job, str]:
    """
    Add somalier and peddy based jobs that infer relatedness and sex, compare that
    to the provided PED file, and attempt to recover it. If unable to recover, cancel
    the further workflow jobs.

    Returns a job, and a bucket path to a fixed PED file if able to recover.
    """

    extract_jobs = []
    fp_file_by_sample = dict()
    for sn, input_path, input_index in zip(
        samples_df['s'], samples_df['file'], samples_df['index']
    ):
        fp_file_by_sample[sn] = join(fingerprints_bucket, f'{sn}.somalier')
        if utils.can_reuse(fp_file_by_sample[sn], overwrite):
            extract_jobs.append(b.new_job(f'Somalier extract, {sn} [reuse]'))
        else:
            j = b.new_job(f'Somalier extract, {sn}')
            j.image(utils.SOMALIER_IMAGE)
            j.memory('standard')
            if input_path.endswith('.bam'):
                j.cpu(4)
                j.storage(f'200G')
            elif input_path.endswith('.cram'):
                j.cpu(4)
                j.storage(f'50G')
            else:
                j.cpu(2)
                j.storage(f'10G')
            if depends_on:
                j.depends_on(*depends_on)

            input_file = b.read_input_group(
                base=input_path,
                index=input_index,
            )

            j.command(
                f"""set -ex
                
                somalier extract -d extracted/ --sites {sites} -f {reference.base} \\
                {input_file['base']}
                
                mv extracted/*.somalier {j.output_file}
                """
            )
            b.write_output(j.output_file, fp_file_by_sample[sn])
            extract_jobs.append(j)

    relate_j = b.new_job(f'Somalier relate')
    relate_j.image(utils.SOMALIER_IMAGE)
    relate_j.cpu(1)
    relate_j.memory('standard')  # ~ 4G/core ~ 4G
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    relate_j.storage(f'{1 + len(extract_jobs) // 4000 * 1}G')
    relate_j.depends_on(*extract_jobs)
    fp_files = [b.read_input(fp) for sn, fp in fp_file_by_sample.items()]
    relate_j.command(
        f"""set -e

        cat {ped_file} | grep -v Family.ID > samples.ped 

        somalier relate \\
        {' '.join(fp_files)} \\
        --ped samples.ped \\
        -o related \\
        --infer

        ls
        mv related.html {relate_j.output_html}
        mv related.pairs.tsv {relate_j.output_pairs}
        mv related.samples.tsv {relate_j.output_samples}
      """
    )

    # Copy somalier outputs to buckets
    sample_hash = utils.hash_sample_names(samples_df['s'])
    prefix = join(fingerprints_bucket, sample_hash, 'somalier')
    somalier_samples_path = f'{prefix}.samples.tsv'
    somalier_pairs_path = f'{prefix}.pairs.tsv'
    b.write_output(relate_j.output_samples, somalier_samples_path)
    b.write_output(relate_j.output_pairs, somalier_pairs_path)
    # Copy somalier HTML to the web bucket
    rel_path = join('loader', sample_hash[:10], 'somalier.html')
    somalier_html_path = join(web_bucket, rel_path)
    somalier_html_url = f'{web_url}/{rel_path}'
    b.write_output(relate_j.output_html, somalier_html_path)

    check_j = b.new_job(f'Check relatedness and sex')
    check_j.image(utils.PEDDY_IMAGE)
    check_j.cpu(1)
    check_j.memory('standard')  # ~ 4G/core ~ 4G
    with open(join(dirname(abspath(__file__)), 'check_pedigree.py')) as f:
        script = f.read()
    check_j.command(
        f"""set -e
cat <<EOT >> check_pedigree.py
{script}
EOT
python check_pedigree.py \
--somalier-samples {relate_j.output_samples} \
--somalier-pairs {relate_j.output_pairs} \
{('--somalier-html ' + somalier_html_url) if somalier_html_url else ''}
    """
    )

    check_j.depends_on(relate_j)
    return check_j, somalier_samples_path


def _make_realign_jobs(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    alignment_input: AlignmentInput,
    reference: hb.ResourceGroup,
    out_bucket: str,
    overwrite: bool,
    depends_on: Optional[List[Job]] = None,
    completed_analysis: Optional[Analysis] = None,
    analysis_project: Optional[str] = None,
) -> Tuple[Job, str]:
    """
    Runs BWA to realign reads back again to hg38.

    When the input is CRAM/BAM, uses Bazam to stream reads to BWA.
    """
    job_name = f'{project_name}/{sample_name}: BWA align'
    output_cram_path = join(out_bucket, f'{sample_name}.cram')
    if completed_analysis:
        if utils.file_exists(output_cram_path):
            j = b.new_job(f'{job_name} [reuse]')
            return j, output_cram_path
        else:
            logger.warning(
                f'Sample {sample_name} does have a completed CRAM analysis, '
                f'but the output file does not exist. Setting status to "failed" '
                f'and rerunning the analysis.'
            )
            if analysis_project:
                aapi.update_analysis_status(
                    completed_analysis.id, AnalysisUpdateModel(status='failed')
                )

    logger.info(
        f'Sample {sample_name} does not have a completed CRAM analysis yet. '
        f'Parsing the "reads" metadata field and submitting the alignmentment'
    )
    am = AnalysisModel(
        type='cram',
        output=output_cram_path,
        status='queued',
        sample_ids=[sample_name],
    )
    if utils.can_reuse(output_cram_path, overwrite):
        # Create a "completed" analysis and return an empty job
        am.status = 'completed'
        if analysis_project:
            aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        return b.new_job(f'{job_name} [reuse]'), output_cram_path

    j = b.new_job(job_name)
    if analysis_project:
        # Interacting with the sample metadata server:
        # 1. Create a "queued" analysis
        aid = aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        # 2. Queue a job that updates the status to "in-progress"
        sm_in_progress_j = utils.make_sm_in_progress_job(
            b,
            project=analysis_project,
            analysis_id=aid,
            analysis_type='cram',
            project_name=project_name,
            sample_name=sample_name,
        )
        # 2. Queue a job that updates the status to "completed"
        sm_completed_j = utils.make_sm_completed_job(
            b,
            project=analysis_project,
            analysis_id=aid,
            analysis_type='cram',
            project_name=project_name,
            sample_name=sample_name,
        )
        # Set up dependencies
        j.depends_on(sm_in_progress_j)
        sm_completed_j.depends_on(j)
        if depends_on:
            sm_in_progress_j.depends_on(*depends_on)
        logger.info(f'Queueing CRAM re-alignment analysis')
    else:
        sm_completed_j = None
        if depends_on:
            j.depends_on(*depends_on)

    j.image(utils.BAZAM_IMAGE)
    total_cpu = 32

    if alignment_input.bam_or_cram_path:
        use_bazam = True
        bazam_cpu = 10
        bwa_cpu = 32
        bamsormadup_cpu = 10

        assert alignment_input.index_path
        assert not alignment_input.fqs1 and not alignment_input.fqs2
        cram = b.read_input_group(
            base=alignment_input.bam_or_cram_path, index=alignment_input.index_path
        )
        r1_param = (
            f'<(bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base}'
            f' -n{bazam_cpu} -bam {cram.base})'
        )
        r2_param = '-'
    else:
        assert alignment_input.fqs1 and alignment_input.fqs2
        use_bazam = False
        bwa_cpu = 32
        bamsormadup_cpu = 10
        files1 = [b.read_input(f1) for f1 in alignment_input.fqs1]
        files2 = [b.read_input(f1) for f1 in alignment_input.fqs2]
        r1_param = f'<(cat {" ".join(files1)})'
        r2_param = f'<(cat {" ".join(files2)})'
        logger.info(f'r1_param: {r1_param}')
        logger.info(f'r2_param: {r2_param}')

    j.cpu(total_cpu)
    j.memory('standard')
    j.storage('300G')
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'crai': '{root}.crai',
        }
    )

    rg_line = f'@RG\\tID:{sample_name}\\tSM:{sample_name}'
    # BWA command options:
    # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
    # -p     smart pairing (ignoring in2.fq)
    # -t16   threads
    # -Y     use soft clipping for supplementary alignments
    # -R     read group header line such as '@RG\tID:foo\tSM:bar'
    command = f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {j.output_cram.cram}); sleep 600; done) &

bwa mem -K 100000000 {'-p' if use_bazam else ''} -t{bwa_cpu} -Y \\
    -R '{rg_line}' {reference.base} {r1_param} {r2_param} | \\
bamsormadup inputformat=sam threads={bamsormadup_cpu} SO=coordinate \\
    M={j.duplicate_metrics} outputformat=sam \\
    tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp | \\
samtools view -T {reference.base} -O cram -o {j.output_cram.cram}

samtools index -@{total_cpu} {j.output_cram.cram} {j.output_cram.crai}

df -h; pwd; du -sh $(dirname {j.output_cram.cram})
    """
    j.command(command)
    b.write_output(j.output_cram, splitext(output_cram_path)[0])
    b.write_output(
        j.duplicate_metrics,
        join(dirname(output_cram_path), f'{sample_name}-duplicate-metrics.csv'),
    )
    if sm_completed_j:
        return sm_completed_j, output_cram_path
    else:
        return j, output_cram_path


def _make_produce_gvcf_jobs(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    cram_path: str,
    intervals_j: Job,
    reference: hb.ResourceGroup,
    noalt_regions: hb.ResourceFile,
    out_bucket: str,
    tmp_bucket: str,  # pylint: disable=unused-argument
    overwrite: bool,  # pylint: disable=unused-argument
    depends_on: Optional[List[Job]] = None,
    analysis_project: str = None,
    completed_analysis: Optional[Analysis] = None,
) -> Tuple[Job, str]:
    """
    Takes all samples with a 'file' of 'type'='bam' in `samples_df`,
    and runs HaplotypeCaller on them, and sets a new 'file' of 'type'='gvcf'

    HaplotypeCaller is run in an interval-based sharded way, with per-interval
    HaplotypeCaller jobs defined in a nested loop.
    """
    job_name = f'{project_name}/{sample_name}: make GVCF'
    out_gvcf_path = join(out_bucket, f'{sample_name}.g.vcf.gz')
    if completed_analysis:
        if utils.file_exists(out_gvcf_path):
            j = b.new_job(f'{job_name} [reuse]')
            logger.info(
                f'Completed analysis and the output for a "{job_name}" '
                f'job exists, reusing: {out_gvcf_path}'
            )
            return j, out_gvcf_path
        else:
            logger.warning(
                f'Sample {sample_name} does have a completed GVCF analysis, '
                f'but the output file {out_gvcf_path} does not exist. Setting status '
                f'to "failed" and rerunning the analysis.'
            )
            aapi.update_analysis_status(
                completed_analysis.id, AnalysisUpdateModel(status='failed')
            )
    else:
        logger.info(
            f'Sample {sample_name} does not have a completed GVCF analysis yet. '
            f'Submitting the variant calling jobs.'
        )

    am = AnalysisModel(
        type='gvcf',
        output=out_gvcf_path,
        status='queued',
        sample_ids=[sample_name],
    )
    if utils.can_reuse(out_gvcf_path, overwrite=overwrite):
        # Create a "completed" analysis and return an empty job
        am.status = 'completed'
        if analysis_project:
            aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        return b.new_job(f'{job_name} [reuse]'), out_gvcf_path

    haplotype_caller_jobs = []
    for idx in range(utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS):
        haplotype_caller_jobs.append(
            _add_haplotype_caller_job(
                b,
                sample_name=sample_name,
                project_name=project_name,
                cram=b.read_input_group(
                    **{
                        'cram': cram_path,
                        'crai': re.sub('.cram$', '.crai', cram_path),
                    }
                ),
                interval=intervals_j.intervals[f'interval_{idx}'],
                reference=reference,
                interval_idx=idx,
                number_of_intervals=utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
                depends_on=depends_on,
            )
        )
    hc_gvcf_path = join(tmp_bucket, 'raw', f'{sample_name}.g.vcf.gz')
    merge_j = _add_merge_gvcfs_job(
        b=b,
        sample_name=sample_name,
        project_name=project_name,
        gvcfs=[j.output_gvcf for j in haplotype_caller_jobs],
        output_gvcf_path=hc_gvcf_path,
    )

    postproc_job = _make_postproc_gvcf_jobs(
        b=b,
        sample_name=sample_name,
        project_name=project_name,
        input_gvcf_path=hc_gvcf_path,
        out_gvcf_path=out_gvcf_path,
        noalt_regions=noalt_regions,
        overwrite=overwrite,
        depends_on=[merge_j],
    )

    if analysis_project:
        # Interacting with the sample metadata server:
        # 1. Create a "queued" analysis
        aid = aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        # 2. Queue a job that updates the status to "in-progress"
        sm_in_progress_j = utils.make_sm_in_progress_job(
            b,
            project=analysis_project,
            analysis_id=aid,
            analysis_type='gvcf',
            project_name=project_name,
            sample_name=sample_name,
        )
        # 2. Queue a job that updates the status to "completed"
        sm_completed_j = utils.make_sm_completed_job(
            b,
            project=analysis_project,
            analysis_id=aid,
            analysis_type='gvcf',
            project_name=project_name,
            sample_name=sample_name,
        )
        # Set up dependencies
        haplotype_caller_jobs[0].depends_on(sm_in_progress_j)
        if depends_on:
            sm_in_progress_j.depends_on(*depends_on)
        logger.info(f'Queueing GVCF analysis')
    else:
        if depends_on:
            haplotype_caller_jobs[0].depends_on(*depends_on)
        sm_completed_j = None

    if sm_completed_j:
        sm_completed_j.depends_on(postproc_job)
        return sm_completed_j, out_gvcf_path
    else:
        return postproc_job, out_gvcf_path


def _make_postproc_gvcf_jobs(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    input_gvcf_path: str,
    out_gvcf_path: str,
    noalt_regions: hb.ResourceFile,
    overwrite: bool,  # pylint: disable=unused-argument
    depends_on: Optional[List[Job]] = None,
) -> Job:
    reblock_gvcf_job = _add_reblock_gvcf_job(
        b,
        sample_name=sample_name,
        project_name=project_name,
        input_gvcf=b.read_input_group(
            **{'g.vcf.gz': input_gvcf_path, 'g.vcf.gz.tbi': input_gvcf_path + '.tbi'}
        ),
        overwrite=overwrite,
    )
    if depends_on:
        reblock_gvcf_job.depends_on(*depends_on)
    subset_to_noalt_job = _add_subset_noalt_step(
        b,
        sample_name=sample_name,
        project_name=project_name,
        input_gvcf=reblock_gvcf_job.output_gvcf,
        noalt_regions=noalt_regions,
        overwrite=overwrite,
        output_gvcf_path=out_gvcf_path,
    )
    return subset_to_noalt_job


def _add_reblock_gvcf_job(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    input_gvcf: hb.ResourceGroup,
    overwrite: bool,
    output_gvcf_path: Optional[str] = None,
) -> Job:
    """
    Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration
    """
    job_name = f'{project_name}/{sample_name}: ReblockGVCF'
    if utils.can_reuse(output_gvcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'30G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )

    j.command(
        f"""
    gatk --java-options "-Xms{mem_gb - 1}g" \\
        ReblockGVCF \\
        -V {input_gvcf['g.vcf.gz']} \\
        -do-qual-approx \\
        -O {j.output_gvcf['g.vcf.gz']} \\
        --create-output-variant-index true"""
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def _add_subset_noalt_step(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    input_gvcf: hb.ResourceGroup,
    noalt_regions: str,
    overwrite: bool,
    output_gvcf_path: Optional[str] = None,
) -> Job:
    """
    1. Subset GVCF to main chromosomes to avoid downstream errors
    2. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    """
    job_name = f'{project_name}/{sample_name}: SubsetToNoalt'
    if utils.can_reuse(output_gvcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.BCFTOOLS_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'30G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )
    j.command(
        f"""set -e

    bcftools view \\
        {input_gvcf['g.vcf.gz']} \\
        -T {noalt_regions} \\
        | bcftools annotate -x INFO/DS \\
        -o {j.output_gvcf['g.vcf.gz']} \\
        -Oz

    bcftools index --tbi {j.output_gvcf['g.vcf.gz']}
        """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def _make_joint_genotype_jobs(
    b: hb.Batch,
    genomicsdb_bucket: str,
    samples: List,
    gvcf_by_sid: Dict[str, str],
    reference: hb.ResourceGroup,
    dbsnp: str,
    out_bucket: str,
    tmp_bucket: str,
    local_tmp_dir: str,
    overwrite: bool,
    depends_on: Optional[List[Job]] = None,
    analysis_project: str = None,
    completed_analysis: Optional[Analysis] = None,
) -> Tuple[Job, str]:
    """
    Assumes all samples have a 'file' of 'type'='gvcf' in `samples_df`.
    Adds samples to the GenomicsDB and runs joint genotyping on them.
    Outputs a multi-sample VCF under `output_vcf_path`.
    """
    is_small_callset = len(samples) < 1000
    # 1. For small callsets, we don't apply the ExcessHet filtering.
    # 2. For small callsets, we gather the VCF shards and collect QC metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics
    # collected from them.
    is_huge_callset = len(samples) >= 100000
    # For huge callsets, we allocate more memory for the SNPs Create Model step

    genomics_gcs_path_per_interval = dict()
    for idx in range(utils.NUMBER_OF_GENOMICS_DB_INTERVALS):
        genomics_gcs_path_per_interval[idx] = join(
            genomicsdb_bucket,
            f'interval_{idx}_outof_{utils.NUMBER_OF_GENOMICS_DB_INTERVALS}',
        )
    # Determining which samples to add. Using the first interval, so the assumption
    # is that all DBs have the same set of samples.
    (
        sample_names_to_add,
        sample_names_will_be_in_db,
        sample_names_already_added,
        updating_existing_db,
        sample_map_bucket_path,
    ) = _samples_to_add_to_db(
        genomicsdb_gcs_path=genomics_gcs_path_per_interval[0],
        interval_idx=0,
        local_tmp_dir=local_tmp_dir,
        samples=samples,
        tmp_bucket=tmp_bucket,
        gvcf_by_sid=gvcf_by_sid,
    )
    if not sample_names_to_add:
        logger.info('')

    samples_hash = utils.hash_sample_names(sample_names_will_be_in_db)
    gathered_vcf_path = join(
        out_bucket, 'jointly-called', 'tmp', 'gathered', samples_hash + '.vcf.gz'
    )
    job_name = 'Joint-calling+VQSR'
    vqsred_vcf_path = join(out_bucket, 'jointly-called', samples_hash + '.vcf.gz')
    if completed_analysis:
        if utils.file_exists(vqsred_vcf_path):
            logger.info(
                f'Completed analysis and the output for these samples exist: '
                f'"{vqsred_vcf_path}", reusing.'
            )
            j = b.new_job('Joint calling+VQSR [reuse]')
            return j, vqsred_vcf_path
        else:
            logger.warning(
                f'Joint-calling completed analysis for these samples exists, '
                f'but the output file is missing. Setting status to "failed" and '
                f'rerunning the joint-calling analysis.'
            )
            aapi.update_analysis_status(
                completed_analysis.id,
                AnalysisUpdateModel(status='failed'),
            )
    else:
        logger.info(
            f'Completed joint-calling analysis does not exist for this set of samples. '
            f'Submitting the joint-calling and VQSR jobs.'
        )

    am = AnalysisModel(
        sample_ids=[s['id'] for s in samples],
        type='joint-calling',
        output=vqsred_vcf_path,
        status='queued',
    )
    if utils.can_reuse(vqsred_vcf_path, overwrite=overwrite):
        # Create a "completed" analysis and return an empty job
        am.status = 'completed'
        if analysis_project:
            aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        return b.new_job(f'{job_name} [reuse]'), vqsred_vcf_path

    intervals_j = _add_split_intervals_job(
        b=b,
        interval_list=utils.UNPADDED_INTERVALS,
        scatter_count=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
        ref_fasta=utils.REF_FASTA,
    )

    if analysis_project:
        # Interacting with the sample metadata server:
        # 1. Create a "queued" analysis
        aid = aapi.create_new_analysis(project=analysis_project, analysis_model=am)
        # 2. Queue a job that updates the status to "in-progress"
        sm_in_progress_j = utils.make_sm_in_progress_job(
            b,
            project=analysis_project,
            analysis_id=aid,
            analysis_type='joint-calling',
        )
        # 2. Queue a job that updates the status to "completed"
        sm_completed_j = utils.make_sm_completed_job(
            b,
            project=analysis_project,
            analysis_id=aid,
            analysis_type='joint-calling',
        )
        # Set up dependencies
        intervals_j.depends_on(sm_in_progress_j)
        if depends_on:
            sm_in_progress_j.depends_on(*depends_on)
        logger.info(f'Queueing {am.type} with analysis ID: {aid}')
    else:
        if depends_on:
            intervals_j.depends_on(*depends_on)
        sm_completed_j = None

    import_gvcfs_job_per_interval = dict()
    if sample_names_to_add:
        for idx in range(utils.NUMBER_OF_GENOMICS_DB_INTERVALS):
            import_gvcfs_job, _ = _add_import_gvcfs_job(
                b=b,
                genomicsdb_gcs_path=genomics_gcs_path_per_interval[idx],
                sample_names_to_add=sample_names_to_add,
                sample_names_to_skip=sample_names_already_added,
                sample_names_will_be_in_db=sample_names_will_be_in_db,
                updating_existing_db=updating_existing_db,
                sample_map_bucket_path=sample_map_bucket_path,
                interval=intervals_j.intervals[f'interval_{idx}'],
                interval_idx=idx,
                number_of_intervals=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
                depends_on=[intervals_j],
            )
            import_gvcfs_job_per_interval[idx] = import_gvcfs_job

    make_site_only_jobs = []
    scattered_vcf_by_interval: Dict[int, hb.ResourceGroup] = dict()

    for idx in range(utils.NUMBER_OF_GENOMICS_DB_INTERVALS):
        samples_hash = utils.hash_sample_names(sample_names_will_be_in_db)
        site_only_vcf_path = join(
            tmp_bucket, 'jointly-called', samples_hash, f'interval_{idx}.vcf.gz'
        )
        if utils.can_reuse(site_only_vcf_path, overwrite):
            make_site_only_jobs.append(b.new_job('Joint genotyping [reuse]'))
            scattered_vcf_by_interval[idx] = b.read_input_group(
                **{
                    'vcf.gz': site_only_vcf_path,
                    'vcf.gz.tbi': site_only_vcf_path + '.tbi',
                }
            )
        else:
            genotype_vcf_job = _add_genotype_gvcfs_job(
                b,
                genomicsdb_path=genomics_gcs_path_per_interval[idx],
                interval=intervals_j.intervals[f'interval_{idx}'],
                reference=reference,
                dbsnp=dbsnp,
                overwrite=overwrite,
                interval_idx=idx,
                number_of_samples=len(sample_names_will_be_in_db),
                number_of_intervals=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
            )
            if import_gvcfs_job_per_interval.get(idx):
                genotype_vcf_job.depends_on(import_gvcfs_job_per_interval.get(idx))

            vcf = genotype_vcf_job.output_vcf
            if not is_small_callset:
                exccess_filter_job = _add_exccess_het_filter(
                    b,
                    input_vcf=vcf,
                    overwrite=overwrite,
                    interval=intervals_j.intervals[f'interval_{idx}'],
                )
                vcf = exccess_filter_job.output_vcf
            make_site_only_job = _add_make_sites_only_job(
                b,
                input_vcf=vcf,
                overwrite=overwrite,
            )
            make_site_only_jobs.append(make_site_only_job)
            scattered_vcf_by_interval[idx] = make_site_only_job.output_vcf

    final_gathered_vcf_job = _add_final_gather_vcf_job(
        b,
        input_vcfs=list(scattered_vcf_by_interval.values()),
        overwrite=overwrite,
        output_vcf_path=gathered_vcf_path,
    )
    tmp_vqsr_bucket = join(tmp_bucket, 'vqsr', samples_hash)
    vqsr_job = make_vqsr_jobs(
        b,
        input_vcf_gathered=gathered_vcf_path,
        input_vcfs_scattered=list(scattered_vcf_by_interval.values()),
        is_small_callset=is_small_callset,
        is_huge_callset=is_huge_callset,
        work_bucket=tmp_vqsr_bucket,
        web_bucket=tmp_vqsr_bucket,
        depends_on=[final_gathered_vcf_job],
        intervals=intervals_j.intervals,
        scatter_count=utils.NUMBER_OF_GENOMICS_DB_INTERVALS,
        output_vcf_path=vqsred_vcf_path,
    )
    if sm_completed_j:
        sm_completed_j.depends_on(vqsr_job)
        return sm_completed_j, vqsred_vcf_path
    else:
        return vqsr_job, vqsred_vcf_path


def _add_split_intervals_job(
    b: hb.Batch,
    interval_list: str,
    scatter_count: int,
    ref_fasta: str,
) -> Job:
    """
    Split genome into intervals to parallelise GnarlyGenotyper.

    Returns: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f'Make {scatter_count} intervals')
    j.image(utils.GATK_IMAGE)
    java_mem = 3
    j.memory('standard')  # ~ 4G/core ~ 4G
    j.storage('16G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )

    j.command(
        f"""set -e

    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number 
    # of intervals. But we have to expect exactly the {scatter_count} number of 
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{java_mem}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta} \\
      -mode INTERVAL_SUBDIVISION
      """
    )
    # Could save intervals to a bucket here to avoid rerunning the job
    return j


def _add_haplotype_caller_job(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    cram: hb.ResourceGroup,
    interval: hb.ResourceFile,
    reference: hb.ResourceGroup,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
    output_gvcf_path: Optional[str] = None,
    overwrite: bool = False,
) -> Job:
    """
    Run HaplotypeCaller on an input BAM or CRAM, and output GVCF
    """
    job_name = f'{project_name}/{sample_name}: HaplotypeCaller'
    if interval_idx is not None:
        job_name += f', {sample_name} {interval_idx}/{number_of_intervals}'
    if utils.can_reuse(output_gvcf_path, overwrite):
        return b.new_job(f'{job_name} [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage('60G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )
    if depends_on:
        j.depends_on(*depends_on)

    j.command(
        f"""set -e
    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m; sleep 300; done) &

    gatk --java-options "-Xms{java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
      HaplotypeCaller \\
      -R {reference.base} \\
      -I {cram['cram']} \\
      -L {interval} \\
      -O {j.output_gvcf['g.vcf.gz']} \\
      -G AS_StandardAnnotation \\
      -GQB 20 \
      -ERC GVCF \\

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m
    """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path)
    return j


def _add_merge_gvcfs_job(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    gvcfs: List[hb.ResourceGroup],
    output_gvcf_path: Optional[str],
) -> Job:
    """
    Combine by-interval GVCFs into a single sample GVCF file
    """

    job_name = f'{project_name}/{sample_name}: merge {len(gvcfs)} GVCFs'
    j = b.new_job(job_name)
    j.image(utils.PICARD_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{len(gvcfs) * 1.5 + 2}G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    input_cmd = ' '.join(f'INPUT={g["g.vcf.gz"]}' for g in gvcfs)

    j.command(
        f"""set -e

    (while true; do df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m; sleep 300; done) &

    java -Xms{java_mem}g -jar /usr/picard/picard.jar \
      MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}

    df -h; pwd; du -sh $(dirname {j.output_gvcf['g.vcf.gz']}); free -m
      """
    )
    if output_gvcf_path:
        b.write_output(j.output_gvcf, output_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def _samples_to_add_to_db(
    genomicsdb_gcs_path,
    interval_idx,
    local_tmp_dir,
    samples,
    tmp_bucket: str,
    gvcf_by_sid: Dict[str, str],
) -> Tuple[Set[str], Set[str], Set[str], bool, str]:
    if utils.file_exists(join(genomicsdb_gcs_path, 'callset.json')):
        # Check if sample exists in the DB already
        genomicsdb_metadata = join(local_tmp_dir, f'callset-{interval_idx}.json')
        # This command will download the DB metadata file locally.
        # The `-O` argument to `tar` means "write the file being extracted to the stdout",
        # and the file to be extracted is specified as a positional argument to `tar`.
        cmd = f'gsutil cp {join(genomicsdb_gcs_path, "callset.json")} {genomicsdb_metadata}'
        logger.info(cmd)
        subprocess.run(cmd, check=False, shell=True)

        with open(genomicsdb_metadata) as f:
            db_metadata = json.load(f)
        sample_names_in_db = set(s['sample_name'] for s in db_metadata['callsets'])
        sample_names_requested = set([s['id'] for s in samples])
        sample_names_to_add = sample_names_requested - sample_names_in_db
        sample_names_to_remove = sample_names_in_db - sample_names_requested
        if sample_names_to_remove:
            # GenomicsDB doesn't support removing, so creating a new DB
            updating_existing_db = False
            sample_names_already_added = set()
            sample_names_to_add = {s['id'] for s in samples}
            sample_names_will_be_in_db = sample_names_to_add
            logger.info(
                f'GenomicDB {genomicsdb_gcs_path} exists, but '
                f'{len(sample_names_to_remove)} samples need '
                f'to be removed: {", ".join(sample_names_to_remove)}, so creating a new '
                f'DB with {len(sample_names_will_be_in_db)} samples: '
                f'{", ".join(sample_names_will_be_in_db)}'
            )
        else:
            updating_existing_db = True
            sample_names_will_be_in_db = sample_names_in_db | sample_names_to_add
            sample_names_already_added = sample_names_requested & sample_names_in_db
            if sample_names_already_added:
                logger.info(
                    f'{len(sample_names_already_added)} samples '
                    f'{", ".join(sample_names_already_added)} already exist in the DB '
                    f'{genomicsdb_gcs_path}, skipping adding them.'
                )
            if sample_names_to_remove:
                logger.info(
                    f'There are {len(sample_names_to_remove)} samples that need to be '
                    f'removed from the DB {genomicsdb_gcs_path}: '
                    f'{", ".join(sample_names_to_remove)}. Re-creating the DB '
                    f'using the updated set of samples'
                )
            elif sample_names_to_add:
                logger.info(
                    f'Will add {len(sample_names_to_add)} samples '
                    f'{", ".join(sample_names_to_add)} into the DB {genomicsdb_gcs_path}'
                )
            else:
                logger.warning(
                    f'Nothing will be added into the DB {genomicsdb_gcs_path}'
                )
    else:
        # Initiate new DB
        sample_names_already_added = set()
        sample_names_to_add = {s['id'] for s in samples}
        sample_names_will_be_in_db = sample_names_to_add
        updating_existing_db = False
        logger.info(
            f'GenomicDB {genomicsdb_gcs_path} doesn\'t exist, so creating a new one '
            f'with {len(sample_names_to_add)} samples: {", ".join(sample_names_to_add)}'
        )

    sample_map_bucket_path = join(tmp_bucket, 'work', 'sample_name.csv')
    local_sample_map_fpath = join(local_tmp_dir, 'sample_name.csv')
    with open(local_sample_map_fpath, 'w') as f:
        for sid in sample_names_to_add:
            f.write('\t'.join([sid, gvcf_by_sid[sid]]) + '\n')
    subprocess.run(
        f'gsutil cp {local_sample_map_fpath} {sample_map_bucket_path}',
        check=False,
        shell=True,
    )

    return (
        sample_names_to_add,
        sample_names_will_be_in_db,
        sample_names_already_added,
        updating_existing_db,
        sample_map_bucket_path,
    )


def _add_import_gvcfs_job(
    b: hb.Batch,
    genomicsdb_gcs_path: str,
    sample_names_to_add: Set[str],
    sample_names_to_skip: Set[str],
    sample_names_will_be_in_db: Set[str],
    updating_existing_db: bool,
    sample_map_bucket_path: str,
    interval: hb.ResourceFile,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[Optional[Job], Set[str]]:
    """
    Add GVCFs to a genomics database (or create a new instance if it doesn't exist)
    Returns a Job, or None if no new samples to add
    """
    if updating_existing_db:
        # Update existing DB
        genomicsdb_param = f'--genomicsdb-update-workspace-path {genomicsdb_gcs_path}'
        job_name = 'Adding to GenomicsDB'
    else:
        # Initiate new DB
        genomicsdb_param = f'--genomicsdb-workspace-path {genomicsdb_gcs_path}'
        job_name = 'Creating GenomicsDB'

    sample_map = b.read_input(sample_map_bucket_path)

    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(16)
    java_mem = 16
    j.memory('lowmem')  # ~ 1G/core ~ 14.4G
    if depends_on:
        j.depends_on(*depends_on)

    j.declare_resource_group(output={'tar': '{root}.tar'})
    j.command(
        f"""set -e

    # We've seen some GenomicsDB performance regressions related to intervals, 
    # so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg. There's no data in between since 
    # we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!

    (while true; do df -h; pwd; du -sh $(dirname {j.output['tar']}); free -m; sleep 300; done) &

    echo "Adding {len(sample_names_to_add)} samples: {', '.join(sample_names_to_add)}"
    {f'echo "Skipping adding {len(sample_names_to_skip)} samples that are already in the DB: '
     f'{", ".join(sample_names_to_skip)}"' if sample_names_to_skip else ''}

    gatk --java-options -Xms{java_mem}g \\
      GenomicsDBImport \\
      {genomicsdb_param} \\
      --batch-size 50 \\
      -L {interval} \\
      --sample-name-map {sample_map} \\
      --reader-threads {java_mem} \\
      --merge-input-intervals \\
      --consolidate

    df -h; pwd; du -sh $(dirname {j.output['tar']}); free -m
    """
    )
    return j, sample_names_will_be_in_db


def _add_genotype_gvcfs_job(
    b: hb.Batch,
    genomicsdb_path: str,
    reference: hb.ResourceGroup,
    dbsnp: str,
    overwrite: bool,
    number_of_samples: int,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    interval: Optional[hb.ResourceFile] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Run joint-calling on all samples in a genomics database
    """
    job_name = 'Joint genotyping: GenotypeGVCFs'
    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(2)
    j.memory('standard')  # ~ 4G/core ~ 8G
    # 4G (fasta+fai+dict) + 4G per sample divided by the number of intervals
    j.storage(f'{4 + number_of_samples * 4 // number_of_intervals}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""
set -o pipefail
set -ex

cd $(dirname {j.output_vcf})

export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS
gsutil -q cp -r {genomicsdb_path} .

(while true; do df -h; pwd; free -m; sleep 300; done) &

df -h; pwd; free -m

gatk --java-options -Xms8g \\
  GenotypeGVCFs \\
  -R {reference.base} \\
  -O {j.output_vcf['vcf.gz']} \\
  -D {dbsnp} \\
  --only-output-calls-starting-in-intervals \\
  -V gendb://{basename(genomicsdb_path)} \\
  {f'-L {interval} ' if interval else ''} \\
  --merge-input-intervals

df -h; pwd; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_gnarly_genotyper_job(
    b: hb.Batch,
    genomicsdb_path: str,
    reference: hb.ResourceGroup,
    dbsnp: str,
    overwrite: bool,
    number_of_samples: int,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Runs GATK GnarlyGenotyper on a combined_gvcf VCF bgzipped file.

    GnarlyGenotyper performs "quick and dirty" joint genotyping on large cohorts,
    pre-called with HaplotypeCaller, and post-processed with ReblockGVCF.

    HaplotypeCaller must be used with `-ERC GVCF` or `-ERC BP_RESOLUTION` to add
    genotype likelihoods.

    ReblockGVCF must be run to add all the annotations necessary for VQSR:
    QUALapprox, VarDP, RAW_MQandDP.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: GnarlyGenotyper'
    if interval_idx is not None:
        job_name += f' {interval_idx}/{number_of_intervals}'

    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(2)
    j.memory('standard')  # ~ 4G/core ~ 8G
    # 4G (fasta+fai+dict) + 4G per sample divided by the number of intervals
    j.storage(f'{4 + number_of_samples * 4 // number_of_intervals}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    j.command(
        f"""
set -o pipefail
set -ex

cd $(dirname {j.output_vcf})

export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS
gsutil -q cp -r {genomicsdb_path} .
    
(while true; do df -h; pwd; free -m; sleep 300; done) &

df -h; pwd; free -m

gatk --java-options -Xms8g \\
  GnarlyGenotyper \\
  -R {reference.base} \\
  -O {j.output_vcf['vcf.gz']} \\
  -D {dbsnp} \\
  --only-output-calls-starting-in-intervals \\
  --keep-all-sites \\
  -V gendb://{basename(genomicsdb_path)} \\
  {f'-L {interval} ' if interval else ''} \\
  --create-output-variant-index

df -h; pwd; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_exccess_het_filter(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    excess_het_threshold: float = 54.69,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Filter a large cohort callset on Excess Heterozygosity.

    The filter applies only to large callsets (`not is_small_callset`)

    Requires all samples to be unrelated.

    ExcessHet estimates the probability of the called samples exhibiting excess
    heterozygosity with respect to the null hypothesis that the samples are unrelated.
    The higher the score, the higher the chance that the variant is a technical artifact
    or that there is consanguinuity among the samples. In contrast to Inbreeding
    Coefficient, there is no minimal number of samples for this annotation.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: ExcessHet filter'
    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'32G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    # Captring stderr to avoid Batch pod from crashing with OOM from millions of
    # warning messages from VariantFiltration, e.g.:
    # > JexlEngine - ![0,9]: 'ExcessHet > 54.69;' undefined variable ExcessHet
    gatk --java-options -Xms3g \\
      VariantFiltration \\
      --filter-expression 'ExcessHet > {excess_het_threshold}' \\
      --filter-name ExcessHet \\
      {f'-L {interval} ' if interval else ''} \\
      -O {j.output_vcf['vcf.gz']} \\
      -V {input_vcf['vcf.gz']} \\
      2> {j.stderr}
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_make_sites_only_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the AS-VQSR modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: MakeSitesOnlyVcf'
    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'32G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      MakeSitesOnlyVcf \\
      -I {input_vcf['vcf.gz']} \\
      -O {j.output_vcf['vcf.gz']}
      """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_final_gather_vcf_job(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    overwrite: bool,
    output_vcf_path: str = None,
) -> Job:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    job_name = f'Gather {len(input_vcfs)} VCFs'
    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{1 + len(input_vcfs) * 1}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    (while true; do df -h; pwd free -m; sleep 300; done) &

    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms{java_mem}g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}
    
    df -h; pwd; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


if __name__ == '__main__':
    main()  # pylint: disable=E1120
