"""
Create jobs to create and apply a VQSR model
"""

import os
from typing import List, Optional, Dict
import logging
import hailtop.batch as hb
from hailtop.batch.job import Job
from seqr_loader import utils
from seqr_loader.utils import BROAD_REF_BUCKET

logger = logging.getLogger('VQSR')
logger.setLevel('INFO')


SNP_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.8,
    99.6,
    99.5,
    99.4,
    99.3,
    99.0,
    98.0,
    97.0,
    90.0,
]
SNP_HARD_FILTER_LEVEL = 99.7
SNP_RECALIBRATION_ANNOTATION_VALUES = [
    'QD',
    'MQRankSum',
    'ReadPosRankSum',
    'FS',
    'MQ',
    'SOR',
    'DP',
]
SNP_RECALIBRATION_ANNOTATION_VALUES_AS = [
    'AS_QD',
    'AS_MQRankSum',
    'AS_ReadPosRankSum',
    'AS_FS',
    'AS_SOR',
    'AS_MQ',
]
INDEL_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.5,
    99.0,
    97.0,
    96.0,
    95.0,
    94.0,
    93.5,
    93.0,
    92.0,
    91.0,
    90.0,
]
INDEL_HARD_FILTER_LEVEL = 99.7
INDEL_RECALIBRATION_ANNOTATION_VALUES = [
    'FS',
    'ReadPosRankSum',
    'MQRankSum',
    'QD',
    'SOR',
    'DP',
]
INDEL_RECALIBRATION_ANNOTATION_VALUES_AS = [
    'AS_FS',
    'AS_SOR',
    'AS_ReadPosRankSum',
    'AS_MQRankSum',
    'AS_QD',
]


def make_vqsr_jobs(
    b: hb.Batch,
    input_vcf_gathered: str,
    is_small_callset: bool,
    is_huge_callset: bool,
    work_bucket: str,
    web_bucket: str,
    depends_on: Optional[List[Job]],
    intervals: Dict,
    scatter_count: int,
    output_vcf_path: str,
    use_as_annotations: bool,
    overwrite: bool,
) -> Job:
    """
    Add jobs that perform the allele-specific VQSR variant QC

    :param b: Batch object to add jobs to
    :param input_vcf_gathered: path to an input gathered VCF
    :param is_small_callset: for small callsets, we gather the VCF shards and collect
        QC metrics directly. For anything larger, we need to keep the VCF sharded and
        gather metrics collected from them
    :param is_huge_callset: For huge callsets, we allocate more memory for the SNPs
        Create Model step
    :param work_bucket: bucket for intermediate files
    :param web_bucket: bucket for plots and evaluation results (exposed via http)
    :param depends_on: job that the created jobs should only run after
    :param intervals: ResourceGroup object with intervals to scatter
    :param scatter_count: number of interavals
    :param output_vcf_path: path to write final recalibrated VCF to
    :param use_as_annotations: use allele-specific annotation for VQSR
    :param overwrite: whether to not reuse intermediate files
    :return: a final Job, and a path to the VCF with VQSR annotations
    """

    # Reference files. All options have defaults.
    dbsnp_vcf = os.path.join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf')
    dbsnp_vcf_index = os.path.join(
        BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf.idx'
    )
    hapmap_resource_vcf = os.path.join(BROAD_REF_BUCKET, 'hapmap_3.3.hg38.vcf.gz')
    hapmap_resource_vcf_index = os.path.join(
        BROAD_REF_BUCKET, 'hapmap_3.3.hg38.vcf.gz.tbi'
    )
    omni_resource_vcf = os.path.join(BROAD_REF_BUCKET, '1000G_omni2.5.hg38.vcf.gz')
    omni_resource_vcf_index = os.path.join(
        BROAD_REF_BUCKET, '1000G_omni2.5.hg38.vcf.gz.tbi'
    )
    one_thousand_genomes_resource_vcf = os.path.join(
        BROAD_REF_BUCKET, '1000G_phase1.snps.high_confidence.hg38.vcf.gz'
    )
    one_thousand_genomes_resource_vcf_index = os.path.join(
        BROAD_REF_BUCKET, '1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi'
    )
    mills_resource_vcf = os.path.join(
        BROAD_REF_BUCKET, 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    )
    mills_resource_vcf_index = os.path.join(
        BROAD_REF_BUCKET, 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
    )
    axiom_poly_resource_vcf = os.path.join(
        BROAD_REF_BUCKET, 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz'
    )
    axiom_poly_resource_vcf_index = os.path.join(
        BROAD_REF_BUCKET,
        'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi',
    )
    dbsnp_vcf = b.read_input_group(base=dbsnp_vcf, index=dbsnp_vcf_index)
    hapmap_resource_vcf = b.read_input_group(
        base=hapmap_resource_vcf, index=hapmap_resource_vcf_index
    )
    omni_resource_vcf = b.read_input_group(
        base=omni_resource_vcf, index=omni_resource_vcf_index
    )
    one_thousand_genomes_resource_vcf = b.read_input_group(
        base=one_thousand_genomes_resource_vcf,
        index=one_thousand_genomes_resource_vcf_index,
    )
    mills_resource_vcf = b.read_input_group(
        base=mills_resource_vcf, index=mills_resource_vcf_index
    )
    axiom_poly_resource_vcf = b.read_input_group(
        base=axiom_poly_resource_vcf, index=axiom_poly_resource_vcf_index
    )
    dbsnp_resource_vcf = dbsnp_vcf

    # To fit only a site-only VCF
    if is_small_callset:
        small_disk = 50
    elif not is_huge_callset:
        small_disk = 100
    else:
        small_disk = 200

    # To fit a joint-called VCF
    if is_small_callset:
        medium_disk = 100
    elif not is_huge_callset:
        medium_disk = 200
    else:
        medium_disk = 500

    if is_small_callset:
        huge_disk = 200
    elif not is_huge_callset:
        huge_disk = 500
    else:
        huge_disk = 2000

    site_only_j = _add_make_sites_only_job(
        b=b,
        input_vcf=b.read_input_group(
            **{
                'vcf.gz': input_vcf_gathered,
                'vcf.gz.tbi': input_vcf_gathered + '.tbi',
            }
        ),
        overwrite=overwrite,
        disk=medium_disk,
    )
    gathered_vcf = site_only_j.output_vcf
    first_job = site_only_j
    if depends_on:
        first_job.depends_on(*depends_on)

    indels_variant_recalibrator_job = add_indels_variant_recalibrator_job(
        b,
        sites_only_variant_filtered_vcf=gathered_vcf,
        mills_resource_vcf=mills_resource_vcf,
        axiom_poly_resource_vcf=axiom_poly_resource_vcf,
        dbsnp_resource_vcf=dbsnp_resource_vcf,
        disk_size=small_disk,
        use_as_annotations=use_as_annotations,
        work_bucket=web_bucket,
    )
    indels_recalibration = indels_variant_recalibrator_job.recalibration
    indels_tranches = indels_variant_recalibrator_job.tranches

    snp_max_gaussians = 6
    if is_small_callset:
        snp_max_gaussians = 4
    elif is_huge_callset:
        snp_max_gaussians = 8

    if is_huge_callset:
        # Run SNP recalibrator in a scattered mode
        model_file = add_snps_variant_recalibrator_create_model_step(
            b,
            sites_only_variant_filtered_vcf=gathered_vcf,
            hapmap_resource_vcf=hapmap_resource_vcf,
            omni_resource_vcf=omni_resource_vcf,
            one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
            dbsnp_resource_vcf=dbsnp_resource_vcf,
            disk_size=small_disk,
            use_as_annotations=use_as_annotations,
            web_bucket=web_bucket,
            work_bucket=work_bucket,
            is_small_callset=is_small_callset,
            is_huge_callset=is_huge_callset,
            max_gaussians=snp_max_gaussians,
        ).model_file

        snps_recalibrator_jobs = [
            add_snps_variant_recalibrator_scattered_step(
                b,
                sites_only_vcf=gathered_vcf,
                interval=intervals[f'interval_{idx}'],
                model_file=model_file,
                hapmap_resource_vcf=hapmap_resource_vcf,
                omni_resource_vcf=omni_resource_vcf,
                one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
                dbsnp_resource_vcf=dbsnp_resource_vcf,
                disk_size=small_disk,
                use_as_annotations=use_as_annotations,
                max_gaussians=snp_max_gaussians,
            )
            for idx in range(scatter_count)
        ]
        snps_recalibrations = [j.recalibration for j in snps_recalibrator_jobs]
        snps_tranches = [j.tranches for j in snps_recalibrator_jobs]
        snps_gathered_tranches = add_snps_gather_tranches_step(
            b,
            tranches=snps_tranches,
            disk_size=small_disk,
        ).out_tranches

        scattered_vcfs = [
            add_apply_recalibration_step(
                b,
                input_vcf=gathered_vcf,
                interval=intervals[f'interval_{idx}'],
                indels_recalibration=indels_recalibration,
                indels_tranches=indels_tranches,
                snps_recalibration=snps_recalibrations[idx],
                snps_tranches=snps_gathered_tranches,
                disk_size=huge_disk,
                use_as_annotations=use_as_annotations,
                snp_filter_level=SNP_HARD_FILTER_LEVEL,
                indel_filter_level=INDEL_HARD_FILTER_LEVEL,
            ).recalibrated_vcf
            for idx in range(scatter_count)
        ]
        recalibrated_gathered_vcf_job = _add_final_gather_vcf_step(
            b,
            input_vcfs=scattered_vcfs,
            disk_size=huge_disk,
            output_vcf_path=output_vcf_path,
        )

    else:
        snps_recalibrator_job = add_snps_variant_recalibrator_step(
            b,
            sites_only_variant_filtered_vcf=gathered_vcf,
            hapmap_resource_vcf=hapmap_resource_vcf,
            omni_resource_vcf=omni_resource_vcf,
            one_thousand_genomes_resource_vcf=one_thousand_genomes_resource_vcf,
            dbsnp_resource_vcf=dbsnp_resource_vcf,
            disk_size=small_disk,
            use_as_annotations=use_as_annotations,
            max_gaussians=snp_max_gaussians,
            web_bucket=web_bucket,
            work_bucket=work_bucket,
        )
        snps_recalibration = snps_recalibrator_job.recalibration
        snps_tranches = snps_recalibrator_job.tranches

        recalibrated_gathered_vcf_job = add_apply_recalibration_step(
            b,
            input_vcf=gathered_vcf,
            indels_recalibration=indels_recalibration,
            indels_tranches=indels_tranches,
            snps_recalibration=snps_recalibration,
            snps_tranches=snps_tranches,
            disk_size=huge_disk,
            use_as_annotations=use_as_annotations,
            indel_filter_level=SNP_HARD_FILTER_LEVEL,
            snp_filter_level=INDEL_HARD_FILTER_LEVEL,
            output_vcf_path=output_vcf_path,
        )

    return recalibrated_gathered_vcf_job


def _add_make_sites_only_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    disk: int,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the AS-VQSR modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    job_name = 'VQSR: MakeSitesOnlyVcf'
    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'{disk}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      MakeSitesOnlyVcf \\
      -I {input_vcf['vcf.gz']} \\
      -O {j.output_vcf['vcf.gz']} \\
      --CREATE_INDEX
      """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def add_tabix_step(
    b: hb.Batch,
    vcf_path: str,
    disk_size: int,
) -> Job:
    """
    Regzip and tabix the combined VCF (for some reason the one output with mt2vcf
    is not block-gzipped)
    """
    j = b.new_job('VQSR: Tabix')
    j.image(utils.BCFTOOLS_IMAGE)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        combined_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    vcf_inp = b.read_input(vcf_path)
    j.command(
        f"""set -e
        gunzip {vcf_inp} -c | bgzip -c > {j.combined_vcf['vcf.gz']}
        tabix -p vcf {j.combined_vcf['vcf.gz']}
        """
    )
    return j


def add_split_intervals_step(
    b: hb.Batch,
    interval_list: hb.ResourceFile,
    scatter_count: int,
    ref_fasta: hb.ResourceGroup,
    disk_size: int,
) -> Job:
    """
    Split genome into intervals to parallelise GnarlyGenotyper.

    Returns: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job('VQSR: SplitIntervals')
    j.image(utils.GATK_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )

    j.command(
        f"""set -e

    # Modes other than INTERVAL_SUBDIVISION will produce an unpredictable number
    # of intervals. But we have to produce exactly {scatter_count} number of
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{mem_gb - 1}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta.base} \\
      -mode INTERVAL_SUBDIVISION
      """
    )
    return j


def add_sites_only_gather_vcf_step(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceFile],
    disk_size: int,
) -> Job:
    """
    Gathers VCF files from scattered operations into a single VCF file

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    j = b.new_job('VQSR: SitesOnlyGatherVcf')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'{disk_size}G')

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in
    # our invocation. This argument disables expensive checks that the file headers
    # contain the same set of genotyped samples and that files are in order by position
    # of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}"""
    )
    return j


def add_indels_variant_recalibrator_job(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    mills_resource_vcf: hb.ResourceGroup,
    axiom_poly_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    work_bucket: str = None,
    max_gaussians: int = 4,
) -> Job:
    """
    Run VariantRecalibrator to calculate VQSLOD tranches for indels

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value. 4 is a
    reasonable default for indels, as their number is smaller than SNPs.

    Returns: a Job object with 3 outputs: j.recalibration (ResourceGroup), j.tranches,
    and j.indel_rscript_file. The latter is usedful to produce the optional tranche plot.
    """
    j = b.new_job('VQSR: IndelsVariantRecalibrator')
    j.image(utils.GATK_IMAGE)
    j.memory('highmem')
    ncpu = 4  # ~ 8G/core ~ 32G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 4
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join(
        [f'-tranche {v}' for v in INDEL_RECALIBRATION_TRANCHE_VALUES]
    )
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
                INDEL_RECALIBRATION_ANNOTATION_VALUES_AS
                if use_as_annotations
                else INDEL_RECALIBRATION_ANNOTATION_VALUES
            )
        ]
    )
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{java_mem}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode INDEL \\
      {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
      --max-gaussians {max_gaussians} \\
      -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiom_poly_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.indel_rscript_file}
      
      ls $(dirname {j.indel_rscript_file})
      """
    )
    if work_bucket:
        b.write_output(
            j.indel_rscript_file,
            os.path.join(work_bucket, 'recalibration-indels-features.Rscript'),
        )
    return j


def add_snps_variant_recalibrator_create_model_step(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    web_bucket: str = None,
    work_bucket: str = None,
    is_small_callset: bool = False,
    is_huge_callset: bool = False,
    max_gaussians: int = 4,
) -> Job:
    """
    First step of VQSR for SNPs: run VariantRecalibrator to subsample variants
    and produce a file of the VQSR model.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    Returns: a Job object with 2 outputs: j.model and j.snp_rscript_file.
    The latter is useful to produce the optional tranche plot.
    """
    j = b.new_job('VQSR: SNPsVariantRecalibratorCreateModel')
    j.image(utils.GATK_IMAGE)
    j.memory('highmem')
    if is_small_callset:
        ncpu = 8  # ~ 8G/core ~ 64G
    else:
        ncpu = 16  # ~ 8G/core ~ 128G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 10
    j.storage(f'{disk_size}G')

    downsample_factor = 75 if is_huge_callset else 10

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
                SNP_RECALIBRATION_ANNOTATION_VALUES_AS
                if use_as_annotations
                else SNP_RECALIBRATION_ANNOTATION_VALUES
            )
        ]
    )
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{java_mem}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
      --sample-every-Nth-variant {downsample_factor} \\
      --output-model {j.model_file} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.snp_rscript}

      ls $(dirname {j.snp_rscript})

      ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}
      ln {j.tranches}.pdf {j.tranches_pdf}
      """
    )
    if work_bucket:
        b.write_output(
            j.snp_rscript,
            os.path.join(work_bucket, 'recalibration-snps-features.RScript'),
        )
    if web_bucket:
        b.write_output(
            j.snp_rscript_pdf,
            os.path.join(web_bucket, 'recalibration-snps-features.pdf'),
        )
        b.write_output(
            j.tranches_pdf,
            os.path.join(web_bucket, 'recalibration-snps-tranches.pdf'),
        )
    return j


def add_snps_variant_recalibrator_scattered_step(
    b: hb.Batch,
    sites_only_vcf: hb.ResourceGroup,
    model_file: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    interval: Optional[hb.ResourceGroup] = None,
    max_gaussians: int = 4,
) -> Job:
    """
    Second step of VQSR for SNPs: run VariantRecalibrator scattered to apply
    the VQSR model file to each genomic interval.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    Returns: a Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches
    """
    j = b.new_job('VQSR: SNPsVariantRecalibratorScattered')

    j.image(utils.GATK_IMAGE)
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f'{mem_gb}G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
                SNP_RECALIBRATION_ANNOTATION_VALUES_AS
                if use_as_annotations
                else SNP_RECALIBRATION_ANNOTATION_VALUES
            )
        ]
    )
    j.command(
        f"""set -euo pipefail

    MODEL_REPORT={model_file}

    gatk --java-options -Xms{mem_gb - 1}g \\
      VariantRecalibrator \\
      -V {sites_only_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {f'-L {interval} ' if interval else ''} \\
      {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
      --input-model {model_file} --output-tranches-for-scatter \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base}"""
    )
    return j


def add_snps_variant_recalibrator_step(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    web_bucket: str,
    work_bucket: str,
    disk_size: int,
    use_as_annotations: bool,
    max_gaussians: int = 4,
) -> Job:
    """
    Recalibrate SNPs in one run (alternative to scatter-gather approach)
    """
    j = b.new_job('VQSR: SNPsVariantRecalibrator')

    j.image(utils.GATK_IMAGE)
    j.memory('highmem')
    ncpu = 8  # ~ 8G/core ~ 64G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 8
    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
                SNP_RECALIBRATION_ANNOTATION_VALUES_AS
                if use_as_annotations
                else SNP_RECALIBRATION_ANNOTATION_VALUES
            )
        ]
    )
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{java_mem}g \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} \\
      --rscript-file {j.snp_rscript}

      ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}
      ln {j.tranches}.pdf {j.tranches_pdf}
      """
    )

    if work_bucket:
        b.write_output(
            j.snp_rscript,
            os.path.join(work_bucket, 'recalibration-snps-features.RScript'),
        )
    if web_bucket:
        b.write_output(
            j.snp_rscript_pdf,
            os.path.join(web_bucket, 'recalibration-snps-features.pdf'),
        )
        b.write_output(
            j.tranches_pdf,
            os.path.join(web_bucket, 'recalibration-snps-tranches.pdf'),
        )
    return j


def add_snps_gather_tranches_step(
    b: hb.Batch,
    tranches: List[hb.ResourceFile],
    disk_size: int,
) -> Job:
    """
    Third step of VQSR for SNPs: run GatherTranches to gather scattered per-interval
    tranches outputs.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    Returns: a Job object with one output j.out_tranches
    """
    j = b.new_job('VQSR: SNPGatherTranches')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.cpu(2)
    j.storage(f'{disk_size}G')

    inputs_cmdl = ' '.join([f'--input {t}' for t in tranches])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      GatherTranches \\
      --mode SNP \\
      {inputs_cmdl} \\
      --output {j.out_tranches}"""
    )
    return j


def add_apply_recalibration_step(
    b: hb.Batch,
    input_vcf: hb.ResourceFile,
    indels_recalibration: hb.ResourceGroup,
    indels_tranches: hb.ResourceFile,
    snps_recalibration: hb.ResourceGroup,
    snps_tranches: hb.ResourceFile,
    disk_size: int,
    use_as_annotations: bool,
    indel_filter_level: float,
    snp_filter_level: float,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Apply a score cutoff to filter variants based on a recalibration table.
    Runs ApplyVQSR twice to apply first indel, then SNP recalibrations.

    Targets indel_filter_level and snp_filter_level sensitivities. The tool matches
    them internally to a VQSLOD score cutoff based on the model's estimated sensitivity
    to a set of true variants.

    The filter determination is not just a pass/fail process. The tool evaluates for
    each variant which "tranche", or slice of the dataset, it falls into in terms of
    sensitivity to the truthset. Variants in tranches that fall below the specified
    truth sensitivity filter level have their FILTER field annotated with the
    corresponding tranche level. This results in a callset that is filtered to the
    desired level but retains the information necessary to increase sensitivity
    if needed.

    Returns: a Job object with one ResourceGroup output j.output_vcf, correponding
    to a VCF with tranche annotated in the FILTER field
    """
    j = b.new_job('VQSR: ApplyRecalibration')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O tmp.indel.recalibrated.vcf \\
      -V {input_vcf['vcf.gz']} \\
      --recal-file {indels_recalibration} \\
      --tranches-file {indels_tranches} \\
      --truth-sensitivity-filter-level {indel_filter_level} \\
      --create-output-variant-index true \\
      {f'-L {interval} ' if interval else ''} \\
      {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
      -mode INDEL
      
    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})

    rm {input_vcf['vcf.gz']} {indels_recalibration} {indels_tranches}

    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})

    gatk --java-options -Xms5g \\
      ApplyVQSR \\
      -O {j.output_vcf['vcf.gz']} \\
      -V tmp.indel.recalibrated.vcf \\
      --recal-file {snps_recalibration} \\
      --tranches-file {snps_tranches} \\
      --truth-sensitivity-filter-level {snp_filter_level} \\
      --create-output-variant-index true \\
      {f'-L {interval} ' if interval else ''} \\
      {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
      -mode SNP

    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})
      """
    )

    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def add_collect_metrics_sharded_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    dbsnp_vcf: hb.ResourceGroup,
    interval_list: hb.ResourceFile,
    ref_dict: hb.ResourceFile,
    disk_size: int,
):
    """
    Run CollectVariantCallingMetrics for site-level evaluation.

    This produces detailed and summary metrics report files. The summary metrics
    provide cohort-level variant metrics and the detailed metrics segment variant
    metrics for each sample in the callset. The detail metrics give the same metrics
    as the summary metrics for the samples plus several additional metrics.

    These are explained in detail at
    https://broadinstitute.github.io/picard/picard-metric-definitions.html.

    Returns: Job object with a single ResourceGroup output j.metrics, with
    j.metrics.detail_metrics and j.metrics.summary_metrics ResourceFiles
    """
    j = b.new_job('VQSR: CollectMetricsSharded')
    j.image(utils.GATK_IMAGE)
    j.memory('8G')
    j.cpu(2)
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        metrics={
            'detail_metrics': '{root}.variant_calling_detail_metrics',
            'summary_metrics': '{root}.variant_calling_summary_metrics',
        }
    )

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      CollectVariantCallingMetrics \\
      --INPUT {input_vcf['vcf.gz']} \\
      --DBSNP {dbsnp_vcf.base} \\
      --SEQUENCE_DICTIONARY {ref_dict} \\
      --OUTPUT {j.metrics} \\
      --THREAD_COUNT 8 \\
      --TARGET_INTERVALS {interval_list}"""
    )
    return j


def _add_final_gather_vcf_step(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    disk_size: int,
    output_vcf_path: str = None,
) -> Job:
    """
    Combines recalibrated VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    j = b.new_job('VQSR: FinalGatherVcf')
    j.image(utils.GATK_IMAGE)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        f"""set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms6g \\
      GatherVcfsCloud \\
      --ignore-safety-checks \\
      --gather-type BLOCK \\
      {input_cmdl} \\
      --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}"""
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def _add_final_filter_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    output_vcf_path: str = None,
) -> Job:
    """
    Hard-filters the VQSR'ed VCF
    """
    j = b.new_job('VQSR: final filter')
    j.image(utils.BCFTOOLS_IMAGE)
    j.memory(f'8G')
    j.storage(f'100G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""
    df -h; pwd; free -m

    bcftools view -f.,PASS {input_vcf['vcf.gz']} -Oz -o {j.output_vcf['vcf.gz']} && \\
    tabix {j.output_vcf['vcf.gz']}

    df -h; pwd; free -m
    """
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def _add_variant_eval_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    ref_fasta: hb.ResourceGroup,
    dbsnp_vcf: hb.ResourceGroup,
    disk_size: int,
    output_path: str = None,
) -> Job:
    """
    Run VariantEval for site-level evaluation.
    Saves the QC to `output_path` bucket
    """
    j = b.new_job('VQSR: VariantEval')
    j.image(utils.GATK_IMAGE)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms6g \\
      VariantEval \\
      --eval {input_vcf['vcf.gz']} \\
      -R {ref_fasta.base} \\
      -D {dbsnp_vcf.base} \\
      --output {j.output}"""
    )
    if output_path:
        b.write_output(j.output, output_path)
    return j


def add_gather_variant_calling_metrics_step(
    b: hb.Batch,
    input_details: List[hb.ResourceGroup],
    input_summaries: List[hb.ResourceGroup],
    disk_size: int,
    output_path_prefix: str = None,
) -> Job:
    """
    Combines metrics from multiple CollectVariantCallingMetrics runs.

    Returns: Job object with a single ResourceGroup output j.metrics, with
    j.metrics.detail_metrics and j.metrics.summary_metrics ResourceFiles

    Saves the QC results to a bucket with the `output_path_prefix` prefix
    """
    j = b.new_job('VQSR: GatherVariantCallingMetrics')
    j.image(utils.GATK_IMAGE)
    j.memory(f'8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        metrics={
            'detail_metrics': '{root}.variant_calling_detail_metrics',
            'summary_metrics': '{root}.variant_calling_summary_metrics',
        }
    )

    input_cmdl = ' '.join('--INPUT {f} ' for f in input_details + input_summaries)
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms2g \\
      AccumulateVariantCallingMetrics \\
      {input_cmdl} \\
      --OUTPUT {j.metrics}"""
    )
    if output_path_prefix:
        b.write_output(j.metrics, output_path_prefix)
    return j
