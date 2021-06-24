"""
Provides find_inputs function that finds input files (.g.vcf.gz, .bam or .cram) 
in provided buckets, along with corresponding indices (tbi, bai, crai).
Compares sample names to the provided PED, and returns DataFrame.
"""

import logging
import os
import re
import subprocess
from os.path import join, basename
from typing import Optional, List, Dict, Tuple
import pandas as pd
from utils import file_exists

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def find_inputs(
    gvcf_buckets: List[str],
    bam_buckets: List[str],
    bam_to_realign_buckets: List[str],
    local_tmp_dir: str,
    work_bucket: str,
    ped_fpath: Optional[str] = None,
) -> Tuple[pd.DataFrame, str]:
    """
    Find input files (.g.vcf.gz, .bam or .cram) in provided buckets,
    along with corresponding indices (tbi, bai, crai).
    Compares sample names to the provided PED, and returns DataFrame.
    :param gvcf_buckets: buckets to find GVCF files
    :param bam_buckets: buckets to find BAM files
        (will be passed to HaplotypeCaller to produce GVCFs)
    :param bam_to_realign_buckets: buckets to find BAM files
        (will be re-aligned with BWA before passing to HaplotypeCaller)
    :param work_bucket: bucket for temporary files
    :param ped_fpath: pedigree file. If not provided, a bare one will be generated
        with the sample names derived from input file names, and with missing/unknown
        pedigree and sex information
    :param local_tmp_dir: temporary local directory
    :return: a tuple of a DataFrame with the pedigree information and paths to
        input files (columns: s, file, index, type), and a path to a PED file
    """
    input_buckets_by_type = dict(
        gvcf=gvcf_buckets,
        bam=bam_buckets,
        bam_to_realign=bam_to_realign_buckets,
    )
    found_files_by_type = _find_files_by_type(input_buckets_by_type)
    found_indices_by_type = _find_file_indices(found_files_by_type)

    if ped_fpath:
        # PED file provided, so creating DataFrame based on sample names in the file,
        # and comparing to found input files
        df = _df_based_on_ped_file(
            ped_fpath,
            found_files_by_type,
            found_indices_by_type,
            local_tmp_dir,
        )

    else:
        # PED file not provided, so creating DataFrame purely from the found input files
        data: Dict[str, List] = dict(
            fam_id=[],
            s=[],
            father_id=[],
            mother_id=[],
            sex=[],
            phenotype=[],
            file=[],
            index=[],
            type=[],
        )
        for input_type in found_files_by_type:
            for fp, index in zip(
                found_files_by_type[input_type], found_indices_by_type[input_type]
            ):
                sn = _get_file_base_name(fp)
                data['s'].append(sn)
                data['file'].append(fp)
                data['index'].append(index)
                data['type'].append(input_type)
                data['fam_id'].append(sn)
                data['father_id'].append('0')
                data['mother_id'].append('0')
                data['sex'].append('0')
                data['phenotype'].append('0')

        df = pd.DataFrame(data=data).set_index('s', drop=False)
        ped_fpath = join(work_bucket, 'bare.ped')
        df[['fam_id', 's', 'father_id', 'mother_id', 'sex', 'phenotype']].write(
            ped_fpath
        )

    return df, ped_fpath


def _find_files_by_type(
    input_buckets_by_type: Dict[str, List[str]]
) -> Dict[str, List[str]]:
    """
    Find input files like .g.vcf.gz, .bam, .cram
    """
    found_files_by_type: Dict[str, List[str]] = {t: [] for t in input_buckets_by_type}
    for input_type, buckets in input_buckets_by_type.items():
        patterns = ['*.g.vcf.gz'] if input_type == 'gvcf' else ['*.bam', '*.cram']
        for ib in buckets:
            for pattern in patterns:
                cmd = f"gsutil ls '{join(ib, pattern)}'"
                try:
                    found_files_by_type[input_type].extend(
                        line.strip()
                        for line in subprocess.check_output(cmd, shell=True)
                        .decode()
                        .split()
                    )
                except subprocess.CalledProcessError:
                    pass
    return found_files_by_type


def _find_file_indices(
    found_files_by_type: Dict[str, List[str]]
) -> Dict[str, List[str]]:
    """
    Find corresponding index files. Will look for:
    '<sample>.g.vcf.gz.tbi' for '<sample>.g.vcf.gz',
    '<sample>.bam.bai' or '<sample>.bai' for '<sample>.bam',
    '<sample>.cram.crai' or '<sample>.crai' for '<sample>.cram',
    """
    found_indices_by_type: Dict[str, List[str]] = {t: [] for t in found_files_by_type}
    for input_type, fpaths in found_files_by_type.items():
        for fp in fpaths:
            index = fp + '.tbi'
            if fp.endswith('.g.vcf.gz'):
                if not file_exists(index):
                    logger.critical(f'Not found TBI index for file {fp}')
            elif fp.endswith('.bam'):
                index = fp + '.bai'
                if not file_exists(index):
                    index = re.sub('.bam$', '.bai', fp)
                    if not file_exists(index):
                        logger.critical(f'Not found BAI index for file {fp}')
            elif fp.endswith('.cram'):
                index = fp + '.crai'
                if not file_exists(index):
                    index = re.sub('.cram$', '.crai', fp)
                    if not file_exists(index):
                        logger.critical(f'Not found CRAI index for file {fp}')
            else:
                logger.critical(f'Unrecognised input file extention {fp}')
            found_indices_by_type[input_type].append(index)
    return found_indices_by_type


def _df_based_on_ped_file(
    ped_fpath: str,
    found_files_by_type: Dict[str, List[str]],
    found_indices_by_type: Dict[str, List[str]],
    local_tmp_dir: str,
) -> pd.DataFrame:
    """
    Compares found input files to the samples provided in the PED file.
    Creates an input DataFrame using the PED file as a base.
    """
    local_sample_list_fpath = join(local_tmp_dir, basename(ped_fpath))
    subprocess.run(
        f'gsutil cat {ped_fpath} | grep -v ^Family.ID > {local_sample_list_fpath}',
        check=False,
        shell=True,
    )
    df = pd.read_csv(
        local_sample_list_fpath,
        delimiter='\t',
        names=['fam_id', 's', 'father_id', 'mother_id', 'sex', 'phenotype'],
    ).set_index('s', drop=False)
    ped_snames = list(df['s'])

    # Checking that file base names have a 1-to-1 match with the samples in PED
    # First, checking the match of PED sample names to input files
    all_input_snames = []
    for fpaths in found_files_by_type.values():
        all_input_snames.extend([_get_file_base_name(fp) for fp in fpaths])

    mismatches_with_ped_found = False

    for ped_sname in ped_snames:
        matching_files = [
            input_sname for input_sname in all_input_snames if ped_sname == input_sname
        ]
        if len(matching_files) > 1:
            logging.warning(
                f'Multiple input files found for the sample {ped_sname}:'
                f'{matching_files}'
            )
            mismatches_with_ped_found = True
        elif len(matching_files) == 0:
            logging.warning(f'No files found for the sample {ped_sname}')
            mismatches_with_ped_found = True

    # Second, checking the match of input files to PED sample names, and filling a dict
    for input_type in found_files_by_type:
        for fp, index in zip(
            found_files_by_type[input_type], found_indices_by_type[input_type]
        ):
            input_sname = _get_file_base_name(fp)
            matching_sn = [sn for sn in ped_snames if sn == input_sname]
            if len(matching_sn) > 1:
                logging.warning(
                    f'Multiple PED records found for the input {input_sname}:'
                    f'{matching_sn}'
                )
                mismatches_with_ped_found = True
            elif len(matching_sn) == 0:
                logging.warning(f'No PED records found for the input {input_sname}')
                mismatches_with_ped_found = True
            else:
                df.loc[matching_sn[0], 'type'] = input_type
                df.loc[matching_sn[0], 'file'] = fp
                df.loc[matching_sn[0], 'index'] = index

    if mismatches_with_ped_found:
        logging.critical(
            'Mismatches found between provided input bucket contents, '
            'and the provided PED file. Make sure all files in the buckets are '
            'named after sample names in the PED file as '
            '{sample}.cram, {sample}.g.vcf.gz, etc'
        )
    df = df[df.file.notnull()]
    return df


def _get_file_base_name(file_path):
    """
    Strips directory path and extention from a file name. Effectively extractss
    assumed sample name (i.e. assumes file named as gs://path/to/{sample}.g.vcf.gz)
    """
    return re.sub('(.bam|.cram|.g.vcf.gz)$', '', os.path.basename(file_path))
