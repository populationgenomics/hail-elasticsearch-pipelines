"""
Provides find_inputs function that finds input files (.g.vcf.gz, .bam or .cram) 
in provided buckets, along with corresponding indices (tbi, bai, crai).
Compares sample names to the provided PED, and returns DataFrame.
"""

import logging
import os
import re
import subprocess
from os.path import join, basename, splitext
from typing import Optional, List, Dict, Tuple
import pandas as pd
from utils import file_exists

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def find_inputs(
    gvcfs: List[str],
    crams: List[str],
    data_to_realign: List[str],
    local_tmp_dir: str,
    work_bucket: str,
    ped_fpath: Optional[str] = None,
) -> Tuple[pd.DataFrame, str]:
    """
    Find input files (.g.vcf.gz, .bam or .cram) in provided buckets,
    along with corresponding indices (tbi, bai, crai).
    Compares sample names to the provided PED, and returns DataFrame.
    :param gvcfs: glob paths to find GVCF files
    :param crams: glob paths to find CRAM/BAM files
        (will be passed to HaplotypeCaller to produce GVCFs)
    :param data_to_realign: buckets to find CRAM/BAM files
        (will be re-aligned with BWA before passing to HaplotypeCaller)
    :param work_bucket: bucket for temporary files
    :param ped_fpath: pedigree file. If not provided, a bare one will be generated
        with the sample names derived from input file names, and with missing/unknown
        pedigree and sex information
    :param local_tmp_dir: temporary local directory
    :return: a tuple of a DataFrame with the pedigree information and paths to
        input files (columns: s, file, index, type), and a path to a PED file
    """
    gvcfs, crams, fastq_to_align, cram_to_realign = _find_files_by_type(
        gvcfs, crams, data_to_realign
    )
    gvcf_with_index_by_by_sn = _find_file_indices(gvcfs)
    cram_with_index_by_by_sn = _find_file_indices(crams)
    fastq_pairs_by_by_sn = _find_fastq_pairs(fastq_to_align)
    cram_with_index_to_realign_by_sn = _find_file_indices(cram_to_realign)

    data: Dict[str, List] = {
        'Family.ID': [],
        's': [],
        'Paternal.ID': [],
        'Maternal.ID': [],
        'Sex': [],
        'Phenotype': [],
        'file': [],
        'file2': [],
        'index': [],
        'type': [],
    }
    all_sns = (
        set(gvcf_with_index_by_by_sn.keys())
        | set(cram_with_index_by_by_sn.keys())
        | set(fastq_pairs_by_by_sn.keys())
        | set(cram_with_index_to_realign_by_sn.keys())
    )
    for sn in all_sns:
        data['s'].append(sn)
        data['Family.ID'].append(sn)
        data['Paternal.ID'].append('0')
        data['Maternal.ID'].append('0')
        data['Sex'].append('0')
        data['Phenotype'].append('0')
        if sn in gvcf_with_index_by_by_sn:
            file, index = gvcf_with_index_by_by_sn[sn]
            data['file'].append(file)
            data['file2'].append(None)
            data['index'].append(index)
            data['type'].append('gvcf')
        if sn in cram_with_index_by_by_sn:
            file, index = cram_with_index_by_by_sn[sn]
            data['file'].append(file)
            data['file2'].append(None)
            data['index'].append(index)
            data['type'].append('cram')
        if sn in cram_with_index_to_realign_by_sn:
            file, index = cram_with_index_to_realign_by_sn[sn]
            data['file'].append(file)
            data['file2'].append(None)
            data['index'].append(index)
            data['type'].append('cram_to_realign')
        if sn in cram_with_index_to_realign_by_sn:
            file1, file2 = fastq_pairs_by_by_sn[sn]
            data['file'].append(file1)
            data['file2'].append(file2)
            data['index'].append(None)
            data['type'].append('fastq_to_realign')

    df = pd.DataFrame(data=data).set_index('s', drop=False)

    # If PED file is provided, adding it and comparing to the found input files
    if ped_fpath:
        df = _add_ped_info(df, ped_fpath, local_tmp_dir)

    new_ped_fpath = join(work_bucket, 'data.ped')
    df.to_csv(new_ped_fpath, sep='\t', header=False, index=False)
    return df, new_ped_fpath


def _find_files_by_type(
    gvcf_globs: List[str],
    crams_globs: List[str],
    data_to_realign_globs: List[str],
) -> Tuple[List[str], List[str], List[str], List[str]]:
    """
    Find input files like .g.vcf.gz, .bam, .cram, .fastq
    Return 4 lists of paths: gvcfs, crams, fastq_to_align, cram_to_realign
    """
    input_globs_by_type = dict(
        gvcf=gvcf_globs,
        cram=crams_globs,
        data_to_realign=data_to_realign_globs,
    )

    expected_ext_by_type = dict(
        gvcf='.g.vcf.gz',
        cram=['.cram', '.bam'],
        data_to_realign=['.cram', '.bam', '.fq', '.fastq', '.fq.gz', '.fastq.gz'],
    )

    gvcfs = []
    crams = []
    fastq_to_align = []
    cram_to_realign = []

    for input_type, glob_paths in input_globs_by_type.items():
        for glob_path in glob_paths:
            assert any(
                glob_path.endswith(ext) for ext in expected_ext_by_type[input_type]
            ), (input_type, glob_path)
            glob_path = glob_path.lstrip("'").rstrip("'")
            cmd = f"gsutil ls '{glob_path}'"
            found_files = [
                line.strip()
                for line in subprocess.check_output(cmd, shell=True).decode().split()
            ]
            if input_type == 'data_to_realign':
                if glob_path.endswith('.cram') or glob_path.endswith('.bam'):
                    cram_to_realign.extend(found_files)
                else:
                    fastq_to_align.extend(found_files)
            elif input_type == 'gvcf':
                gvcfs.extend(found_files)
            elif input_type == 'cram':
                crams.extend(found_files)

    return gvcfs, crams, fastq_to_align, cram_to_realign


def _find_file_indices(fpaths: List[str]) -> Dict[str, Tuple[str, str]]:
    """
    Find corresponding index files. Will look for:
    '<sample>.g.vcf.gz.tbi' for '<sample>.g.vcf.gz',
    '<sample>.bam.bai' or '<sample>.bai' for '<sample>.bam',
    '<sample>.cram.crai' or '<sample>.crai' for '<sample>.cram',

    Returns a dict mapping sample name to a pair of file paths (file, index)
    """
    result: Dict[str, Tuple[str, str]] = dict()
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
        result[_get_file_base_name(fp)] = fp, index
    return result


def _add_ped_info(
    df: pd.DataFrame,
    ped_fpath: str,
    local_tmp_dir: str,
) -> pd.DataFrame:
    """
    Compares found input files to the samples provided in the PED file.
    Creates an input DataFrame using the PED file as a base.
    """
    local_sample_list_fpath = join(local_tmp_dir, basename(ped_fpath))
    subprocess.run(
        f'gsutil cat {ped_fpath} | grep -v ^Family.ID | cut -f1-6 > {local_sample_list_fpath}',
        check=False,
        shell=True,
    )
    ped_df = pd.read_csv(
        local_sample_list_fpath,
        delimiter='\t',
        names=[
            'Family.ID',
            's',
            'Paternal.ID',
            'Maternal.ID',
            'Sex',
            'Phenotype',
        ],
    ).set_index('s', drop=False)

    snames = set(df['s'])
    ped_snames = set(ped_df['s'])

    # Checking that file base names have a 1-to-1 match with the samples in PED
    # First, checking the match of PED sample names to input files
    if ped_snames - snames:
        logging.warning(
            f'No files found for the PED samples: {", ".join(ped_snames - snames)}'
        )
    if snames - ped_snames:
        logging.warning(
            f'No PED records found for the inputs: {", ".join(snames - ped_snames)}'
        )
    if ped_snames != snames:
        logging.warning(
            'Mismatches found between provided inputs and the provided PED file. '
            'Make sure all input files are named after sample names in the PED file: '
            '{sample}.cram, {sample}.g.vcf.gz, etc.'
        )

    for sn in ped_snames:
        if sn in snames:
            df.loc[sn, 'Family.ID'] = ped_df.loc[sn, 'Family.ID']
            df.loc[sn, 'Paternal.ID'] = ped_df.loc[sn, 'Paternal.ID']
            df.loc[sn, 'Maternal.ID'] = ped_df.loc[sn, 'Maternal.ID']
            df.loc[sn, 'Sex'] = ped_df.loc[sn, 'Sex']
            df.loc[sn, 'Phenotype'] = ped_df.loc[sn, 'Phenotype']

    return df


def _get_file_base_name(file_path):
    """
    Strips directory path and extention from a file name. Effectively extractss
    assumed sample name (i.e. assumes file named as gs://path/to/{sample}.g.vcf.gz)
    """
    return re.sub('(.bam|.cram|.g.vcf.gz)$', '', os.path.basename(file_path))


def splitext_gz(fname: str) -> Tuple[str, str]:
    """
    Split on file extensions, allowing for zipped extensions.
    """
    base, ext = splitext(fname)
    if ext in ['.gz', '.bz2', '.zip']:
        base, ext2 = splitext(base)
        ext = ext2 + ext
    return base, ext


def _find_fastq_pairs(fpaths: List[str]) -> Dict[str, Tuple[str, str]]:
    """
    Find pairs of FASTQ files
    """
    logger.info('Finding FASTQ pairs...')

    fastqs_by_sample_name: Dict[str, Tuple[Optional[str], Optional[str]]] = dict()
    for fpath in fpaths:
        fn, ext = splitext_gz(basename(fpath))
        if ext in ['.fq', '.fq.gz', '.fastq', '.fastq.gz']:
            sname, l_fpath, r_fpath = None, None, None
            if fn.endswith('_1'):
                sname = fn[:-2]
                l_fpath = fpath
            elif fn.endswith('_R1'):
                sname = fn[:-3]
                l_fpath = fpath
            elif fn.endswith('_2'):
                sname = fn[:-2]
                r_fpath = fpath
            elif fn.endswith('_R2'):
                sname = fn[:-3]
                r_fpath = fpath

            if sname:
                m = re.match(r'(.*)_S\d+', sname)
                if m:
                    sname = m.group(1)
                sname = sname.replace('-', '_')
            else:
                sname = fn
                logger.info('Cannot detect file for ' + sname)

            l, r = fastqs_by_sample_name.get(sname, (None, None))
            if l and l_fpath:
                logger.critical(
                    'Duplicated left FASTQ files for '
                    + sname
                    + ': '
                    + l
                    + ' and '
                    + l_fpath
                )
            if r and r_fpath:
                logger.critical(
                    'Duplicated right FASTQ files for '
                    + sname
                    + ': '
                    + r
                    + ' and '
                    + r_fpath
                )
            fastqs_by_sample_name[sname] = l or l_fpath, r or r_fpath

    paired_fastqs_by_sample_name: Dict[str, Tuple[str, str]] = dict()
    for sname, (l, r) in fastqs_by_sample_name.items():
        if not l:
            logger.error(f'ERROR: for sample {sname} left reads not found')
        if not r:
            logger.error(f'ERROR: for sample {sname} right reads not found')
        if l and r:
            paired_fastqs_by_sample_name[sname] = str(l), str(r)

    return paired_fastqs_by_sample_name
