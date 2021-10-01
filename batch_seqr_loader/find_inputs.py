"""
Provides find_inputs function that finds input files (.g.vcf.gz, .bam or .cram) 
in provided buckets, along with corresponding indices (tbi, bai, crai).
Compares sample names to the provided PED, and returns DataFrame.
"""

import logging
import os
import re
import subprocess
from dataclasses import dataclass
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
    gvcf_with_index_by_by_sn = find_file_indices(gvcfs)
    cram_with_index_by_by_sn = find_file_indices(crams)
    fastq_pairs_by_sn = find_fastq_pairs(fastq_to_align)
    cram_with_index_to_realign_by_sn = find_file_indices(cram_to_realign)

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
        | set(fastq_pairs_by_sn.keys())
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
        if sn in fastq_pairs_by_sn:
            file1, file2 = fastq_pairs_by_sn[sn]
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
            glob_path = glob_path.strip("'").strip('"')
            assert any(
                glob_path.endswith(ext) for ext in expected_ext_by_type[input_type]
            ), (input_type, glob_path)
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


def find_file_indices(fpaths: List[str]) -> Dict[str, Tuple[str, Optional[str]]]:
    """
    Find corresponding index files. Will look for:
    '<sample>.g.vcf.gz.tbi' for '<sample>.g.vcf.gz',
    '<sample>.bam.bai' or '<sample>.bai' for '<sample>.bam',
    '<sample>.cram.crai' or '<sample>.crai' for '<sample>.cram',

    Returns a dict mapping sample name to a pair of file paths (file, index).
    if it can't find an index file, it will add None, which will trigger a re-index job
    """
    result: Dict[str, Tuple[str, Optional[str]]] = dict()
    for fp in fpaths:
        index: Optional[str] = fp + '.tbi'
        if fp.endswith('.g.vcf.gz'):
            if not file_exists(index):
                logger.critical(f'Not found TBI index for file {fp}, will create')
                index = None
        elif fp.endswith('.bam'):
            index = fp + '.bai'
            if not file_exists(index):
                index = re.sub('.bam$', '.bai', fp)
                if not file_exists(index):
                    logger.warning(f'Not found BAI index for file {fp}, will create')
                    index = None
        elif fp.endswith('.cram'):
            index = fp + '.crai'
            if not file_exists(index):
                index = re.sub('.cram$', '.crai', fp)
                if not file_exists(index):
                    logger.warning(f'Not found CRAI index for file {fp}, will create')
                    index = None
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
        f'gsutil cat {ped_fpath} | grep -v ^Family.ID | cut -f1-6 > '
        f'{local_sample_list_fpath}',
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


@dataclass
class AlignmentInput:
    """
    Sort of a union type for possible alignment inputs
    """

    bam_or_cram_path: Optional[str] = None
    index_path: Optional[str] = None
    fqs1: Optional[List[str]] = None
    fqs2: Optional[List[str]] = None


def sm_get_reads_data(  # pylint: disable=too-many-return-statements
    meta: Dict,
    check_existence: bool = True,
) -> Optional[AlignmentInput]:
    """
    Verify the meta.reads object in a sample db entry
    """
    reads_data = meta.get('reads')
    reads_type = meta.get('reads_type')
    if not reads_data:
        logger.error(f'No "meta/reads" field in meta')
        return None
    if not reads_type:
        logger.error(f'No "meta/reads_type" field in meta')
        return None
    supported_types = ('fastq', 'bam', 'cram')
    if reads_type not in supported_types:
        logger.error(f'ERROR: "reads_type" is expected to be one of {supported_types}')
        return None

    if reads_type in ('bam', 'cram'):
        if len(reads_data) > 1:
            logger.error('Supporting only single bam/cram input')
            return None

        bam_path = reads_data[0]['location']
        if not (bam_path.endswith('.cram') or bam_path.endswith('.bam')):
            logger.error(
                f'ERROR: expected the file to have an extention .cram or .bam,'
                f'got: {bam_path}'
            )
            return None
        if check_existence and not file_exists(bam_path):
            logger.error(f'ERROR: index file doesn\'t exist: {bam_path}')
            return None

        # Index:
        if not reads_data[0].get('secondaryFiles'):
            logger.error(
                f'ERROR: bam/cram input is expected to have '
                f'a non-empty list field "secondaryFile" section with indices'
            )
            return None
        index_path = reads_data[0]['secondaryFiles'][0]['location']
        if (
            bam_path.endswith('.cram')
            and not index_path.endswith('.crai')
            or bam_path.endswith('.bai')
            and not index_path.endswith('.bai')
        ):
            logger.error(
                f'ERROR: expected the index file to have an extention '
                f'.crai or .bai, got: {index_path}'
            )
        if check_existence and not file_exists(index_path):
            logger.error(f'ERROR: index file doesn\'t exist: {index_path}')
            return None

        return AlignmentInput(bam_or_cram_path=bam_path, index_path=index_path)

    else:
        fqs1 = []
        fqs2 = []
        for lane_data in reads_data:
            assert len(lane_data) == 2, lane_data
            if check_existence and not file_exists(lane_data[0]['location']):
                logger.error(
                    f'ERROR: read 1 file doesn\'t exist: {lane_data[0]["location"]}'
                )
                return None
            if check_existence and not file_exists(lane_data[1]['location']):
                logger.error(
                    f'ERROR: read 2 file doesn\'t exist: {lane_data[1]["location"]}'
                )
                return None

            fqs1.append(lane_data[0]['location'])
            fqs2.append(lane_data[1]['location'])
        return AlignmentInput(fqs1=fqs1, fqs2=fqs2)


def find_fastq_pairs(fpaths: List[str]) -> Dict[str, Tuple[str, str]]:
    """
    Find pairs of FASTQ files for each sample
    """
    logger.info('Finding FASTQ pairs...')

    fastqs_by_sample_name: Dict[str, Tuple[List[str], List[str]]] = dict()
    for fpath in fpaths:
        fn, ext = splitext_gz(basename(fpath))
        if ext not in ['.fq', '.fq.gz', '.fastq', '.fastq.gz']:
            continue
        sname, l_fpath, r_fpath = None, None, None

        # Parsing the file name according to the Illumina spec
        # https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source
        # /Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
        # Example: SampleName_S1_L001_R1_001.fastq.gz

        # Stripping the segment number
        m = re.match(r'(.*)_\d\d\d$', fn)
        if m:
            fn = m.group(1)

        # Parsing the number in pair
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
        else:
            logger.critical(
                f'Fastq base file name is expected to have a _1/_2/_R1/_R2 '
                f'suffix. Found: {fn}'
            )

        if sname:
            # Stripping different combinations of S_ (sample number) and L_ (lane)
            for suf in [r'_L\d+', r'_S\d+']:
                m = re.match(r'(.*)' + suf, sname)
                if m:
                    sname = m.group(1)
                sname = sname.replace('-', '_')
        else:
            sname = fn
            logger.info('Cannot detect file for ' + sname)

        if sname not in fastqs_by_sample_name:
            fastqs_by_sample_name[sname] = ([], [])
        ls, rs = fastqs_by_sample_name[sname]

        if l_fpath:
            if ls and l_fpath in ls:
                logger.critical(
                    'Duplicated left FASTQ files for ' + sname + ': ' + l_fpath
                )
            ls.append(str(l_fpath))

        if r_fpath:
            if rs and r_fpath in rs:
                logger.critical(
                    'Duplicated right FASTQ files for ' + sname + ': ' + r_fpath
                )
            rs.append(str(r_fpath))
        fastqs_by_sample_name[sname] = ls, rs

    for sname, (ls, rs) in fastqs_by_sample_name.items():
        if len(ls) == 0:
            logger.error(f'ERROR: for sample {sname} left reads not found')
        if len(rs) == 0:
            logger.error(f'ERROR: for sample {sname} right reads not found')
        if len(ls) != len(rs):
            logger.error(
                f'ERROR: for sample {sname}, the number of '
                f'left fastqs ({len(ls)}) != the number of '
                f'right fastqs ({len(rs)})'
            )

    # Joining muiltiple pairs with comma
    return {
        sn: (','.join(rs), ','.join(ls))
        for sn, (rs, ls) in fastqs_by_sample_name.items()
    }
