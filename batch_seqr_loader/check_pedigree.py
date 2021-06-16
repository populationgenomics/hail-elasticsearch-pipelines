#!/usr/bin/env python3

"""
This script parses `somalier relate` (https://github.com/brentp/somalier) outputs,
and returns a non-zero code if either sex or pedigree is mismatching
"""

import logging
import subprocess
import sys
import tempfile
from itertools import combinations
from os.path import join, basename
import pandas as pd
import click
from peddy import Ped

logger = logging.getLogger('check-pedigree')
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--somalier-samples',
    'somalier_samples_fpath',
    required=True,
    help='Path to somalier {prefix}.samples.tsv output file',
)
@click.option(
    '--somalier-pairs',
    'somalier_pairs_fpath',
    required=True,
    help='Path to somalier {prefix}.pairs.tsv output file',
)
@click.option(
    '--somalier-html',
    'somalier_html_fpath',
    help='Path to somalier {prefix}.html output file',
)
def main(
    somalier_samples_fpath: str,
    somalier_pairs_fpath: str,
    somalier_html_fpath: str,
):  # pylint: disable=missing-function-docstring
    if somalier_samples_fpath.startswith('gs://'):
        local_tmp_dir = tempfile.mkdtemp()
        subprocess.run(
            f'gsutil cp {somalier_samples_fpath} {somalier_pairs_fpath} {local_tmp_dir}/',
            check=False,
            shell=True,
        )
        somalier_pairs_fpath = join(local_tmp_dir, basename(somalier_pairs_fpath))
        somalier_samples_fpath = join(local_tmp_dir, basename(somalier_samples_fpath))
    
    # Checking sex
    df = pd.read_csv(somalier_samples_fpath, delimiter='\t')
    mismatching_sex = (df['sex'] == 2) & (df['original_pedigree_sex'] != 'female') | (
        df['sex'] == 1
    ) & (df['original_pedigree_sex'] != 'male')

    if mismatching_sex.any():
        logger.info(
            f'Found samples with mismatched sex: {df[mismatching_sex].sample_id}. ' +
            (f'Review the somalier results for more detail: {somalier_html_fpath}'
             if somalier_html_fpath else '')
        )

    # Checking relatedness
    ped = Ped(somalier_samples_fpath)
    sample_by_id = {s.sample_id: s for s in ped.samples()}
    mismatching_pairs = []
    df = pd.read_csv(somalier_pairs_fpath, delimiter='\t')
    for _, row in df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        inferred_rel = row['relatedness']
        provided_rel = row['expected_relatedness']
        if (provided_rel > 0.2 and inferred_rel < 0.1 or
             inferred_rel > 0.2 and provided_rel < 0.1):
            mismatching_pairs.append(
                f'{s1}, {s2}, '
                f'provided relatedness: {provided_rel} '
                f'({ped.relation(sample_by_id[s1], sample_by_id[s2])}), '
                f'inferred: {inferred_rel}')
    if mismatching_pairs:
        logger.info(f'Found sample pairs with mismatched relatedness:')
        for pair in mismatching_pairs:
            logger.info(pair)
        if somalier_html_fpath:
            logger.info(f'Review the somalier results for more detail: {somalier_html_fpath}')
    
    if mismatching_sex.any() or mismatching_pairs:
        sys.exit(1)
    
    logger.info(
        f'\nInferred sex and pedigree matches for all samples with the data in the PED file. ' +
        (f'Review the somalier results for more detail: {somalier_html_fpath}'
         if somalier_html_fpath else '')
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
