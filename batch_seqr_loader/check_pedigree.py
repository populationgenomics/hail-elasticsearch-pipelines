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

logger = logging.getLogger()


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
        logger.error(
            f'Found samples with mismatched sex: {df[mismatching_sex].sample_id}. '
            f'Review the somalier results for more detail: {somalier_html_fpath}'
        )
        sys.exit(1)
    else:
        logger.info(
            f'All sample sex match the provided one. '
            f'Review the somalier results for more detail: {somalier_html_fpath}'
        )

    # Checking relatedness
    related_pairs = set()  # pylint: disable=unused-variable
    ped = Ped(somalier_samples_fpath)  # pylint: disable=unused-variable
    for s1, s2 in combinations(ped.samples(), 2):
        relation = ped.relation(s1, s2)
        if relation not in [
            'unrelated',
            'unknown',
            'related at unknown level',
            'mom-dad',
        ]:
            related_pairs.add(tuple(sorted([s1.sample_id, s2.sample_id])))

    mismatching_pairs = set()
    df = pd.read_csv(somalier_pairs_fpath, delimiter='\t')
    for row in df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        rel = row['relatedness']
        exp_rel = row['expected_relatedness']
        if exp_rel > 0.2 and rel < 0.1:
            mismatching_pairs.add((s1, s2, rel, exp_rel))
    if mismatching_pairs:
        logger.error(
            f'Found sample pairs with mismatched relatedness: {mismatching_pairs}. '
            f'Review the somalier results for more detail: {somalier_html_fpath}'
        )
        sys.exit(1)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
