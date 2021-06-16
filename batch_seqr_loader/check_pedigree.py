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
    '--somalier-prefix',
    'somalier_prefix',
    required=True,
    help='Prefix to somalier output files, with required {prefix}.pairs.tsv',
)
@click.option('--local-tmp-dir', 'local_tmp_dir')
def main(
    somalier_prefix: str,
    local_tmp_dir: str,
):  # pylint: disable=missing-function-docstring

    somalier_pairs_tsv_path = f'{somalier_prefix}.pairs.tsv'
    somalier_samples_tsv_path = f'{somalier_prefix}.samples.tsv'
    somalier_html_path = f'{somalier_prefix}.html'

    if somalier_prefix.startswith('gs://'):
        local_tmp_dir = local_tmp_dir or tempfile.mkdtemp()
        subprocess.run(
            f'gsutil cp {somalier_samples_tsv_path} {somalier_pairs_tsv_path} {local_tmp_dir}/',
            check=False,
            shell=True,
        )
        somalier_pairs_tsv_path = join(local_tmp_dir, basename(somalier_pairs_tsv_path))
        somalier_samples_tsv_path = join(
            local_tmp_dir, basename(somalier_samples_tsv_path)
        )

    # Checking sex
    df = pd.read_csv(somalier_samples_tsv_path, delimiter='\t')
    mismatching_sex = (df['sex'] == 2) & (df['original_pedigree_sex'] != 'female') | (
        df['sex'] == 1
    ) & (df['original_pedigree_sex'] != 'male')

    if mismatching_sex.any():
        logger.error(
            f'Found samples with mismatched sex: {df[mismatching_sex].sample_id}. '
            f'Review the somalier results for more detail: {somalier_html_path}'
        )
        sys.exit(1)
    else:
        logger.info(
            f'All sample sex match the provided one. '
            f'Review the somalier results for more detail: {somalier_html_path}'
        )

    # Checking relatedness
    related_pairs = set()  # pylint: disable=unused-variable
    ped = Ped(somalier_samples_tsv_path)  # pylint: disable=unused-variable
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
    df = pd.read_csv(somalier_pairs_tsv_path, delimiter='\t')
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
            f'Review the somalier results for more detail: {somalier_html_path}'
        )
        sys.exit(1)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
