"""
Let this be the entrypoint / driver for loading data into SEQR for the CPG
See the README for more information. This is WIP.

    - 2021/04/16 Michael Franklin
"""

from batch.load import ElasticSearchCredentials

from typing import Optional


def main(
    source_paths: str,
    dest_path: str,
    genome_version: str,
    reference_path: str,
    clinvar_path: str,
    hgmd_path: Optional[str] = None,
    disable_validation: bool = False,
    dataset_type: str = "VARIANTS",
    remap_path: str = None,
    subset_path: str = None,
    vep_config_json_path: Optional[str] = None,
):
    """
    Args:
        source_paths: Path or list of paths of VCFs to be loaded.
        dest_path: Path to write the matrix table.
        genome_version: Reference Genome Version (37 or 38)
        reference_path: Path to the Hail table storing the reference variants.'
        clinvar_path: Path to the Hail table storing the clinvar variants.
        hgmd_path: [OPTIONAL] Path to the Hail table storing the hgmd variants.
        disable_validation: Disable checking whether the dataset matches the specified genome version and WGS vs. WES sample type.
        dataset_type: VARIANTS or SV.
        remap_path: Path to a tsv file with two columns: s and seqr_id.
        subset_path: Path to a tsv file with one column of sample IDs: s.
        vep_config_json_path: Path of hail vep config .json file

    Returns:

    """
    pass