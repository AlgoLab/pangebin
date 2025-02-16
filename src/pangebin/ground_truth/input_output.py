"""Ground truth input-output module."""

from collections.abc import Iterator
from pathlib import Path

import yaml  # type: ignore[import-untyped]

import pangebin.ground_truth.items as gt_items


def plasmid_biosamples_from_file(
    yaml_file: Path,
) -> Iterator[gt_items.PlasmidBioSampleNCBI]:
    """Get plasmid biosamples from file.

    Parameters
    ----------
    yaml_file : Path
        Path to the YAML file

    Yield
    -----
    PlasmidBioSampleNCBI
        Plasmid BioSample from NCBI
    """
    with yaml_file.open() as f_in:
        for dict_from_yaml in yaml.safe_load_all(f_in):
            yield gt_items.PlasmidBioSampleNCBI.from_dict(dict_from_yaml)
