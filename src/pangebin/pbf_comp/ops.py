"""PlasBin-flow conversion operations."""

import pangebin.pbf_comp.items as pbf_items
import pangebin.plasbin.bins.input_output as bin_io
import pangebin.plasbin.bins.items as bin_item
import pangebin.plasbin.input_output as pb_io


def pbf_to_pg_plasmidness(pbf_plasmidness: float) -> float:
    """Convert PBF plasmidness to PG plasmidness."""
    return 2 * pbf_plasmidness - 1


def pg_bin_dir_to_pbf_bininfo(
    io_manager: pb_io.Manager,
    bin_number: int,
) -> pbf_items.PBFBinInfo:
    """Convert PangeBin bin directory to PlasBin-flow bin info."""
    with bin_io.Reader.open(io_manager.bin_seq_normcov_path(bin_number)) as bin_fin:
        seq_mults = [pbf_items.ContigMult(seq_id, mult) for seq_id, mult in bin_fin]
    bin_stats = bin_item.Stats.from_yaml(io_manager.bin_stats_path(bin_number))
    return pbf_items.PBFBinInfo(
        f"P{bin_number + 1}",
        bin_stats.coverage_flow(),
        bin_stats.gc_content_interval(),
        seq_mults,
    )
