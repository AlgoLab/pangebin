"""GFA operations."""

from __future__ import annotations

import contextlib
import itertools
import logging
from typing import TYPE_CHECKING, cast

import gfapy  # type: ignore[import-untyped]
import networkx as nx

import pangebin.gfa.header as gfa_header
import pangebin.gfa.line as gfa_line
import pangebin.gfa.segment as gfa_segment
import pangebin.gfa.tag as gfa_tag
import pangebin.input_output as io

from . import link as gfa_link

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path

    from gfapy.line.edge import (  # type: ignore[import-untyped]
        Containment as GfaContainment,
    )
    from gfapy.line.edge import Link as GfaLink
    from gfapy.line.group.path.path import (  # type: ignore[import-untyped]
        Path as GfaPath,
    )
    from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]


_LOGGER = logging.getLogger(__name__)


def is_skesa_gfa_fixed(gfa_path: Path) -> bool:
    """Check ig the SKESA GFA file is fixed."""
    yes_fix_tag = (
        f"{gfa_header.Tag.SKESA_FIX}"
        f":{gfa_header.TagType.SKESA_FIX}"
        f":{gfa_header.SKESAFixTagValue.YES}"
    )
    with io.open_file_read(gfa_path) as f_in:
        for line in f_in:
            if (
                line.startswith(str(gfa_line.Type.HEADER))
                and yes_fix_tag in line.split()
            ):
                _LOGGER.debug("SKESA GFA file is fixed.")
                return True
    _LOGGER.debug("SKESA GFA file is not fixed.")
    return False


def is_unicycler_gfa_fixed(gfa_path: Path) -> bool:
    """Check ig the Unicycler GFA file is fixed."""
    yes_fix_tag = (
        f"{gfa_header.Tag.UNICYCLER_FIX}"
        f":{gfa_header.TagType.UNICYCLER_FIX}"
        f":{gfa_header.UnicyclerFixTagValue.YES}"
    )
    with io.open_file_read(gfa_path) as f_in:
        for line in f_in:
            if (
                line.startswith(str(gfa_line.Type.HEADER))
                and yes_fix_tag in line.split()
            ):
                _LOGGER.debug("Unicycler GFA file is fixed.")
                return True
    _LOGGER.debug("Unicycler GFA file is not fixed.")
    return False


def fix_skesa_gfa(
    in_gfa_path: Path,
    out_gfa_path: Path | None = None,
) -> Path:
    """Fix a SKESA GFA file.

    Parameters
    ----------
    in_gfa_path : Path
        Path to input GFA file
    out_gfa_path : Path | None, optional
        Path to output GFA file, must be different from input if used, by default None

    Returns
    -------
    Path
        Path to output GFA file

    Warnings
    --------
    This function modifies the input GFA file if out_gfa_path is None.

    Raises
    ------
    ValueError
        If input and output files are the same

    """
    if out_gfa_path is None:
        _LOGGER.debug("Fixing SKESA GFA file %s", in_gfa_path)
    else:
        _LOGGER.debug("Fixing SKESA GFA file %s to %s.", in_gfa_path, out_gfa_path)

    yes_fix_tag = (
        f"{gfa_header.Tag.SKESA_FIX}"
        f":{gfa_header.TagType.SKESA_FIX}"
        f":{gfa_header.SKESAFixTagValue.YES}"
    )
    with (
        io.possible_tmp_file(in_gfa_path, out_gfa_path) as (use_in_path, use_out_path),
        io.open_file_read(use_in_path) as f_in,
        io.open_file_write(use_out_path) as f_out,
    ):
        f_out.write(f"{gfa_line.Type.HEADER}\t{yes_fix_tag}\n")
        for line in f_in:
            if line.startswith(gfa_line.Type.SEGMENT):
                split_line = line.split()
                f_out.write(split_line[0])  # S
                f_out.write("\t")
                f_out.write(split_line[1])  # segment name
                f_out.write("\t")
                f_out.write(split_line[2])  # sequence
                for optional_field in split_line[3:]:
                    tag_name, tag_type, value = optional_field.split(":")
                    if tag_type == gfa_tag.FieldType.SIGNED_INT:
                        f_out.write(f"\t{tag_name}:{tag_type}:{int(float(value))}")
                    else:
                        f_out.write(f"\t{tag_name}:{tag_type}:{value}")
                f_out.write("\n")
            else:
                f_out.write(line)

    return use_out_path


def fix_unicycler_gfa(
    in_gfa_path: Path,
    out_gfa_path: Path | None = None,
) -> Path:
    """Fix a Unicycler GFA file.

    Because of https://github.com/rrwick/Unicycler/issues/355

    Parameters
    ----------
    in_gfa_path : Path
        Path to input GFA file
    out_gfa_path : Path | None, optional
        Path to output GFA file, must be different from input if used, by default None

    Returns
    -------
    Path
        Path to output GFA file

    Warnings
    --------
    This function modifies the input GFA file if out_gfa_path is None.

    Raises
    ------
    ValueError
        If input and output files are the same

    """
    if out_gfa_path is None:
        _LOGGER.debug("Fixing Skeza GFA file %s", in_gfa_path)
    else:
        _LOGGER.debug("Fixing Skeza GFA file %s to %s.", in_gfa_path, out_gfa_path)

    entry_point_to_orient = {
        "-": "+",
        "+": "-",
    }

    yes_fix_tag = (
        f"{gfa_header.Tag.UNICYCLER_FIX}"
        f":{gfa_header.TagType.UNICYCLER_FIX}"
        f":{gfa_header.UnicyclerFixTagValue.YES}"
    )
    with (
        io.possible_tmp_file(in_gfa_path, out_gfa_path) as (use_in_path, use_out_path),
        io.open_file_read(use_in_path) as f_in,
        io.open_file_write(use_out_path) as f_out,
    ):
        f_out.write(f"{gfa_line.Type.HEADER}\t{yes_fix_tag}\n")
        for line in f_in:
            if line.startswith(gfa_line.Type.LINK):
                split_line = line.split()
                f_out.write(split_line[0])  # L
                f_out.write("\t")
                f_out.write(split_line[1])  # segment 1 name
                f_out.write("\t")
                f_out.write(split_line[2])  # from orient
                f_out.write("\t")
                f_out.write(split_line[3])  # segment 2 name
                f_out.write("\t")
                f_out.write(entry_point_to_orient[split_line[4]])  # to orient
                f_out.write("\t")
                f_out.write(split_line[5])  # Overlap
                for optional_field in split_line[6:]:
                    f_out.write(f"\t{optional_field}")
                f_out.write("\n")
            else:
                f_out.write(line)

    return use_out_path


def set_segment_length_tags(graph: gfapy.Gfa) -> None:
    """Set the segment length attribute in a GFA graph.

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting segment length tags.")
    segment: GfaSegment
    for segment in graph.segments:
        segment.set_datatype(
            gfa_segment.Tag.LENGTH,
            gfa_segment.TagType.LENGTH,
        )
        if len(segment.sequence) == 0:
            segment.sequence = gfapy.Placeholder()
            gfa_segment.set_length(segment, 0)
        else:
            gfa_segment.set_length(segment)


def rename_segments(
    graph: gfapy.Gfa,
    segment_prefix_name: str,
) -> None:
    """Rename segments in a GFA graph.

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph
    segment_prefix_name : str
        Prefix for segment names

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Renaming segments.")
    segment: GfaSegment
    for counter, segment in enumerate(graph.segments):
        segment.name = gfa_segment.format_name(segment_prefix_name, counter + 1)


def convert_kmer_coverage_to_normalized_coverage(
    graph: gfapy.Gfa,
) -> None:
    """Convert k-mer coverage to normalized coverage.

    Parameters
    ----------
    graph : gfapy.Gfa
        GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Converting k-mer coverage to normalized coverage.")
    total_coverage = sum(gfa_segment.kmer_coverage(seg) for seg in graph.segments)
    total_length = sum(gfa_segment.length(seg) for seg in graph.segments)

    seg: GfaSegment
    for seg in graph.segments:
        seg.set_datatype(
            gfa_segment.Tag.NORMALIZED_COVERAGE,
            gfa_segment.TagType.NORMALIZED_COVERAGE,
        )
        gfa_segment.set_normalized_coverage(
            seg,
            float(
                (gfa_segment.kmer_coverage(seg) * total_length)
                / (gfa_segment.length(seg) * total_coverage),
            ),
        )
        seg.delete(gfa_segment.Tag.KMER_COVERAGE)


def set_standardized_header_tag(gfa: gfapy.Gfa) -> None:
    """Set the standardized header tag in the GFA graph."""
    _LOGGER.info("Setting standardized header tag.")
    gfa.header.set_datatype(
        gfa_header.Tag.STANDARDIZED,
        gfa_header.TagType.STANDARDIZED,
    )
    gfa_header.set_standardized_header_tag(gfa, is_standardized=True)


def sub_radius_graph(
    gfa_graph: gfapy.Gfa,
    segment_names: Iterable[str],
    radius: int,
) -> gfapy.Gfa:
    """Get a subgraph of a GFA graph."""
    kept_segment_names = set()

    lifo_distanced_segments: list[tuple[str, int]] = [
        (segment_name, 0) for segment_name in segment_names
    ]
    while lifo_distanced_segments:
        segment_name, distance = lifo_distanced_segments.pop()
        if distance == radius:
            kept_segment_names.add(segment_name)
        elif distance < radius and segment_name not in kept_segment_names:
            # OPTIMIZE could do better (if only we have a better API...)
            kept_segment_names.add(segment_name)
            link_line: GfaLink
            for link_line in gfa_graph.segment(segment_name).dovetails:
                lifo_distanced_segments.extend(
                    (neighbor, distance + 1)
                    for neighbor in (
                        link_line.field_to_s("from_segment"),
                        link_line.field_to_s("to_segment"),
                    )
                    if neighbor not in kept_segment_names
                )
    sub_graph = gfapy.Gfa()
    for line in gfa_graph.lines:
        sub_graph.add_line(str(line))

    #
    # Remove segments
    #
    segment_to_remove: GfaSegment
    for segment_to_remove in [
        s for s in sub_graph.segments if s.name not in kept_segment_names
    ]:
        with contextlib.suppress(gfapy.RuntimeError):
            sub_graph.rm(segment_to_remove)
    #
    # Remove remaining links
    #
    links_to_remove: GfaLink
    for links_to_remove in [
        link_line
        for link_line in sub_graph.dovetails
        if cast("GfaLink", link_line).field_to_s("from_segment")
        not in kept_segment_names
        or cast("GfaLink", link_line).field_to_s("to_segment") not in kept_segment_names
    ]:
        with contextlib.suppress(gfapy.RuntimeError):
            sub_graph.rm(links_to_remove)

    #
    # Remove remaining paths
    #
    path_to_remove: GfaPath
    for path_to_remove in [
        path_line
        for path_line in sub_graph.paths
        if any(
            cast("gfapy.OrientedLine", oriented_identifier).name
            not in kept_segment_names
            for oriented_identifier in cast("GfaPath", path_line).segment_names
        )
    ]:
        with contextlib.suppress(gfapy.RuntimeError):
            sub_graph.rm(path_to_remove)

    return sub_graph


def transform_small_contigs_into_links(
    gfa: gfapy.Gfa,
    min_contig_length: int,
) -> Iterator[str]:
    """Remove small contigs and transform them into links.

    Parameters
    ----------
    gfa : gfapy.Gfa
        GFA graph
    min_contig_length : int
        Minimum contig length

    Warnings
    --------
    This function mutates the GFA graph

    Yields
    ------
    str
        Segment lines removed
    """
    segment_lines_to_remove: list[gfa_segment.GfaSegment] = []
    for segment in gfa.segments:
        if gfa_segment.length(segment) < min_contig_length:
            left_edge_lines = list(segment.dovetails_L)
            right_edge_lines = list(segment.dovetails_R)
            predecessors: list[gfa_segment.OrientedFragment] = []
            successors: list[gfa_segment.OrientedFragment] = []

            link_lines_to_remove: list[gfa_link.GfaLink] = []

            for left_link_line in left_edge_lines:
                pred = gfa_segment.OrientedFragment.from_left_dovetail_line(
                    left_link_line,
                    segment.name,
                )
                if pred.identifier() != segment.name:
                    predecessors.append(pred)
                    link_lines_to_remove.append(left_link_line)

            for right_link_line in right_edge_lines:
                succ = gfa_segment.OrientedFragment.from_right_dovetail_line(
                    right_link_line,
                    segment.name,
                )
                if succ.identifier() != segment.name:
                    successors.append(succ)
                    link_lines_to_remove.append(right_link_line)

            for link_line in link_lines_to_remove:
                gfa.rm(link_line)
            gfa.validate()

            segment_lines_to_remove.append(segment)

            for link in (
                gfa_link.Link(pred, succ)
                for pred, succ in itertools.product(predecessors, successors)
            ):
                with contextlib.suppress(gfapy.NotUniqueError):
                    # XXX the exception should not happen
                    gfa.add_line(link.to_gfa_line())
            gfa.validate()

    for segment in segment_lines_to_remove:
        gfa.rm(segment)

    gfa.validate()

    return (str(line) for line in segment_lines_to_remove)


def connected_components(graph: gfapy.Gfa) -> list[gfapy.Gfa]:
    """Get connected GFA (GFA1) graphs.

    Warning
    -------
    If one sequence is null, it is replaced by a placeholder (*) is the resulting GFA.
    """
    if graph.version != "gfa1":
        raise NotImplementedError

    # Build naive nx graph
    nx_graph: nx.Graph = nx.Graph()
    edge_line: GfaLink | GfaContainment
    for edge_line in graph.edges:
        nx_graph.add_edge(edge_line.from_segment, edge_line.to_segment)

    # Get connected components on nx graph
    d_seg_cc: dict[str, int] = {}
    num_cc = 0
    for k, component in enumerate(nx.connected_components(nx_graph)):
        num_cc += 1
        for segment_name in component:
            d_seg_cc[segment_name] = k

    connected_graphs: list[gfapy.Gfa] = [
        gfapy.Gfa(version="gfa1") for _ in range(num_cc)
    ]

    # Segment lines
    segment_line: GfaSegment
    for segment_line in graph.segments:
        if not segment_line.sequence:
            segment_line.sequence = gfapy.Placeholder()
        connected_graphs[d_seg_cc[segment_line.name]].add_line(str(segment_line))

    # Link and containment lines
    for edge_line in graph.edges:
        connected_graphs[d_seg_cc[edge_line.from_segment]].add_line(str(edge_line))

    # Path
    for path_line in graph.paths:
        connected_graphs[d_seg_cc[path_line.segment_names[0].name]].add_line(
            str(path_line),
        )

    return connected_graphs
