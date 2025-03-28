"""Pan-assembly create module."""

from typing import TYPE_CHECKING

import gfapy  # type: ignore[import-untyped]

import pangebin.assembly.items as asm_items
import pangebin.gfa.link as gfa_link
import pangebin.gfa.panassembly.link as gfa_pan_link
import pangebin.gfa.panassembly.path as gfa_pan_path
import pangebin.gfa.panassembly.segment as gfa_pan_segment
import pangebin.gfa.path as gfa_path
import pangebin.gfa.segment as gfa_segment

if TYPE_CHECKING:
    from gfapy.line.group.path.path import (  # type: ignore[import-untyped]
        Path as GfaPath,
    )

import logging

_LOGGER = logging.getLogger(__name__)


def pangenome_graph_into_panassembly_graph(
    pangenome_gfa: gfapy.Gfa,
    standardized_skesa_gfa: gfapy.Gfa,
    standardized_unicycler_gfa: gfapy.Gfa,
) -> None:
    """Transform a pangenome graph into a panassembly graph.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph
    standardized_skesa_gfa : gfapy.Gfa
        SKESA GFA standardized graph
    standardized_unicycler_gfa : gfapy.Gfa
        Unicycler GFA standardized graph

    Warnings
    --------
    This function mutates the GFA graph.

    """
    # Paths

    set_path_lengths(pangenome_gfa, standardized_skesa_gfa, standardized_unicycler_gfa)
    set_path_normalized_coverages(
        pangenome_gfa,
        standardized_skesa_gfa,
        standardized_unicycler_gfa,
    )

    # Segments

    set_segment_lengths(pangenome_gfa)
    set_segment_contig_lists_and_occurences(pangenome_gfa)

    remove_unused_segments(pangenome_gfa)

    set_segment_contig_percentages(pangenome_gfa)
    set_segment_natures(pangenome_gfa)
    rename_segments(pangenome_gfa)

    # Links

    set_pangenome_link_origins(pangenome_gfa)

    for assembler_gfa, assmbler_id in (
        (standardized_skesa_gfa, asm_items.Identifier.SKESA),
        (standardized_unicycler_gfa, asm_items.Identifier.UNICYCLER),
    ):
        add_assembly_links_to_pangenome(pangenome_gfa, assembler_gfa, assmbler_id)

    set_link_types(pangenome_gfa)

    set_subcontig_normalized_converages(pangenome_gfa)
    set_subcontig_assembly_penalty(pangenome_gfa)

    pangenome_gfa.validate()


def set_path_lengths(
    pangenome_gfa: gfapy.Gfa,
    standardized_skesa_gfa: gfapy.Gfa,
    standardized_unicycler_gfa: gfapy.Gfa,
) -> None:
    """Set path lengths.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph
    standardized_skesa_gfa : gfapy.Gfa
        SKESA GFA standardized graph
    standardized_unicycler_gfa : gfapy.Gfa
        Unicycler GFA standardized graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting path lengths")

    for assembler_gfa in (standardized_skesa_gfa, standardized_unicycler_gfa):
        asm_segment_line: gfa_segment.GfaSegment
        for asm_segment_line in assembler_gfa.segments:
            pan_path_line = gfa_path.get_path_line_by_name(
                pangenome_gfa,
                asm_segment_line.name,
            )
            pan_path_line.set_datatype(
                gfa_pan_path.Tag.LENGTH,
                gfa_pan_path.TagType.LENGTH,
            )
            gfa_pan_path.set_length(pan_path_line, gfa_segment.length(asm_segment_line))


def set_path_normalized_coverages(
    pangenome_gfa: gfapy.Gfa,
    standardized_skesa_gfa: gfapy.Gfa,
    standardized_unicycler_gfa: gfapy.Gfa,
) -> None:
    """Set the normalized coverage tag to path in the pangenome GFA graph.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph
    standardized_skesa_gfa : gfapy.Gfa
        SKESA GFA standardized graph
    standardized_unicycler_gfa : gfapy.Gfa
        Unicycler GFA standardized graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting path normalized coverages")

    for assembler_gfa in (standardized_skesa_gfa, standardized_unicycler_gfa):
        asm_segment_line: gfa_segment.GfaSegment
        for asm_segment_line in assembler_gfa.segments:
            try:
                pan_path_line = gfa_path.get_path_line_by_name(
                    pangenome_gfa,
                    asm_segment_line.name,
                )
            except ValueError as err:
                _LOGGER.exception(
                    "Could not find path for segment %s",
                    asm_segment_line.name,
                )
                raise err from err
            gfa_pan_path.set_normalized_coverage(
                pan_path_line,
                gfa_segment.normalized_coverage(asm_segment_line),
            )


def set_segment_lengths(pangenome_gfa: gfapy.Gfa) -> None:
    """Set segment lengths.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting segment lengths")

    segment_line: gfa_segment.GfaSegment
    for segment_line in pangenome_gfa.segments:
        segment_line.set_datatype(
            gfa_segment.Tag.LENGTH,
            gfa_segment.TagType.LENGTH,
        )
        gfa_segment.set_length(segment_line)


def set_segment_contig_lists_and_occurences(pangenome_gfa: gfapy.Gfa) -> None:
    """Set segment contig list and occurences.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting segment contig lists and occurences")

    segment_line: gfa_segment.GfaSegment
    for segment_line in pangenome_gfa.segments:
        segment_line.set_datatype(
            gfa_pan_segment.Tag.CONTIG_LIST,
            gfa_pan_segment.TagType.FROM_CONTIGS,
        )
        gfa_pan_segment.set_contig_list(segment_line, [])

    path_line: GfaPath
    for path_line in pangenome_gfa.paths:
        oriented_identifier: gfapy.OrientedLine
        for oriented_identifier in path_line.segment_names:
            segment_line = pangenome_gfa.segment(
                oriented_identifier.name,
            )
            gfa_pan_segment.append_to_contig_list(segment_line, path_line.name)

    for segment_line in pangenome_gfa.segments:
        segment_line.set_datatype(
            gfa_pan_segment.Tag.OCCURENCE_IN_PATHS,
            gfa_pan_segment.TagType.OCCURENCE_IN_PATH,
        )
        gfa_pan_segment.set_occurence_in_paths(
            segment_line,
            len(gfa_pan_segment.contig_list(segment_line)),
        )


def remove_unused_segments(pangenome_gfa: gfapy.Gfa) -> None:
    """Clean pangenome graph.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Removing unused segments")

    _number_of_disconnected_segments = 0

    segment_line: gfa_segment.GfaSegment
    for segment_line in pangenome_gfa.segments:
        if (
            gfa_segment.length(segment_line) < 1
            or gfa_pan_segment.occurence_in_paths(segment_line) == 0
        ):
            _number_of_disconnected_segments += 1
            segment_line.disconnect()

    _LOGGER.debug(
        "Number of disconnected segments in the pangenome graph: %i",
        _number_of_disconnected_segments,
    )


def set_segment_contig_percentages(pangenome_gfa: gfapy.Gfa) -> None:
    """Set segment contig percentages.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting segment contig percentages")

    segment_line: gfa_segment.GfaSegment
    for segment_line in pangenome_gfa.segments:
        segment_line.set_datatype(
            gfa_pan_segment.Tag.CONTIG_PERCENTAGES,
            gfa_pan_segment.TagType.CONTIG_PERCENTAGES,
        )
        gfa_pan_segment.set_contig_percentages(segment_line, [])

    path_line: GfaPath
    for path_line in pangenome_gfa.paths:
        for oriented_identifier in path_line.segment_names:
            segment_line = pangenome_gfa.segment(oriented_identifier.name)
            segment_perc: float = gfa_segment.length(segment_line) / float(
                gfa_pan_path.length(path_line),
            )
            gfa_pan_segment.append_contig_percentage(
                segment_line,
                round(segment_perc, 6),
            )


def set_segment_natures(pangenome_gfa: gfapy.Gfa) -> None:
    """Set segment nature tags.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting segment natures")

    segment_line: gfa_segment.GfaSegment
    for segment_line in pangenome_gfa.segments:
        segment_line.set_datatype(
            gfa_pan_segment.Tag.SEGMENT_NATURE,
            gfa_pan_segment.TagType.SEGMENT_NATURE,
        )
        distinct_nature_tag_values: set[gfa_pan_segment.NatureTagValue] = {
            gfa_pan_segment.NatureTagValue.from_contig_name(contig_name)
            for contig_name in gfa_pan_segment.contig_list(segment_line)
        }

        if distinct_nature_tag_values:
            if distinct_nature_tag_values == {
                gfa_pan_segment.NatureTagValue.SUB_SKESA_CONTIG,
                gfa_pan_segment.NatureTagValue.SUB_UNICYCLER_CONTIG,
            }:
                gfa_pan_segment.set_segment_nature(
                    segment_line,
                    gfa_pan_segment.NatureTagValue.SUB_SKESA_AND_UNICYCLER_CONTIGS,
                )
            elif len(distinct_nature_tag_values) == 1:
                # XXX tmp condition to check there is no share between only contigs from the same assembler # noqa: E501
                if gfa_pan_segment.occurence_in_paths(segment_line) > 1:
                    _LOGGER.exception(
                        "Segment is shared between contigs only from the assembler: %s",
                        segment_line,
                    )
                    raise ValueError
                gfa_pan_segment.set_segment_nature(
                    segment_line,
                    distinct_nature_tag_values.pop(),
                )
            else:
                _LOGGER.exception(
                    "Unexpected distinct from assembler tag values: %s",
                    distinct_nature_tag_values,
                )
                raise ValueError
        else:
            _LOGGER.exception(
                "Unexpected empty distinct from assembler tag values for segment: %s",
                segment_line,
            )
            raise ValueError


def rename_segments(pangenome_gfa: gfapy.Gfa) -> None:
    """Rename pangenome segments.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Renaming segments")

    sub_skesa_suffix = 1
    sub_unicycler_suffix = 1
    both_suffix = 1

    segment_line: gfa_segment.GfaSegment
    for segment_line in pangenome_gfa.segments:
        segment_nature = gfa_pan_segment.segment_nature(segment_line)
        prefix = gfa_pan_segment.NamePrefix.from_segment_nature(segment_nature)
        match segment_nature:
            case gfa_pan_segment.NatureTagValue.SUB_SKESA_CONTIG:
                segment_line.name = gfa_segment.format_name(prefix, sub_skesa_suffix)
                sub_skesa_suffix += 1
            case gfa_pan_segment.NatureTagValue.SUB_UNICYCLER_CONTIG:
                segment_line.name = gfa_segment.format_name(
                    prefix,
                    sub_unicycler_suffix,
                )
                sub_unicycler_suffix += 1
            case gfa_pan_segment.NatureTagValue.SUB_SKESA_AND_UNICYCLER_CONTIGS:
                segment_line.name = gfa_segment.format_name(prefix, both_suffix)
                both_suffix += 1


def set_pangenome_link_origins(pangenome_gfa: gfapy.Gfa) -> None:
    """Set pangenome link origin tags.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting pangenome link origins")

    link_line: gfa_link.GfaLink
    for link_line in pangenome_gfa.dovetails:
        link_line.set_datatype(
            gfa_pan_link.Tag.LINK_ORIGIN,
            gfa_pan_link.TagType.LINK_ORIGIN,
        )
        gfa_pan_link.set_link_origin(
            link_line,
            gfa_pan_link.LinkOriginTagValue.PANGENOME_LINK,
        )


def add_assembly_links_to_pangenome(
    pangenome_gfa: gfapy.Gfa,
    assembler_gfa: gfapy.Gfa,
    assembler_id: asm_items.Identifier,
) -> None:
    """Add assembly links to pangenome graph.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph
    assembler_gfa : gfapy.Gfa
        Assembler GFA graph
    assembler_id : asm_items.Identifier
        Assembler ID

    Warnings
    --------
    This function mutates the GFA graph

    Note
    ----
    Two link lines in the GFA cannot have the same (or reverse-equivalent) definition.

    """
    _LOGGER.info("Adding %s assembly links to pangenome", assembler_id)

    for assembler_link in (
        gfa_link.Link.from_link_line(link_line) for link_line in assembler_gfa.dovetails
    ):
        assembler_link.simplify()

        try:
            path_from_contig = gfa_path.get_path_line_by_name(
                pangenome_gfa,
                assembler_link.predecessor().identifier(),
            )
            path_to_contig = gfa_path.get_path_line_by_name(
                pangenome_gfa,
                assembler_link.successor().identifier(),
            )
        except ValueError as err:
            _LOGGER.exception("Should continue the loop?")
            raise err from err

        if assembler_link.predecessor().is_forward():
            pred_frag = gfa_segment.OrientedFragment.from_segment_line_to_forward(
                path_from_contig.segment_names[-1],
            )
        else:
            pred_frag = gfa_segment.OrientedFragment.from_segment_line_to_reverse(
                path_from_contig.segment_names[0],
            )

        if assembler_link.successor().is_forward():
            succ_frag = gfa_segment.OrientedFragment.from_segment_line_to_forward(
                path_to_contig.segment_names[0],
            )
        else:
            succ_frag = gfa_segment.OrientedFragment.from_segment_line_to_reverse(
                path_to_contig.segment_names[-1],
            )

        new_panassembly_link = gfa_link.Link(pred_frag, succ_frag)
        new_panassembly_link_origin = gfa_pan_link.LinkOriginTagValue.from_assembler(
            assembler_id,
        )
        add_new_panassembly_link = True
        #
        # Convert multi-edges to single edge
        #
        # XXX See also unexpected multi-edges in GFA issue https://github.com/ggonnella/gfapy/issues/35
        _existing_link_with_orient = gfa_link.link_line_or_its_reversed_from_link(
            pangenome_gfa,
            new_panassembly_link,
        )
        if _existing_link_with_orient is not None:
            if (
                new_panassembly_link_origin
                != gfa_pan_link.LinkOriginTagValue.PANGENOME_LINK
            ):
                _LOGGER.debug("Previous multi-edge: %s", _existing_link_with_orient[0])
                _existing_link_with_orient[0].disconnect()
            else:
                _LOGGER.debug("New multi-edge: %s", new_panassembly_link)
                add_new_panassembly_link = False

        if add_new_panassembly_link:
            new_panassembly_link_line = new_panassembly_link.to_gfa_line()
            new_panassembly_link_line.set_datatype(
                gfa_pan_link.Tag.LINK_ORIGIN,
                gfa_pan_link.TagType.LINK_ORIGIN,
            )
            gfa_pan_link.set_link_origin(
                new_panassembly_link_line,
                new_panassembly_link_origin,
            )
            pangenome_gfa.add_line(new_panassembly_link_line)


def set_link_types(pangenome_gfa: gfapy.Gfa) -> None:
    """Set link type tags.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting link types")

    link_line: gfa_link.GfaLink
    for link_line in pangenome_gfa.dovetails:
        link_line.set_datatype(
            gfa_pan_link.Tag.LINK_TYPE,
            gfa_pan_link.TagType.LINK_TYPE,
        )

        gfa_pan_link.set_link_type(
            link_line,
            gfa_pan_segment.NatureTagValue.from_name(
                link_line.from_name,
            ),
            gfa_pan_segment.NatureTagValue.from_name(
                link_line.to_name,
            ),
        )


def set_subcontig_normalized_converages(pangenome_gfa: gfapy.Gfa) -> None:
    """Set subcontig normalized coverages.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting subcontig normalized coverages")

    segment_line: gfa_segment.GfaSegment
    for segment_line in pangenome_gfa.segments:
        segment_line.set_datatype(
            gfa_segment.Tag.NORMALIZED_COVERAGE,
            gfa_segment.TagType.NORMALIZED_COVERAGE,
        )
        coverage_list: list[float] = [
            gfa_pan_path.normalized_coverage(
                gfa_path.get_path_line_by_name(
                    pangenome_gfa,
                    contig_name,
                ),
            )
            for contig_name in gfa_pan_segment.contig_list(segment_line)
        ]

        coverage_mean = sum(coverage_list) / len(coverage_list)
        gfa_segment.set_normalized_coverage(segment_line, coverage_mean)


def set_subcontig_assembly_penalty(pangenome_gfa: gfapy.Gfa) -> None:
    """Set assembly penalty.

    Parameters
    ----------
    pangenome_gfa : gfapy.Gfa
        Pangenome GFA graph

    Warnings
    --------
    This function mutates the GFA graph

    """
    _LOGGER.info("Setting assembly penalty")

    segment_line: gfa_segment.GfaSegment
    for segment_line in pangenome_gfa.segments:
        segment_line.set_datatype(
            gfa_pan_segment.Tag.PANGENOME_PENALTY,
            gfa_pan_segment.TagType.PANGENOME_PENALTY,
        )
        penalty: float
        match gfa_pan_segment.segment_nature(segment_line):
            case gfa_pan_segment.NatureTagValue.SUB_SKESA_AND_UNICYCLER_CONTIGS:
                penalty = 0.0
            case _:
                penalty = min(1, gfa_segment.length(segment_line) / 1000)
        gfa_pan_segment.set_pangenome_penalty(segment_line, penalty)
