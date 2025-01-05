"""Utilities to manipulate assemblis and pangenomes in GFA format."""

from __future__ import annotations

import sys

import gfapy  # type: ignore[import-untyped]
from gfapy.line.segment import Segment as GfaSegment  # type: ignore[import-untyped]


def get_path_by_name(gfa: gfapy.Gfa, name):
    if gfa.line(str(name)) is not None:
        return gfa.line(str(name))
    print("path", name, "not found!!!!", type(name))
    return None


def get_edge_by_def(gfa: gfapy.Gfa, def_: list):
    for edge in gfa.dovetails:
        if def_[0] == edge.from_segment.name and def_[2] == edge.to_segment.name:
            if edge.from_orient == def_[1] and edge.to_orient == def_[3]:
                return edge
    return None


def add_gfa_to_pangenome(gfa: gfapy.Gfa, pangenome: gfapy.Gfa):
    _type = gfa.segments[0].name[0]

    for seg in gfa.segments:
        match = pangenome.line(seg.name)
        if match is None:
            seg.disconnect()
        else:
            assert match.aa == _type
            if _type in ("u", "s"):
                dp = seg.dp
                assert seg.dp is not None
                match.set_datatype("dp", "f")
                match.dp = float(dp)
            else:
                _err_msg = f"unknow assemler type: {_type}"
                raise TypeError(_err_msg)

    for edge in gfa.dovetails:
        from_contig = edge.from_segment
        to_contig = edge.to_segment
        from_orient = edge.from_orient
        to_orient = edge.to_orient

        try:
            path_from_contig = get_path_by_name(pangenome, from_contig.name)
            path_to_contig = get_path_by_name(pangenome, to_contig.name)
        except:
            print("path not found")
            sys.exit(1)

        if path_from_contig is None or path_to_contig is None:
            continue

        if (from_orient == "+") and (to_orient == "+"):
            left_frag = [x.name for x in path_from_contig.segment_names][-1]
            left_orient = "+"
            right_frag = next(x.name for x in path_to_contig.segment_names)
            right_orient = "+"

        elif (from_orient == "+") and (to_orient == "-"):
            left_frag = [x.name for x in path_from_contig.segment_names][-1]
            left_orient = "+"
            right_frag = [x.name for x in path_to_contig.segment_names][-1]
            right_orient = "-"

        elif (from_orient == "-") and (to_orient == "+"):
            left_frag = next(x.name for x in path_from_contig.segment_names)
            left_orient = "-"
            right_frag = next(x.name for x in path_to_contig.segment_names)
            right_orient = "+"

        elif (from_orient == "-") and (to_orient == "-"):
            # swapping left-right and using + +
            left_frag = [x.name for x in path_to_contig.segment_names][-1]
            left_orient = "+"
            right_frag = next(x.name for x in path_from_contig.segment_names)
            right_orient = "+"

        new_edge = gfapy.Line(
            f"L\t{left_frag}\t{left_orient}\t{right_frag}\t{right_orient}\t0M\taa:A:{_type}\tlt:Z:{_type}{_type}",
        )

        try:
            edge_to_update = get_edge_by_def(
                pangenome,
                [left_frag, left_orient, right_frag, right_orient],
            )
            if edge_to_update is None:
                pangenome.add_line(new_edge)
            else:
                pass
                # docu: print("edge", edge_to_update,"already present")
                # docu: edge_to_update.set_datatype("ls", "a")
        except gfapy.error.notuniqueerror:
            pass


def compute_scores(pangenome: gfapy.Gfa):
    for edge in pangenome.dovetails:
        if pangenome.segment(edge.from_segment.name).name == "1":
            print(pangenome.segment(edge.from_segment.name))
        if edge.lt is None:
            edge.set_datatype("aa", "A")
            edge.aa = "p"
            edge.set_datatype("lt", "Z")
            inc = pangenome.segment(edge.from_segment.name).aa
            out = pangenome.segment(edge.to_segment.name).aa
            if inc is None:
                inc = "n"
            if out is None:
                out = "n"

            edge.lt = inc + out
        inc = pangenome.segment(edge.from_segment.name).aa
        out = pangenome.segment(edge.to_segment.name).aa
        if (
            pangenome.segment(edge.from_segment.name).OC == 1
            and float(pangenome.segment(edge.from_segment.name).ll.split(",")[0]) == 1.0
        ):
            inc = inc.upper()
        if (
            pangenome.segment(edge.to_segment.name).OC == 1
            and float(pangenome.segment(edge.to_segment.name).ll.split(",")[0]) == 1.0
        ):
            out = out.upper()
        edge.lt = inc + out

    for path in pangenome.paths:
        path.set_datatype("cv", "f")  # normalized coverage
        try:
            assert path.dp is not None
            path.cv = path.dp
        except:
            print("in path: path", path.name, "has no coverage:", path)
            sys.exit(1)

    # annotate mean coverage to fragments, based on paths (contigs) coverage
    seg: GfaSegment
    for seg in pangenome.segments:
        seg.set_datatype("cv", "f")
        coverage_list = []

        if seg.cl is None:
            _err_msg = f"segment {seg} has no paths associated"
            raise ValueError(_err_msg)

        for contig in seg.cl.split(","):
            path = pangenome.line(contig)
            assert path is not None
            if contig == path.name:
                try:
                    assert path.cv is not None
                except:
                    print("in seg: path", path.name, "has no coverage:", path)
                    sys.exit(1)
                coverage_list.append(path.cv)
            else:
                _err_msg = f"segment {seg} has no paths associated"
                raise KeyError(_err_msg)

        coverage_mean = sum(coverage_list) / len(coverage_list)
        assert coverage_mean is not None
        seg.cv = float(coverage_mean)

    ## annotate assembly penalty to segments

    for seg in pangenome.segments:
        seg.set_datatype("ap", "f")
        if seg.aa in ("u", "s"):
            penalty_value = seg.LN / 1000
            penalty_value = min(penalty_value, 1)
            seg.ap = penalty_value
        else:
            seg.ap = 0


def clean_pangenome(gfa: gfapy.Gfa):
    gfa.validate()  # REFACTOR GFA validation here
    seg: GfaSegment
    for seg in gfa.segments:
        seg.set_datatype("OC", "i")  # occurences in paths
        seg.OC = 0
        seg.set_datatype("cl", "Z")  # segment list
        seg.cl = ""
        seg.set_datatype("ll", "Z")  # length of segment list
        seg.ll = ""
        seg.LN = len(seg.sequence)
        seg.set_datatype("aa", "A")  # assembly type

    for path in gfa.paths:
        segs = path.segment_names
        _type = path.name[0]
        path_len = 0
        for seg in segs:
            segment = gfa.segment(seg.name)
            segment.OC += 1
            path_len += segment.LN
            if segment.cl == "":
                segment.cl = path.name
            else:
                segment.cl += "," + path.name

        path.set_datatype("LN", "i")
        path.LN = path_len

        path.set_datatype("aa", "A")
        path.aa = _type

    for path in gfa.paths:
        segs = path.segment_names
        for seg in segs:
            segment = gfa.segment(seg.name)
            segment_perc = segment.LN / float(path.LN)
            if segment.ll == "":
                segment.ll = f"{segment_perc:.6f}"
            else:
                segment.ll += f",{segment_perc:.6f}"

    for seg in gfa.segments:
        if seg.LN < 1:
            seg.disconnect()
            continue
        if seg.OC == 0:
            seg.disconnect()
            continue

        contig_list = seg.cl.split(",")
        contig_list = list(set(contig_list))
        sorted_clist = sorted(contig_list)
        if len(sorted_clist) >= 1:
            if sorted_clist[0] != "" and sorted_clist[0][0] != sorted_clist[-1][0]:
                seg.aa = "b"
            elif sorted_clist[0] != "":
                seg.aa = sorted_clist[0][0]

    gfa.validate()
    return gfa
