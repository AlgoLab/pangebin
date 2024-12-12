"""Utilities to manipulate assemblis and pangenomes in GFA format."""

import contextlib
import subprocess
import sys
from itertools import product
from pathlib import Path

import gfapy
from Bio import SeqIO as Seq


def invert(sign):
    if sign == "+":
        return "-"
    if sign == "-":
        return "+"
    if sign == "l":
        return "r"
    if sign == "r":
        return "l"
    return "?"


def edge_exists(edge, collection):
    for e in collection:
        if edge_compare(edge, e) or edge_compare(edge_reverse(edge), e):
            return True
    return True


def edge_compare(edge1, edge2):
    return bool(
        str(edge1[0][0]) == str(edge2[0][0])
        and str(edge1[1][0]) == str(edge2[1][0])
        and str(edge1[0][1]) == str(edge2[0][1])
        and str(edge1[1][1]) == str(edge2[1][1]),
    )


def edge_reverse(edge):
    assert edge[0][2] == "l", "edge[0][2] should be left, 'l'"
    assert edge[1][2] == "r", "edge[1][2] should be right, 'r'"

    return (
        (edge[1][0], invert(edge[1][1]), "l"),
        (edge[0][0], invert(edge[0][1]), "r"),
    )


def stoe(string):
    estr = str(string).split("\t")
    return ((estr[1], estr[2], "l"), (estr[3], estr[4], "r"))


def etos(edge):
    l_orient = "?"
    r_orient = "?"

    if edge[0][2] == "r":
        l_orient = invert(edge[0][1])
    elif edge[0][2] == "l":
        l_orient = edge[0][1]

    if edge[1][2] == "l":
        r_orient = invert(edge[1][1])
    elif edge[1][2] == "r":
        r_orient = edge[1][1]

    return gfapy.Line(f"L\t{edge[0][0]}\t{l_orient}\t{edge[1][0]}\t{r_orient}\t0M")


def extract_node(edge, orient, name):
    node = stoe(edge)

    if orient == "r":
        if node[0][0] == name and node[1][0] == name:
            return None
        if node[1][1] == "-" and node[1][0] == name:
            node = edge_reverse(node)
        elif node[0][1] == "+" and node[0][0] == name:
            assert node[0][0] == name
        else:
            assert False
        return (node[1][0], node[1][1], orient)

    if orient == "l":
        if node[0][0] == name and node[1][0] == name:
            return None
        if node[1][1] == "+" and node[1][0] == name:
            pass
        elif node[0][1] == "-" and node[0][0] == name:
            node = edge_reverse(node)
        else:
            assert False
        return (node[0][0], node[0][1], orient)
    return None


def rename_contigs(
    gfa,
    prefix,
) -> gfapy.gfa:
    graph = gfa
    # graph = gfapy.gfa.Gfa.from_file(gfa, vlevel=0)
    graph.validate()
    graph.vlevel = 3
    counter = 1
    for seg in graph.segments:
        if seg.LN is not None:
            if seg.LN == 0:
                seg.sequence = "*"
                seg.LN = 0
                graph.validate()
        elif len(seg.sequence) == 0:
            seg.sequence = "*"
            seg.LN = 0
            graph.validate()
        seg.name = f"{prefix}{counter}"
        counter += 1
    graph.validate()
    # graph.to_file(gfa)
    return graph


def convert_kc_to_dp(
    gfa,
) -> gfapy.gfa:
    graph = gfa
    graph.validate()
    total_coverage = 0
    total_length = 0
    for seg in graph.segments:
        total_coverage += seg.KC
        if seg.LN is None:
            seg.set_datatype("LN", "i")
            seg.LN = len(seg.sequence)
        total_length += seg.LN

    for seg in graph.segments:
        seg.set_datatype("dp", "f")
        seg.dp = float((seg.KC * total_length) / (seg.LN * total_coverage))
        seg.KC = None
    graph.validate()
    return graph


def remove_nodes_fasta(fasta, threshold, output):
    fasta_input = fasta
    fasta_dict = dict()
    with Path.open(fasta_input) as spe:
        for record in Seq.parse(spe, "fasta"):
            if len(record) >= threshold:
                fasta_dict[record.id] = record.seq
    with Path.open(output, "w") as out:
        for x, rec in fasta_dict.items():
            out.write(f">{x}\n")
            out.write(f"{rec!s}\n")


def remove_nodes(gfa, threshold, output):
    for seg in gfa.segments:
        if (seg.LN is not None and threshold >= seg.LN) or (
            len(seg.sequence) <= threshold
        ):
            pass
        else:
            continue
        right_edges = list(seg.dovetails_R)
        left_edges = list(seg.dovetails_L)
        right_nodes = set()  # edge(node, orient, 'r')
        left_nodes = set()

        for e in right_edges:
            n = extract_node(e, "r", seg.name)
            if n is not None:
                right_nodes.add(n)
                gfa.rm(e)
        for e in left_edges:
            n = extract_node(e, "l", seg.name)
            if n is not None:
                left_nodes.add(n)
                gfa.rm(e)

        gfa.rm(seg)
        pairs = list(product(left_nodes, right_nodes))
        gfa.validate()
        for edge in pairs:
            new_edge = etos(edge)
            with contextlib.suppress(gfapy.error.NotUniqueError):
                gfa.add_line(new_edge)
        gfa.validate()
    gfa.to_file(output)
    return output


def gfa_to_fasta(gfa, fasta):
    command = f'awk \'/^S/\u007bprint ">"$2; print ""$3\u007d\' {gfa} > {fasta}'
    subprocess.call(command, shell=True)
    return fasta


def mix_fasta(fasta_list, output):
    fasta1 = fasta_list[0]
    fasta2 = fasta_list[1]
    command = f"cat {fasta1} {fasta2} > {output}"
    subprocess.call(command, shell=True)


def get_path_by_name(gfa: gfapy.gfa, name):
    if gfa.line(str(name)) is not None:
        return gfa.line(str(name))
    print("path", name, "not found!!!!", type(name))
    return None


def get_segment_by_name(gfa: gfapy.gfa, name):
    if gfa.segment(name) is not None:
        return gfa.segment(name)
    return None


def get_edge_by_def(gfa: gfapy.gfa, def_: list):
    for edge in gfa.dovetails:
        if def_[0] == edge.from_segment.name and def_[2] == edge.to_segment.name:
            if edge.from_orient == def_[1] and edge.to_orient == def_[3]:
                return edge
    return None


def add_gfa_to_pangenome(gfa, pangenome):
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
            right_frag = [x.name for x in path_to_contig.segment_names][0]
            right_orient = "+"

        elif (from_orient == "+") and (to_orient == "-"):
            left_frag = [x.name for x in path_from_contig.segment_names][-1]
            left_orient = "+"
            right_frag = [x.name for x in path_to_contig.segment_names][-1]
            right_orient = "-"

        elif (from_orient == "-") and (to_orient == "+"):
            left_frag = [x.name for x in path_from_contig.segment_names][0]
            left_orient = "-"
            right_frag = [x.name for x in path_to_contig.segment_names][0]
            right_orient = "+"

        elif (from_orient == "-") and (to_orient == "-"):
            # swapping left-right and using + +
            left_frag = [x.name for x in path_to_contig.segment_names][-1]
            left_orient = "+"
            right_frag = [x.name for x in path_from_contig.segment_names][0]
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


def compute_scores(pangenome):
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

        coverage_mean = sum([i for i in coverage_list]) / len(coverage_list)
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


def clean_pangenome(gfa):
    gfa.validate()

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


def extract_gfagz(gfagz) -> gfapy.Gfa:
    command = f"bgzip -dkf {gfagz}"

    subprocess.call(command, shell=True)
    filename = str(gfagz)[:-3]
    gfa = gfapy.Gfa.from_file(filename, vlevel=0)
    gfa.validate()
    return gfa
