"""MILP model setup module."""

from gurobipy import GRB, LinExpr, Model

from pangebin.plasbin.data_utils import (
    COV_KEY,
    DEFAULT_HEAD_STR,
    DEFAULT_SINK,
    DEFAULT_SOURCE,
    DEFAULT_TAIL_STR,
    SEED_KEY,
)


# -----------------------------------------------------------
# INITIALIZING VARIABLES
# -----------------------------------------------------------
# contigs[c] == 1 if contig 'c' belongs to solution, else 0
def contig_vars(m: Model, contigs_dict, contigs):
    contigs = {}
    for c in contigs_dict:
        contigs[c] = m.addVar(vtype=GRB.BINARY, name="contig-" + c)
    return contigs


# links[e] == 1 if link 'e' belongs to solution, else 0
# For each link e, we also add the edge variable in the reverse direction
def link_vars(m: Model, links_list, links, contigs):
    links = {}
    for e in links_list:
        links[e] = m.addVar(
            vtype=GRB.BINARY,
            name="link-" + e[0][0] + e[0][1] + e[1][0] + e[1][1],
        )
        links[e[::-1]] = m.addVar(
            vtype=GRB.BINARY,
            name="link-" + e[1][0] + e[1][1] + e[0][0] + e[0][1],
        )

    for c in contigs:
        h_ext = (c, DEFAULT_HEAD_STR)
        t_ext = (c, DEFAULT_TAIL_STR)

        links[(DEFAULT_SOURCE, h_ext)] = m.addVar(
            vtype=GRB.BINARY,
            name="link-S-to-" + c + "-h",
        )
        links[(h_ext, DEFAULT_SINK)] = m.addVar(
            vtype=GRB.BINARY,
            name="link-" + c + "-h-to-T",
        )
        links[(DEFAULT_SOURCE, t_ext)] = m.addVar(
            vtype=GRB.BINARY,
            name="link-S-to-" + c + "-t",
        )
        links[(t_ext, DEFAULT_SINK)] = m.addVar(
            vtype=GRB.BINARY,
            name="link-" + c + "-t-to-T",
        )
    return links


def GC_vars(m: Model, gc_probs, plas_GC, contig_GC):
    for c in gc_probs:
        contig_GC[c] = {}
        for b in gc_probs[c]:
            contig_GC[c][b] = m.addVar(
                vtype=GRB.BINARY,
                name="contig-" + c + "-GC-bin-" + str(b),
            )
            if b not in plas_GC:
                plas_GC[b] = m.addVar(vtype=GRB.BINARY, name="plasmid-GC-bin-" + str(b))
    return plas_GC, contig_GC


def seed_vars(m: Model, contigs_dict, counted_seed):
    for c in contigs_dict:
        counted_seed[c] = m.addVar(vtype=GRB.BINARY, name="counted-seed_contig-" + c)
    return counted_seed


def len_vars(m: Model, contigs_dict, counted_len):
    counted_ln = {}
    for c in contigs_dict:
        counted_ln[c] = m.addVar(vtype=GRB.INTEGER, name="counted-ln_contig-" + c)
    return counted_ln


def flow_vars(m: Model, links, flows, counted_overall_flow):
    for e in links:
        flows[e] = m.addVar(
            vtype=GRB.CONTINUOUS,
            name="flow_edge-" + str(e[0]) + "-to-" + str(e[1]),
        )
        counted_overall_flow[e] = m.addVar(
            vtype=GRB.CONTINUOUS,
            name="counted-overall-flow_edge-" + str(e[0]) + "-to-" + str(e[1]),
        )
    return flows, counted_overall_flow


# -----------------------------------------------------------
# SETTING CONSTRAINTS
# -----------------------------------------------------------
def link_inclusion_constr(m: Model, links, contigs, constraint_count):
    # print("contigs:", contigs)
    for e in links:
        # print("link:", e)
        end1, end2 = e[0], e[1]
        # print("end1:", end1)
        # print("end2:", end2)
        c1, c2 = end1[0], end2[0]
        # print("c1:", c1)
        # print("c2:", c2)

        ## problem with 't' SINK, check if it is a problem with the data or the code
        # if c1 != DEFAULT_SOURCE and (c1 != DEFAULT_SINK):
        if c1 != DEFAULT_SOURCE and not (c1 == DEFAULT_SINK or c1 == DEFAULT_TAIL_STR):
            m.addConstr(
                links[e] <= contigs[c1],
                "link_ubd-" + str(e) + "-by-" + str(c1),
            )
            constraint_count += 1
        # if c2 != DEFAULT_SOURCE and (c2 != DEFAULT_SINK):
        if c2 != DEFAULT_SOURCE and not (c2 == DEFAULT_SINK or c2 == DEFAULT_TAIL_STR):
            m.addConstr(
                links[e] <= contigs[c2],
                "link_ubd-" + str(e) + "-by-" + str(c2),
            )
            constraint_count += 1
    return m, constraint_count


def extr_inclusion_constr(m: Model, links, contigs, extr_dict, constraint_count):
    for c in contigs:
        x1, x2 = (c, DEFAULT_HEAD_STR), (c, DEFAULT_TAIL_STR)
        expr = LinExpr()
        link_count = 0
        for link in extr_dict[x1]:
            link_count += 1
            expr.addTerms(1, links[link])
            for link in extr_dict[x2]:
                link_count += 1
                expr.addTerms(1, links[link])
        m.addConstr(contigs[c] <= expr, "extr_ubd-" + str(c) + "-by-expr")
        m.addConstr(contigs[c] <= 1, "extr_ubd-" + str(c) + "-by-1")
        constraint_count += 2
        if link_count > 0:
            m.addConstr(
                contigs[c] >= expr / link_count,
                "extr_lbd-" + str(c) + "-by-expr",
            )
            constraint_count += 1
    return m, constraint_count


def seed_inclusion_constr(m: Model, contigs, contigs_dict, constraint_count):
    seed_expr = LinExpr()
    for c in contigs:
        seed_expr.addTerms(contigs_dict[c][SEED_KEY], contigs[c])
    m.addConstr(seed_expr >= 1, "at-least-one-seed")
    constraint_count += 1
    return m, constraint_count


def min_flow_constr(m: Model, links, flows, F, LBD_rd, constraint_count):
    xS_expr = LinExpr()
    xT_expr = LinExpr()
    fS_expr = LinExpr()
    fT_expr = LinExpr()
    for e in links:
        if e[0] == DEFAULT_SOURCE:
            xS_expr.addTerms(1, links[e])
            fS_expr.addTerms(1, flows[e])
        if e[1] == DEFAULT_SINK:
            xT_expr.addTerms(1, links[e])
            fT_expr.addTerms(1, flows[e])
    m.addConstr(xS_expr == 1, "one-edge-from-S")
    m.addConstr(xT_expr == 1, "one-edge-into-T")
    m.addConstr(fS_expr == F, "flow-from-S")
    m.addConstr(fT_expr == F, "flow-into-T")
    constraint_count += 4
    m.addConstr(LBD_rd <= F, "minimum flow (to avoid zero flow)")
    return m, constraint_count


def flow_conservation_constraints(
    m: Model,
    links,
    contigs,
    flows,
    incoming,
    outgoing,
    capacities,
    contigs_dict,
    constraint_count,
):
    for c in contigs:
        ext1, ext2 = (c, DEFAULT_HEAD_STR), (c, DEFAULT_TAIL_STR)
        hin_flow = LinExpr()  # flow into head extremity
        hout_flow = LinExpr()  # flow out of head extremity
        tin_flow = LinExpr()  # flow into tail extremity
        tout_flow = LinExpr()  # flow out of tail extremity

        # TODO: incoming and outgoing should be updated with graph update to avoid 'if' conditions
        if ext1 in incoming:
            for e in incoming[ext1]:
                hin_flow.addTerms(1, flows[e])
        if ext1 in outgoing:
            for e in outgoing[ext1]:
                hout_flow.addTerms(1, flows[e])
        if ext2 in incoming:
            for e in incoming[ext2]:
                tin_flow.addTerms(1, flows[e])
        if ext2 in outgoing:
            for e in outgoing[ext2]:
                tout_flow.addTerms(1, flows[e])

        m.addConstr(hin_flow == tout_flow, "contig-" + c + "-flow-conservation-h2t")
        m.addConstr(tin_flow == hout_flow, "contig-" + c + "-flow-conservation-t2h")
        m.addConstr(hin_flow + tin_flow <= contigs_dict[c][COV_KEY], "cap-" + c)
        constraint_count += 3

    for e in links:
        m.addConstr(flows[e] <= capacities[e] * links[e], "cap-" + str(e))
        constraint_count += 1

    return m, constraint_count


def counted_flow_constr(m: Model, links, flows, counted_F, F, UBD_rd, constraint_count):
    for e in links:
        m.addConstr(counted_F[e] <= UBD_rd * links[e], "xF1-" + str(e))
        m.addConstr(counted_F[e] <= F, "xF2-" + str(e))
        m.addConstr(counted_F[e] >= F - (1 - links[e]) * UBD_rd, "xF3-" + str(e))
        m.addConstr(counted_F[e] >= 0, "xF4-" + str(e))
        m.addConstr(counted_F[e] <= flows[e], "xF5-" + str(e))
        constraint_count += 5

    return m, constraint_count


def GC_constr(m: Model, contig_GC, plas_GC, contigs, constraint_count):
    for c in contig_GC:
        for b in plas_GC:
            m.addConstr(contig_GC[c][b] <= contigs[c], "xGC1")
            m.addConstr(contig_GC[c][b] <= plas_GC[b], "xGC2")
            m.addConstr(contig_GC[c][b] >= contigs[c] + plas_GC[b] - 1, "xGC3")
            constraint_count += 3
    expr = LinExpr()
    for b in plas_GC:
        expr.addTerms(1, plas_GC[b])
    m.addConstr(expr == 1)
    constraint_count += 1

    return m, constraint_count


def test_constr(m: Model, a, b, j, k, n):
    m.addConstr(j * a + k * b <= n, "constr_" + str(j) + str(k) + str(n))
    return m
