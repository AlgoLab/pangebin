"""PangeBin-flow once MILP objectives."""

import gurobipy as gp

import pangebin.plasbin.classbin.milp.variables as lp_vars
import pangebin.plasbin.milp.objectives as pb_lp_obj
import pangebin.plasbin.milp.variables as pb_lp_var
import pangebin.plasbin.network as net


def plasmidness_score(
    network: net.Network,
    flow_vars: pb_lp_var.Flow,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Get linear expression for coverage score."""
    frag_set_fn = pb_lp_obj.ObjectiveFunctionDomain.to_fn(obj_fun_domain)
    return gp.quicksum(
        net.length(network, frag_id)
        * network.plasmidness(frag_id)
        * flow_vars.incoming_forward_reverse(network, frag_id)
        for frag_id in frag_set_fn(network)
    )


def classify_objective(
    var: lp_vars.Classify,
    network: net.Network,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> gp.LinExpr:
    """Set Classify objective."""
    return plasmidness_score(network, var.flow(), obj_fun_domain)


def set_classify_objective(
    m: gp.Model,
    var: lp_vars.Classify,
    network: net.Network,
    obj_fun_domain: pb_lp_obj.ObjectiveFunctionDomain,
) -> None:
    """Set MGCLB objective."""
    m.setObjective(
        classify_objective(var, network, obj_fun_domain),
        gp.GRB.MAXIMIZE,
    )


# REFACTOR potentially useless
# class MultiFlowObjectiveValues:
#     """Multi-flow objective values."""

#     def __init__(
#         self,
#         circular_obj_value: float,
#         partially_circular_obj_value: float,
#     ) -> None:
#         self.__circular_obj_value = circular_obj_value
#         self.__partially_circular_obj_value = partially_circular_obj_value

#     def circular_obj_value(self) -> float:
#         """Get objective value for circular bins."""
#         return self.__circular_obj_value

#     def set_circular_obj_value(self, circular_obj_value: float) -> None:
#         """Set objective value for circular bins."""
#         self.__circular_obj_value = circular_obj_value

#     def partially_circular_obj_value(self) -> float:
#         """Get objective value for partially circular bins."""
#         return self.__partially_circular_obj_value

#     def set_partially_circular_obj_value(
#         self,
#         partially_circular_obj_value: float,
#     ) -> None:
#         """Set objective value for partially circular bins."""
#         self.__partially_circular_obj_value = partially_circular_obj_value
