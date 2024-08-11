#!/usr/bin/env julia
###### AC-OPF using JuMP ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
# only the built-in AD library is supported
#

include("parse_case.jl")
import Ipopt
using JuMP
using AmplNLWriter
import MathOptSymbolicAD

function jump_opf_model(model, ref)
    JuMP.@variable(model, va[i in keys(ref[:bus])])
    JuMP.@variable(
        model,
        ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"],
        start = 1.0
    )

    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    JuMP.@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    JuMP.@variable(
        model,
        -ref[:branch][l]["rate_a"] <= p[(l, i, j) in ref[:arcs]] <= ref[:branch][l]["rate_a"]
    )
    JuMP.@variable(
        model,
        -ref[:branch][l]["rate_a"] <= q[(l, i, j) in ref[:arcs]] <= ref[:branch][l]["rate_a"]
    )

    JuMP.@objective(
        model,
        Min,
        sum(
            gen["cost"][1] * pg[i]^2 + gen["cost"][2] * pg[i] + gen["cost"][3] for
            (i, gen) in ref[:gen]
        )
    )

    for (i, bus) in ref[:ref_buses]
        JuMP.@constraint(model, 1.0 * va[i] == 0)
    end

    for (i, bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        JuMP.@constraint(
            model,
            sum(p[a] for a in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) - sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts) * vm[i]^2
        )

        JuMP.@constraint(
            model,
            sum(q[a] for a in ref[:bus_arcs][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) - sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts) * vm[i]^2
        )
    end

    # Branch power flow physics and limit constraints
    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        # From side of the branch flow
        JuMP.@constraint(
            model,
            p_fr ==
            (g + g_fr) / ttm * vm_fr^2 +
            (-g * tr + b * ti) / ttm * (vm_fr * vm_to * cos(va_fr - va_to)) +
            (-b * tr - g * ti) / ttm * (vm_fr * vm_to * sin(va_fr - va_to))
        )
        JuMP.@constraint(
            model,
            q_fr ==
            -(b + b_fr) / ttm * vm_fr^2 -
            (-b * tr - g * ti) / ttm * (vm_fr * vm_to * cos(va_fr - va_to)) +
            (-g * tr + b * ti) / ttm * (vm_fr * vm_to * sin(va_fr - va_to))
        )

        # To side of the branch flow
        JuMP.@constraint(
            model,
            p_to ==
            (g + g_to) * vm_to^2 +
            (-g * tr - b * ti) / ttm * (vm_to * vm_fr * cos(va_to - va_fr)) +
            (-b * tr + g * ti) / ttm * (vm_to * vm_fr * sin(va_to - va_fr))
        )
        JuMP.@constraint(
            model,
            q_to ==
            -(b + b_to) * vm_to^2 -
            (-b * tr + g * ti) / ttm * (vm_to * vm_fr * cos(va_to - va_fr)) +
            (-g * tr - b * ti) / ttm * (vm_to * vm_fr * sin(va_to - va_fr))
        )

        # Voltage angle difference limit
        JuMP.@constraint(model, branch["angmin"] <= va_fr - va_to <= branch["angmax"])

        # Apparent power limit, from side and to side
        JuMP.@constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        JuMP.@constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
    end
end

function jump_opf_model_rectangular(model, ref)
    # JuMP.@variable(model, va[i in keys(ref[:bus])])
    JuMP.@variable(
        model,
        -ref[:bus][i]["vmax"] <= vr[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"],
        start = 1.0
    )
    JuMP.@variable(
        model,
        -ref[:bus][i]["vmax"] <= vi[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"],
        start = 0.0
    )

    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    JuMP.@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    JuMP.@variable(
        model,
        -ref[:branch][l]["rate_a"] <= p[(l, i, j) in ref[:arcs]] <= ref[:branch][l]["rate_a"]
    )
    JuMP.@variable(
        model,
        -ref[:branch][l]["rate_a"] <= q[(l, i, j) in ref[:arcs]] <= ref[:branch][l]["rate_a"]
    )

    JuMP.@objective(
        model,
        Min,
        sum(
            gen["cost"][1] * pg[i]^2 + gen["cost"][2] * pg[i] + gen["cost"][3] for
            (i, gen) in ref[:gen]
        )
    )

    for (i, bus) in ref[:ref_buses]
        JuMP.@constraint(model, vi[i] == 0)
    end

    for (i, bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        JuMP.@constraint(model, bus["vmin"]^2 <= vr[i]^2 + vi[i]^2)
        JuMP.@constraint(model, bus["vmax"]^2 >= vr[i]^2 + vi[i]^2)

        JuMP.@constraint(
            model,
            sum(p[a] for a in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) - sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts) * (vr[i]^2 + vi[i]^2)
        )

        JuMP.@constraint(
            model,
            sum(q[a] for a in ref[:bus_arcs][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) - sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts) * (vr[i]^2 + vi[i]^2)
        )
    end

    # Branch power flow physics and limit constraints
    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]                     # p_fr is a reference to the optimization variable p[f_idx]
        q_fr = q[f_idx]                     # q_fr is a reference to the optimization variable q[f_idx]
        p_to = p[t_idx]                     # p_to is a reference to the optimization variable p[t_idx]
        q_to = q[t_idx]                     # q_to is a reference to the optimization variable q[t_idx]

        vr_fr = vr[branch["f_bus"]]         # vm_fr is a reference to the optimization variable vm on the from side of the branch
        vr_to = vr[branch["t_bus"]]         # vm_to is a reference to the optimization variable vm on the to side of the branch
        vi_fr = vi[branch["f_bus"]]         # va_fr is a reference to the optimization variable va on the from side of the branch
        vi_to = vi[branch["t_bus"]]         # va_fr is a reference to the optimization variable va on the to side of the branch

        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        # From side of the branch flow
        JuMP.@constraint(
            model,
            p_fr ==
            (g + g_fr) / ttm * (vr_fr^2 + vi_fr^2) +
            (-g * tr + b * ti) / ttm * (vr_fr * vr_to + vi_fr * vi_to) +
            (-b * tr - g * ti) / ttm * (vi_fr * vr_to - vr_fr * vi_to)
        )
        JuMP.@constraint(
            model,
            q_fr ==
            -(b + b_fr) / ttm * (vr_fr^2 + vi_fr^2) -
            (-b * tr - g * ti) / ttm * (vr_fr * vr_to + vi_fr * vi_to) +
            (-g * tr + b * ti) / ttm * (vi_fr * vr_to - vr_fr * vi_to)
        )

        # To side of the branch flow
        JuMP.@constraint(
            model,
            p_to ==
            (g + g_to) * (vr_to^2 + vi_to^2) +
            (-g * tr - b * ti) / ttm * (vr_fr * vr_to + vi_fr * vi_to) +
            (-b * tr + g * ti) / ttm * (-(vi_fr * vr_to - vr_fr * vi_to))
        )
        JuMP.@constraint(
            model,
            q_to ==
            -(b + b_to) * (vr_to^2 + vi_to^2) -
            (-b * tr + g * ti) / ttm * (vr_fr * vr_to + vi_fr * vi_to) +
            (-g * tr - b * ti) / ttm * (-(vi_fr * vr_to - vr_fr * vi_to))
        )

        # Voltage angle difference limit
        JuMP.@constraint(
            model,
            tan(branch["angmin"]) <=
            (vi_fr * vr_to - vr_fr * vi_to) / (vr_fr * vr_to + vi_fr * vi_to) <=
            tan(branch["angmax"])
        )

        # Apparent power limit, from side and to side
        JuMP.@constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        JuMP.@constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
    end
end

function set_ipopt_parameters(m)
    set_optimizer_attribute(m, "print_level", 5)
    set_optimizer_attribute(m, "linear_solver", "ma27")
    set_optimizer_attribute(m, "hsllib", "libhsl.dll")
    set_optimizer_attribute(m, "print_timing_statistics", "yes")
    set_optimizer_attribute(m, "max_iter", 200)
end

function jump_main(logdir, filename, method)
    logpath = joinpath(logdir, "$(filename)_jump.log")

    mpath = joinpath(@__DIR__, "case", "$(filename).m")
    ref = parse_case(mpath)

    m = JuMP.direct_model(Ipopt.Optimizer())
    if method == "polar"
        jump_opf_model(m, ref)
    else
        jump_opf_model_rectangular(m, ref)
    end
    set_ipopt_parameters(m)
    set_optimizer_attribute(m, "output_file", logpath)

    t0 = time()
    optimize!(m)
    t1 = time()

    return t1 - t0, logpath
end

function jump_symbolicad_main(logdir, filename, method)
    logpath = joinpath(logdir, "$(filename)_jump_symbolicad.log")

    mpath = joinpath(@__DIR__, "case", "$(filename).m")
    ref = parse_case(mpath)

    m = JuMP.direct_model(Ipopt.Optimizer())
    if method == "polar"
        jump_opf_model(m, ref)
    else
        jump_opf_model_rectangular(m, ref)
    end
    set_ipopt_parameters(m)
    set_optimizer_attribute(
        m,
        JuMP.MOI.AutomaticDifferentiationBackend(),
        MathOptSymbolicAD.DefaultBackend(),
    )
    set_optimizer_attribute(m, "output_file", logpath)

    t0 = time()
    optimize!(m)
    t1 = time()

    return t1 - t0, logpath
end

function ampl_main(logdir, filename, method)
    logpath = joinpath(logdir, "$(filename)_ampl.log")

    mpath = joinpath(@__DIR__, "case", "$(filename).m")
    ref = parse_case(mpath)

    m = JuMP.Model()
    set_optimizer(m, () -> AmplNLWriter.Optimizer(Ipopt.Ipopt_jll.amplexe))
    if method == "polar"
        jump_opf_model(m, ref)
    else
        jump_opf_model_rectangular(m, ref)
    end
    set_ipopt_parameters(m)
    set_optimizer_attribute(m, "output_file", logpath)

    t0 = time()
    optimize!(m)
    t1 = time()

    return t1 - t0, logpath
end
