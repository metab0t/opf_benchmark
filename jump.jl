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

function jump_opf_model(model, ref)
    JuMP.@variable(model, va[i in keys(ref[:bus])])
    JuMP.@variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start = 1.0)

    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    JuMP.@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= p[(l, i, j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= q[(l, i, j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    JuMP.@objective(model, Min, sum(gen["cost"][1] * pg[i]^2 + gen["cost"][2] * pg[i] + gen["cost"][3] for (i, gen) in ref[:gen]))

    for (i, bus) in ref[:ref_buses]
        JuMP.@constraint(model, va[i] == 0)
    end

    for (i, bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        JuMP.@constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts) * vm[i]^2
        )

        JuMP.@constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads) +
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
        JuMP.@constraint(model, p_fr == (g + g_fr) / ttm * vm_fr^2 + (-g * tr + b * ti) / ttm * (vm_fr * vm_to * cos(va_fr - va_to)) + (-b * tr - g * ti) / ttm * (vm_fr * vm_to * sin(va_fr - va_to)))
        JuMP.@constraint(model, q_fr == -(b + b_fr) / ttm * vm_fr^2 - (-b * tr - g * ti) / ttm * (vm_fr * vm_to * cos(va_fr - va_to)) + (-g * tr + b * ti) / ttm * (vm_fr * vm_to * sin(va_fr - va_to)))

        # To side of the branch flow
        JuMP.@constraint(model, p_to == (g + g_to) * vm_to^2 + (-g * tr - b * ti) / ttm * (vm_to * vm_fr * cos(va_to - va_fr)) + (-b * tr + g * ti) / ttm * (vm_to * vm_fr * sin(va_to - va_fr)))
        JuMP.@constraint(model, q_to == -(b + b_to) * vm_to^2 - (-b * tr + g * ti) / ttm * (vm_to * vm_fr * cos(va_to - va_fr)) + (-g * tr - b * ti) / ttm * (vm_to * vm_fr * sin(va_to - va_fr)))

        # Voltage angle difference limit
        JuMP.@constraint(model, branch["angmin"] <= va_fr - va_to <= branch["angmax"])

        # Apparent power limit, from side and to side
        JuMP.@constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        JuMP.@constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
    end
end

function jump_main(filename)
    logpath = joinpath(@__DIR__, "log", "$(filename)_jump.log")

    mpath = joinpath(@__DIR__, "case", "$(filename).m")
    ref = parse_case(mpath)

    m = JuMP.direct_model(Ipopt.Optimizer())
    jump_opf_model(m, ref)
    set_optimizer_attribute(m, "print_level", 5)
    set_optimizer_attribute(m, "linear_solver", "ma27")
    set_optimizer_attribute(m, "hsllib", "D:\\Ipopt\\hsl\\libhsl.dll")
    set_optimizer_attribute(m, "print_timing_statistics", "yes")
    set_optimizer_attribute(m, "output_file", logpath)

    t0 = time()
    optimize!(m)
    t1 = time()

    return t1 - t0, logpath
end

function ampl_main(filename)
    logpath = joinpath(@__DIR__, "log", "$(filename)_ampl.log")

    mpath = joinpath(@__DIR__, "case", "$(filename).m")
    ref = parse_case(mpath)

    m = JuMP.Model()
    set_optimizer(
        m,
        () -> AmplNLWriter.Optimizer(Ipopt.Ipopt_jll.amplexe)
    )
    jump_opf_model(m, ref)
    set_optimizer_attribute(m, "print_level", 5)
    set_optimizer_attribute(m, "linear_solver", "ma27")
    set_optimizer_attribute(m, "hsllib", "libhsl.dll")
    set_optimizer_attribute(m, "print_timing_statistics", "yes")
    set_optimizer_attribute(m, "output_file", logpath)

    t0 = time()
    optimize!(m)
    t1 = time()

    return t1 - t0, logpath
end
