#!/usr/bin/env julia
###### AC-OPF using ExaModels ######
#
# implementation reference: https://exanauts.github.io/ExaModels.jl/stable/guide/
# only the built-in AD library is supported
#

include("parse_case.jl")
import ExaModels
import NLPModelsIpopt
import LinearAlgebra

function solve_opf(data, ipopt_options)
    w = ExaModels.ExaCore()

    va = ExaModels.variable(w, length(data.bus);)

    vm = ExaModels.variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus, Float64), 1.0),
        lvar = data.vmin,
        uvar = data.vmax,
    )
    pg = ExaModels.variable(w, length(data.gen); lvar = data.pmin, uvar = data.pmax)

    qg = ExaModels.variable(w, length(data.gen); lvar = data.qmin, uvar = data.qmax)

    p = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    q = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    o = ExaModels.objective(w, g.cost1 * pg[g.i]^2 + g.cost2 * pg[g.i] + g.cost3 for g in data.gen)

    c1 = ExaModels.constraint(w, va[i] for i in data.ref_buses)

    c2 = ExaModels.constraint(
        w,
        p[b.f_idx] - b.c5 * vm[b.f_bus]^2 -
        b.c3 * (vm[b.f_bus] * vm[b.t_bus] * cos(va[b.f_bus] - va[b.t_bus])) -
        b.c4 * (vm[b.f_bus] * vm[b.t_bus] * sin(va[b.f_bus] - va[b.t_bus])) for b in data.branch
    )

    c3 = ExaModels.constraint(
        w,
        q[b.f_idx] +
        b.c6 * vm[b.f_bus]^2 +
        b.c4 * (vm[b.f_bus] * vm[b.t_bus] * cos(va[b.f_bus] - va[b.t_bus])) -
        b.c3 * (vm[b.f_bus] * vm[b.t_bus] * sin(va[b.f_bus] - va[b.t_bus])) for b in data.branch
    )

    c4 = ExaModels.constraint(
        w,
        p[b.t_idx] - b.c7 * vm[b.t_bus]^2 -
        b.c1 * (vm[b.t_bus] * vm[b.f_bus] * cos(va[b.t_bus] - va[b.f_bus])) -
        b.c2 * (vm[b.t_bus] * vm[b.f_bus] * sin(va[b.t_bus] - va[b.f_bus])) for b in data.branch
    )

    c5 = ExaModels.constraint(
        w,
        q[b.t_idx] +
        b.c8 * vm[b.t_bus]^2 +
        b.c2 * (vm[b.t_bus] * vm[b.f_bus] * cos(va[b.t_bus] - va[b.f_bus])) -
        b.c1 * (vm[b.t_bus] * vm[b.f_bus] * sin(va[b.t_bus] - va[b.f_bus])) for b in data.branch
    )

    c6 = ExaModels.constraint(
        w,
        va[b.f_bus] - va[b.t_bus] for b in data.branch;
        lcon = data.angmin,
        ucon = data.angmax,
    )
    c7 = ExaModels.constraint(
        w,
        p[b.f_idx]^2 + q[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )
    c8 = ExaModels.constraint(
        w,
        p[b.t_idx]^2 + q[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c9 = ExaModels.constraint(w, b.pd + b.gs * vm[b.i]^2 for b in data.bus)

    c10 = ExaModels.constraint(w, b.qd - b.bs * vm[b.i]^2 for b in data.bus)

    c11 = ExaModels.constraint!(w, c9, a.bus => p[a.i] for a in data.arc)
    c12 = ExaModels.constraint!(w, c10, a.bus => q[a.i] for a in data.arc)

    c13 = ExaModels.constraint!(w, c9, g.bus => -pg[g.i] for g in data.gen)
    c14 = ExaModels.constraint!(w, c10, g.bus => -qg[g.i] for g in data.gen)

    model = ExaModels.ExaModel(w)

    result = NLPModelsIpopt.ipopt(model; ipopt_options...)
end

function solve_opf_rectangular(data, ipopt_options)
    w = ExaModels.ExaCore()

    vr = ExaModels.variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus, Float64), 1.0),
        lvar = -data.vmax,
        uvar = data.vmax,
    )
    vi = ExaModels.variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus, Float64), 0.0),
        lvar = -data.vmax,
        uvar = data.vmax,
    )
    pg = ExaModels.variable(w, length(data.gen); lvar = data.pmin, uvar = data.pmax)

    qg = ExaModels.variable(w, length(data.gen); lvar = data.qmin, uvar = data.qmax)

    p = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    q = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    o = ExaModels.objective(w, g.cost1 * pg[g.i]^2 + g.cost2 * pg[g.i] + g.cost3 for g in data.gen)

    c1 = ExaModels.constraint(w, vi[i] for i in data.ref_buses)
    c_vmag = ExaModels.constraint(
        w,
        vr[i]^2 + vi[i]^2 for i in 1:length(data.bus);
        lcon = data.vmin_sq,
        ucon = data.vmax_sq,
    )

    c2 = ExaModels.constraint(
        w,
        p[b.f_idx] - b.c5 * (vr[b.f_bus]^2 + vi[b.f_bus]^2) -
        b.c3 * (vr[b.f_bus] * vr[b.t_bus] + vi[b.f_bus] * vi[b.t_bus]) -
        b.c4 * (vi[b.f_bus] * vr[b.t_bus] - vr[b.f_bus] * vi[b.t_bus]) for b in data.branch
    )

    c3 = ExaModels.constraint(
        w,
        q[b.f_idx] +
        b.c6 * (vr[b.f_bus]^2 + vi[b.f_bus]^2) +
        b.c4 * (vr[b.f_bus] * vr[b.t_bus] + vi[b.f_bus] * vi[b.t_bus]) -
        b.c3 * (vi[b.f_bus] * vr[b.t_bus] - vr[b.f_bus] * vi[b.t_bus]) for b in data.branch
    )

    c4 = ExaModels.constraint(
        w,
        p[b.t_idx] - b.c7 * (vr[b.t_bus]^2 + vi[b.t_bus]^2) -
        b.c1 * (vr[b.t_bus] * vr[b.f_bus] + vi[b.t_bus] * vi[b.f_bus]) -
        b.c2 * (vi[b.t_bus] * vr[b.f_bus] - vr[b.t_bus] * vi[b.f_bus]) for b in data.branch
    )

    c5 = ExaModels.constraint(
        w,
        q[b.t_idx] +
        b.c8 * (vr[b.t_bus]^2 + vi[b.t_bus]^2) +
        b.c2 * (vr[b.t_bus] * vr[b.f_bus] + vi[b.t_bus] * vi[b.f_bus]) -
        b.c1 * (vi[b.t_bus] * vr[b.f_bus] - vr[b.t_bus] * vi[b.f_bus]) for b in data.branch
    )

    c_angle = ExaModels.constraint(
        w,
        (vi[b.f_bus] * vr[b.t_bus] - vr[b.f_bus] * vi[b.t_bus]) /
        (vr[b.f_bus] * vr[b.t_bus] + vi[b.f_bus] * vi[b.t_bus]) for b in data.branch;
        lcon = tan.(data.angmin),
        ucon = tan.(data.angmax),
    )
    c7 = ExaModels.constraint(
        w,
        p[b.f_idx]^2 + q[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )
    c8 = ExaModels.constraint(
        w,
        p[b.t_idx]^2 + q[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c9 = ExaModels.constraint(w, b.pd + b.gs * (vr[b.i]^2 + vi[b.i]^2) for b in data.bus)

    c10 = ExaModels.constraint(w, b.qd - b.bs * (vr[b.i]^2 + vi[b.i]^2) for b in data.bus)

    c11 = ExaModels.constraint!(w, c9, a.bus => p[a.i] for a in data.arc)
    c12 = ExaModels.constraint!(w, c10, a.bus => q[a.i] for a in data.arc)

    c13 = ExaModels.constraint!(w, c9, g.bus => -pg[g.i] for g in data.gen)
    c14 = ExaModels.constraint!(w, c10, g.bus => -qg[g.i] for g in data.gen)

    model = ExaModels.ExaModel(w)

    result = NLPModelsIpopt.ipopt(model; ipopt_options...)
end

function exa_main(logdir, filename, method)
    logpath = joinpath(logdir, "$(filename)_examodels.log")
    ipopt_options = (
        print_level = 5,
        linear_solver = "ma27",
        hsllib = "libhsl.dll",
        print_timing_statistics = "yes",
        max_iter = 200,
        # output_file=logpath,
    )

    mpath = joinpath(@__DIR__, "case", "$(filename).m")
    data = parse_case_data(mpath)

    t0 = time()
    redirect_stdio(; stdout = logpath) do
        if method == "polar"
            solve_opf(data, ipopt_options)
        else
            solve_opf_rectangular(data, ipopt_options)
        end
    end
    t1 = time()

    return t1 - t0, logpath
end

# main("pglib_opf_case5_pjm")
# main("pglib_opf_case30000_goc")
# main("pglib_opf_case78484_epigrids")
