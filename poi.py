import pyoptinterface as poi
from pyoptinterface import ipopt

sin = poi.sin
cos = poi.cos


from pathlib import Path
import json
import time
import math

libipopt_path = r"D:\Ipopt\bin\ipopt-3.dll"

ipopt.load_library(libipopt_path)


def get_ipopt_model():
    model = ipopt.Model()
    model.set_raw_parameter("print_level", 5)
    model.set_raw_parameter("hsllib", "libhsl.dll")
    model.set_raw_parameter("linear_solver", "ma27")
    model.set_raw_parameter("print_timing_statistics", "yes")
    model.set_raw_parameter("max_iter", 200)

    return model


def solve_opf(model, data, jit_engine="LLVM"):
    Nbus = len(data["bus"])
    va = [model.add_variable() for _ in range(Nbus)]

    vmin = data["vmin"]
    vmax = data["vmax"]
    vm = [model.add_variable(lb=vmin[i], ub=vmax[i], start=1.0) for i in range(Nbus)]

    Ngen = len(data["gen"])
    pmin = data["pmin"]
    pmax = data["pmax"]
    pg = [model.add_variable(lb=pmin[i], ub=pmax[i]) for i in range(Ngen)]
    qmin = data["qmin"]
    qmax = data["qmax"]
    qg = [model.add_variable(lb=qmin[i], ub=qmax[i]) for i in range(Ngen)]

    Narc = len(data["arc"])
    rate_a = data["rate_a"]
    p = [model.add_variable(lb=-rate_a[i], ub=rate_a[i]) for i in range(Narc)]
    q = [model.add_variable(lb=-rate_a[i], ub=rate_a[i]) for i in range(Narc)]

    for g in data["gen"]:
        expr = poi.ExprBuilder()
        v = pg[g["i"] - 1]
        expr += g["cost1"] * v * v + g["cost2"] * v + g["cost3"]
        model.add_objective(expr)

    for i in data["ref_buses"]:
        model.add_linear_constraint(1.0 * va[i - 1], poi.Eq, 0.0)

    def branch_flow(x, p):
        p_from, p_to = x["p_from"], x["p_to"]
        q_from, q_to = x["q_from"], x["q_to"]
        vm_from, vm_to = x["vm_from"], x["vm_to"]
        va_from, va_to = x["va_from"], x["va_to"]
        c1, c2, c3, c4, c5, c6, c7, c8 = p

        va_delta = va_from - va_to
        sin_ft = sin(va_delta)
        cos_ft = cos(va_delta)
        sin_tf = -sin_ft
        cos_tf = cos_ft

        pij = (
            p_from
            - c5 * vm_from * vm_from
            - c3 * vm_from * vm_to * cos_ft
            - c4 * vm_from * vm_to * sin_ft
        )
        qij = (
            q_from
            + c6 * vm_from * vm_from
            + c4 * vm_from * vm_to * cos_ft
            - c3 * vm_from * vm_to * sin_ft
        )
        pji = (
            p_to
            - c7 * vm_to * vm_to
            - c1 * vm_from * vm_to * cos_tf
            - c2 * vm_from * vm_to * sin_tf
        )
        qji = (
            q_to
            + c8 * vm_to * vm_to
            + c2 * vm_from * vm_to * cos_tf
            - c1 * vm_from * vm_to * sin_tf
        )

        return [pij, qij, pji, qji]

    branch_flow_f = model.register_function(
        branch_flow,
        var=[
            "p_from",
            "q_from",
            "p_to",
            "q_to",
            "vm_from",
            "vm_to",
            "va_from",
            "va_to",
        ],
        param=8,
        name="flow",
    )

    branch_parameters = []
    for b in data["branch"]:
        d = []
        for i in range(1, 9):
            value = b[f"c{i}"]
            d.append(value)
        branch_parameters.append(d)

    for b, bp in zip(data["branch"], branch_parameters):
        f_idx = b["f_idx"] - 1
        t_idx = b["t_idx"] - 1
        f_bus = b["f_bus"] - 1
        t_bus = b["t_bus"] - 1
        x = [
            p[f_idx],
            q[f_idx],
            p[t_idx],
            q[t_idx],
            vm[f_bus],
            vm[t_bus],
            va[f_bus],
            va[t_bus],
        ]
        param = bp
        con = model.add_nl_constraint(
            branch_flow_f, x, param, poi.Eq, [0.0, 0.0, 0.0, 0.0]
        )

    angmin = data["angmin"]
    angmax = data["angmax"]
    for i, b in enumerate(data["branch"]):
        f_idx = b["f_idx"] - 1
        t_idx = b["t_idx"] - 1
        f_bus = b["f_bus"] - 1
        t_bus = b["t_bus"] - 1
        rate_a_sq = b["rate_a_sq"]
        model.add_linear_constraint(va[f_bus] - va[t_bus], poi.In, angmin[i], angmax[i])
        model.add_quadratic_constraint(
            p[f_idx] * p[f_idx] + q[f_idx] * q[f_idx], poi.Leq, rate_a_sq
        )
        model.add_quadratic_constraint(
            p[t_idx] * p[t_idx] + q[t_idx] * q[t_idx], poi.Leq, rate_a_sq
        )

    # bus balance constraint
    p_balance_expr = [poi.ExprBuilder() for _ in range(Nbus)]
    q_balance_expr = [poi.ExprBuilder() for _ in range(Nbus)]
    for i, b in enumerate(data["bus"]):
        p_balance_expr[i] += b["pd"] + b["gs"] * vm[i] * vm[i]
        q_balance_expr[i] += b["qd"] - b["bs"] * vm[i] * vm[i]
    for a in data["arc"]:
        bus = a["bus"] - 1
        i = a["i"] - 1
        p_balance_expr[bus] += p[i]
        q_balance_expr[bus] += q[i]
    for g in data["gen"]:
        bus = g["bus"] - 1
        i = g["i"] - 1
        p_balance_expr[bus] -= pg[i]
        q_balance_expr[bus] -= qg[i]

    for i in range(Nbus):
        model.add_quadratic_constraint(p_balance_expr[i], poi.Eq, 0.0)
        model.add_quadratic_constraint(q_balance_expr[i], poi.Eq, 0.0)

    model.optimize(jit_engine=jit_engine)


def solve_opf_rectangular(model, data, jit_engine="LLVM"):
    Nbus = len(data["bus"])

    vmin = data["vmin"]
    vmax = data["vmax"]
    vr = [model.add_variable(lb=-vmax[i], ub=vmax[i], start=1.0) for i in range(Nbus)]
    vi = [model.add_variable(lb=-vmax[i], ub=vmax[i], start=0.0) for i in range(Nbus)]

    Ngen = len(data["gen"])
    pmin = data["pmin"]
    pmax = data["pmax"]
    pg = [model.add_variable(lb=pmin[i], ub=pmax[i]) for i in range(Ngen)]
    qmin = data["qmin"]
    qmax = data["qmax"]
    qg = [model.add_variable(lb=qmin[i], ub=qmax[i]) for i in range(Ngen)]

    Narc = len(data["arc"])
    rate_a = data["rate_a"]
    p = [model.add_variable(lb=-rate_a[i], ub=rate_a[i]) for i in range(Narc)]
    q = [model.add_variable(lb=-rate_a[i], ub=rate_a[i]) for i in range(Narc)]

    for g in data["gen"]:
        expr = poi.ExprBuilder()
        v = pg[g["i"] - 1]
        expr += g["cost1"] * v * v + g["cost2"] * v + g["cost3"]
        model.add_objective(expr)

    for i in data["ref_buses"]:
        model.add_linear_constraint(1.0 * vi[i - 1], poi.Eq, 0.0)

    for i, b in enumerate(data["bus"]):
        model.add_quadratic_constraint(
            vr[i] * vr[i] + vi[i] * vi[i], poi.In, vmin[i] * vmin[i], vmax[i] * vmax[i]
        )

    def branch_flow(x, p):
        p_fr, p_to = x["p_from"], x["p_to"]
        q_fr, q_to = x["q_from"], x["q_to"]
        vr_fr, vr_to = x["vr_from"], x["vr_to"]
        vi_fr, vi_to = x["vi_from"], x["vi_to"]
        c1, c2, c3, c4, c5, c6, c7, c8 = p

        vfr_sq = vr_fr * vr_fr + vi_fr * vi_fr
        vto_sq = vr_to * vr_to + vi_to * vi_to
        re_vfr_vto = vr_fr * vr_to + vi_fr * vi_to
        im_vfr_vto = vi_fr * vr_to - vr_fr * vi_to
        re_vto_vfr = re_vfr_vto
        im_vto_vfr = -im_vfr_vto

        pij = p_fr - c5 * vfr_sq - c3 * re_vfr_vto - c4 * im_vfr_vto
        qij = q_fr + c6 * vfr_sq + c4 * re_vfr_vto - c3 * im_vfr_vto
        pji = p_to - c7 * vto_sq - c1 * re_vto_vfr - c2 * im_vto_vfr
        qji = q_to + c8 * vto_sq + c2 * re_vto_vfr - c1 * im_vto_vfr

        return [pij, qij, pji, qji]

    branch_flow_f = model.register_function(
        branch_flow,
        var=[
            "p_from",
            "q_from",
            "p_to",
            "q_to",
            "vr_from",
            "vr_to",
            "vi_from",
            "vi_to",
        ],
        param=8,
        name="flow",
    )

    branch_parameters = []
    for k, b in enumerate(data["branch"]):
        d = []
        for i in range(1, 9):
            value = b[f"c{i}"]
            d.append(value)
        branch_parameters.append(d)

    for b, bp in zip(data["branch"], branch_parameters):
        f_idx = b["f_idx"] - 1
        t_idx = b["t_idx"] - 1
        f_bus = b["f_bus"] - 1
        t_bus = b["t_bus"] - 1
        x = [
            p[f_idx],
            q[f_idx],
            p[t_idx],
            q[t_idx],
            vr[f_bus],
            vr[t_bus],
            vi[f_bus],
            vi[t_bus],
        ]
        param = bp
        con = model.add_nl_constraint(
            branch_flow_f, x, param, poi.Eq, [0.0, 0.0, 0.0, 0.0]
        )

    def branch_angle(x):
        vr_fr, vr_to = x["vr_from"], x["vr_to"]
        vi_fr, vi_to = x["vi_from"], x["vi_to"]

        re_vfr_vto = vr_fr * vr_to + vi_fr * vi_to
        im_vfr_vto = vi_fr * vr_to - vr_fr * vi_to

        angdiff = im_vfr_vto / re_vfr_vto

        return [angdiff]

    branch_angle_f = model.register_function(
        branch_angle,
        var=[
            "vr_from",
            "vr_to",
            "vi_from",
            "vi_to",
        ],
        name="angle",
    )
    angmin = data["angmin"]
    angmax = data["angmax"]
    for i, b in enumerate(data["branch"]):
        f_idx = b["f_idx"] - 1
        t_idx = b["t_idx"] - 1
        f_bus = b["f_bus"] - 1
        t_bus = b["t_bus"] - 1
        x = [
            vr[f_bus],
            vr[t_bus],
            vi[f_bus],
            vi[t_bus],
        ]
        con = model.add_nl_constraint(
            branch_angle_f, x, poi.In, [math.tan(angmin[i])], [math.tan(angmax[i])]
        )

    for i, b in enumerate(data["branch"]):
        f_idx = b["f_idx"] - 1
        t_idx = b["t_idx"] - 1
        f_bus = b["f_bus"] - 1
        t_bus = b["t_bus"] - 1
        rate_a_sq = b["rate_a_sq"]
        model.add_quadratic_constraint(
            p[f_idx] * p[f_idx] + q[f_idx] * q[f_idx], poi.Leq, rate_a_sq
        )
        model.add_quadratic_constraint(
            p[t_idx] * p[t_idx] + q[t_idx] * q[t_idx], poi.Leq, rate_a_sq
        )

    # bus balance constraint
    p_balance_expr = [poi.ExprBuilder() for _ in range(Nbus)]
    q_balance_expr = [poi.ExprBuilder() for _ in range(Nbus)]
    for i, b in enumerate(data["bus"]):
        p_balance_expr[i] += b["pd"] + b["gs"] * (vr[i] * vr[i] + vi[i] * vi[i])
        q_balance_expr[i] += b["qd"] - b["bs"] * (vr[i] * vr[i] + vi[i] * vi[i])
    for a in data["arc"]:
        bus = a["bus"] - 1
        i = a["i"] - 1
        p_balance_expr[bus] += p[i]
        q_balance_expr[bus] += q[i]
    for g in data["gen"]:
        bus = g["bus"] - 1
        i = g["i"] - 1
        p_balance_expr[bus] -= pg[i]
        q_balance_expr[bus] -= qg[i]

    for i in range(Nbus):
        model.add_quadratic_constraint(p_balance_expr[i], poi.Eq, 0.0)
        model.add_quadratic_constraint(q_balance_expr[i], poi.Eq, 0.0)

    model.optimize(jit_engine=jit_engine)


def poi_main(logdir, filename, method, jit_engine="LLVM"):
    # filename = "pglib_opf_case30000_goc"
    # filename = "pglib_opf_case78484_epigrids"
    # filename = "pglib_opf_case10480_goc"
    # filename = "pglib_opf_case19402_goc"
    model = get_ipopt_model()
    logpath = str(logdir / f"{filename}_poi_{jit_engine}.log")
    model.set_raw_parameter("output_file", logpath)

    json_path = Path(__file__).parent / "json" / f"{filename}.json"
    json_path = json_path.resolve()
    with open(json_path, "r") as f:
        data = json.load(f)

    t0 = time.time()
    if method == "polar":
        solve_opf(model, data, jit_engine=jit_engine)
    else:
        solve_opf_rectangular(model, data, jit_engine=jit_engine)
    t1 = time.time()

    # with open("jit.c", "w") as f:
    #     f.write(model.jit_compiler.source_code)

    return t1 - t0, logpath


if __name__ == "__main__":
    poi_main("pglib_opf_case5_pjm", "C")
