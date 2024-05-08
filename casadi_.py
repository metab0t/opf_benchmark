# conda create --name casadi python=3.8 -y
# conda activate casadi
# pip install casadi

import casadi
from pathlib import Path
import json
import time


def solve_opf(data, options):
    x, x0, lbx, ubx, cons, lbg, ubg = [], [], [], [], [], [], []

    va, vm = {}, {}
    Nbus = len(data["bus"])
    vmin = data["vmin"]
    vmax = data["vmax"]
    for k in range(Nbus):
        va[k] = casadi.SX.sym(f"va{k}")
        x.append(va[k])
        x0.append(0.0)
        lbx.append(-casadi.inf)
        ubx.append(casadi.inf)
        vm[k] = casadi.SX.sym(f"vm{k}")
        x.append(vm[k])
        x0.append(1.0)
        lbx.append(vmin[k])
        ubx.append(vmax[k])

    pg, qg = {}, {}
    Ngen = len(data["gen"])
    pmin = data["pmin"]
    pmax = data["pmax"]
    qmin = data["qmin"]
    qmax = data["qmax"]
    for k in range(Ngen):
        pg[k] = casadi.SX.sym(f"pg{k}")
        x.append(pg[k])
        x0.append(0.0)
        lbx.append(pmin[k])
        ubx.append(pmax[k])
        qg[k] = casadi.SX.sym(f"qg{k}")
        x.append(qg[k])
        x0.append(0.0)
        lbx.append(qmin[k])
        ubx.append(qmax[k])

    p, q = {}, {}
    Narc = len(data["arc"])
    rate_a = data["rate_a"]
    for k in range(Narc):
        a = rate_a[k]
        p[k] = casadi.SX.sym(f"p{k}")
        x.append(p[k])
        x0.append(0.0)
        lbx.append(-a)
        ubx.append(a)
        q[k] = casadi.SX.sym(f"q{k}")
        x.append(q[k])
        x0.append(0.0)
        lbx.append(-a)
        ubx.append(a)

    f = sum(
        g["cost1"] * pg[g["i"] - 1] * pg[g["i"] - 1]
        + g["cost2"] * pg[g["i"] - 1]
        + g["cost3"]
        for g in data["gen"]
    )

    for k in data["ref_buses"]:
        cons.append(va[k - 1])
        lbg.append(0)
        ubg.append(0)

    p_balance_expr = [0.0 for _ in range(Nbus)]
    q_balance_expr = [0.0 for _ in range(Nbus)]
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
        cons.append(p_balance_expr[i])
        lbg.append(0)
        ubg.append(0)
        cons.append(q_balance_expr[i])
        lbg.append(0)
        ubg.append(0)

    angmin = data["angmin"]
    angmax = data["angmax"]
    for i, br in enumerate(data["branch"]):
        f_idx = br["f_idx"] - 1
        t_idx = br["t_idx"] - 1
        f_bus = br["f_bus"] - 1
        t_bus = br["t_bus"] - 1
        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]
        vm_fr = vm[f_bus]
        vm_to = vm[t_bus]
        va_fr = va[f_bus]
        va_to = va[t_bus]
        g, b = br["g"], br["b"]
        tr, ti = br["tr"], br["ti"]
        ttm = tr**2 + ti**2
        g_fr = br["g_fr"]
        b_fr = br["b_fr"]
        g_to = br["g_to"]
        b_to = br["b_to"]
        cons.append(
            (g + g_fr) / ttm * vm_fr**2
            + (-g * tr + b * ti) / ttm * (vm_fr * vm_to * casadi.cos(va_fr - va_to))
            + (-b * tr - g * ti) / ttm * (vm_fr * vm_to * casadi.sin(va_fr - va_to))
            - p_fr
        )
        cons.append(
            -(b + b_fr) / ttm * vm_fr**2
            - (-b * tr - g * ti) / ttm * (vm_fr * vm_to * casadi.cos(va_fr - va_to))
            + (-g * tr + b * ti) / ttm * (vm_fr * vm_to * casadi.sin(va_fr - va_to))
            - q_fr
        )
        cons.append(
            (g + g_to) * vm_to**2
            + (-g * tr - b * ti) / ttm * (vm_to * vm_fr * casadi.cos(va_to - va_fr))
            + (-b * tr + g * ti) / ttm * (vm_to * vm_fr * casadi.sin(va_to - va_fr))
            - p_to
        )
        cons.append(
            -(b + b_to) * vm_to**2
            - (-b * tr + g * ti) / ttm * (vm_to * vm_fr * casadi.cos(va_to - va_fr))
            + (-g * tr - b * ti) / ttm * (vm_to * vm_fr * casadi.sin(va_to - va_fr))
            - q_to
        )
        for i in range(4):
            lbg.append(0)
            ubg.append(0)

        cons.append(va_fr - va_to)
        lbg.append(angmin[i])
        ubg.append(angmax[i])

        rate_a_sq = br["rate_a_sq"]
        cons.append(p_fr**2 + q_fr**2)
        lbg.append(-casadi.inf)
        ubg.append(rate_a_sq)
        cons.append(p_to**2 + q_to**2)
        lbg.append(-casadi.inf)
        ubg.append(rate_a_sq)

    model = casadi.nlpsol(
        "model",
        "ipopt",
        {"x": casadi.vcat(x), "f": f, "g": casadi.vcat(cons)},
        options,
    )
    solution = model(
        lbx=lbx,
        ubx=ubx,
        lbg=lbg,
        ubg=ubg,
        x0=x0,
    )
    return solution


def casadi_main(filename, jit = False):
    # filename = "pglib_opf_case30000_goc"
    # filename = "pglib_opf_case78484_epigrids"
    # filename = "pglib_opf_case10480_goc"
    # filename = "pglib_opf_case19402_goc"

    json_path = Path(__file__).parent / "json" / f"{filename}.json"
    json_path = json_path.resolve()
    with open(json_path, "r") as f:
        data = json.load(f)

    logpath =  str(Path(__file__).parent / "log" / f"{filename}_casadi.log")
    options = {
        "ipopt.print_level": 5,
        "ipopt.hsllib": "libhsl.dll",
        "ipopt.linear_solver": "ma27",
        "ipopt.print_timing_statistics": "yes",
        "ipopt.output_file": logpath,
    }

    if jit:
        flags = ["/O2"]
        compiler = "cl"
        jit_options = {"flags": flags, "verbose": True, "compiler": compiler}
        options.update({"jit": True, "compiler": "shell", "jit_options": jit_options})

    t0 = time.time()
    solve_opf(data, options)
    t1 = time.time()
    
    return t1 - t0, logpath


if __name__ == "__main__":
    casadi_main("pglib_opf_case6495_rte")
