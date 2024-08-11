from itertools import product
import json
from pathlib import Path

import pandas as pd


def extract_timing(f):
    total_time = None
    ipopt_time = None
    ad_time = None
    n_iter = None
    with open(f, "r") as f:
        lines = f.readlines()
    for line in reversed(lines):
        if total_time and ipopt_time and ad_time and n_iter:
            break
        if "OverallAlgorithm" in line:
            ipopt_time = float(line.split()[1])
        elif "Function Evaluations" in line:
            ad_time = float(line.split()[2])
        elif "Total time" in line:
            total_time = float(line.split()[2])
        elif "Number of Iterations....:" in line:
            n_iter = int(line.split()[3])
    return total_time, ipopt_time, ad_time, n_iter


def main():
    json_path = Path(__file__).parent / "test_cases.json"

    with open(json_path, "r") as f:
        test_cases = json.load(f)

    # we do not care the first case5 result
    test_cases = test_cases[1:]
    # case30000 does not converge under cartesian coordinates
    test_cases = [s for s in test_cases if "case30000" not in s]

    solvers = [
        "ampl",
        "jump",
        "casadi",
        "gravity",
        "jump_symbolicad",
        "examodels",
        "poi",
    ]

    solvers_details = []
    for solver in solvers:
        for i in range(3):
            solvers_details.append(f"{solver}.{i}")
    wide_df_polar = pd.DataFrame(columns=solvers_details, index=test_cases)
    wide_df_cartesian = pd.DataFrame(columns=solvers_details, index=test_cases)

    logdir_polar = Path(__file__).parent / "log"
    logdir_cartesian = Path(__file__).parent / "log_rect"

    def f(logdir, wide_df):
        for case, solver in product(test_cases, solvers):
            suffix = solver
            if solver == "poi":
                suffix = "poi_LLVM"
            
            logpath = logdir / f"{case}_{suffix}.log"
            with open(logpath, "r") as f:
                total_time, ipopt_time, ad_time, n_iter = extract_timing(logpath)

            if total_time:
                wide_df.at[case, f"{solver}.0"] = f"{n_iter}"
                wide_df.at[case, f"{solver}.1"] = f"{ad_time:.1f}"
                wide_df.at[case, f"{solver}.2"] = f"{ipopt_time:.1f}"
            else:
                wide_df.at[case, f"{solver}.0"] = "NA"
                wide_df.at[case, f"{solver}.1"] = "NA"
                wide_df.at[case, f"{solver}.2"] = "NA"

    f(logdir_polar, wide_df_polar)
    f(logdir_cartesian, wide_df_cartesian)

    wide_df_polar.to_csv("wide_df_polar.csv")
    wide_df_cartesian.to_csv("wide_df_cartesian.csv")


if __name__ == "__main__":
    main()
