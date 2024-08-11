from itertools import product
import json
from pathlib import Path

import pandas as pd


def extract_timing(f):
    total_time = None
    ipopt_time = None
    ad_time = None
    with open(f, "r") as f:
        lines = f.readlines()
    for line in reversed(lines):
        if total_time and ipopt_time and ad_time:
            break
        if "OverallAlgorithm" in line:
            ipopt_time = float(line.split()[1])
        elif "Function Evaluations" in line:
            ad_time = float(line.split()[2])
        elif "Total time" in line:
            total_time = float(line.split()[2])
    return total_time, ipopt_time, ad_time


def main():
    json_path = Path(__file__).parent / "test_cases.json"

    with open(json_path, "r") as f:
        test_cases = json.load(f)

    # we do not care the first case5 result
    test_cases = test_cases[1:]

    solvers = [
        "ampl",
        "jump",
        "casadi",
        "gravity",
        "jump_symbolicad",
        "examodels",
        "poi",
    ]
    # solvers = ["examodels", "poi"]

    # construct a dataframe, case as row, solver as column
    df = pd.DataFrame(columns=solvers, index=test_cases)

    solvers_details = []
    for solver in solvers:
        for i in range(3):
            solvers_details.append(f"{solver}.{i}")
    wide_df = pd.DataFrame(columns=solvers_details, index=test_cases)

    logdir = Path(__file__).parent / "log"

    for case, solver in product(test_cases, solvers):
        suffix = solver
        if solver == "poi":
            suffix = "poi_LLVM"
        # print(f"Get {case} with {solver}")
        logpath = logdir / f"{case}_{suffix}.log"
        with open(logpath, "r") as f:
            total_time, ipopt_time, ad_time = extract_timing(logpath)

        if total_time:
            content = f"{ad_time:.1f}/{ipopt_time:.1f}/{total_time:.1f}"
        else:
            content = "NA"
        df.at[case, solver] = content

        if total_time:
            wide_df.at[case, f"{solver}.0"] = f"{ad_time:.1f}"
            wide_df.at[case, f"{solver}.1"] = f"{ipopt_time:.1f}"
            wide_df.at[case, f"{solver}.2"] = f"{total_time:.1f}"
        else:
            wide_df.at[case, f"{solver}.0"] = "NA"
            wide_df.at[case, f"{solver}.1"] = "NA"
            wide_df.at[case, f"{solver}.2"] = "NA"

    print(df)
    df.to_csv("result.csv")

    wide_df.to_csv("result_wide.csv")


if __name__ == "__main__":
    main()
