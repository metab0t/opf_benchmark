from itertools import product
import json
from pathlib import Path

import pandas as pd


def extract_timing(f):
    total_time = None
    ipopt_time = None
    ad_time = None
    number_of_iterations = None
    objective_value = None
    with open(f, "r") as f:
        lines = f.readlines()
    for line in reversed(lines):
        if total_time and ipopt_time and ad_time and number_of_iterations and objective_value:
            break
        if "OverallAlgorithm" in line:
            ipopt_time = float(line.split()[1])
        elif "Function Evaluations" in line:
            ad_time = float(line.split()[2])
        elif "Total time" in line:
            total_time = float(line.split()[2])
        elif "Number of Iterations....:" in line:
            number_of_iterations = int(line.split()[3])
        elif "Objective...............:" in line:
            objective_value = float(line.split()[2])
    return total_time, ipopt_time, ad_time, number_of_iterations, objective_value


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

    # construct a dataframe, case as row, solver as column
    df = pd.DataFrame(columns=solvers, index=test_cases)

    solvers_details = []
    for solver in solvers:
        for i in range(2):
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
            total_time, ipopt_time, ad_time, number_of_iterations, objective_value = extract_timing(logpath)

        if total_time:
            content = f"{ad_time:.1f}/{ipopt_time:.1f}/{total_time:.1f}"
            noad_time = ipopt_time - ad_time
            content = f"{noad_time:.1f}"
        else:
            content = "NA"
        df.at[case, solver] = content

        if total_time:
            wide_df.at[case, f"{solver}.0"] = f"{number_of_iterations}"
            wide_df.at[case, f"{solver}.1"] = f"{objective_value:.3e}"
        else:
            wide_df.at[case, f"{solver}.0"] = "NA"
            wide_df.at[case, f"{solver}.1"] = "NA"

    print(df)
    df.to_csv("result_noadnlp.csv")

    wide_df.to_csv("result_iter_obj.csv")


if __name__ == "__main__":
    main()
