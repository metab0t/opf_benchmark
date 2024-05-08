from itertools import product
import json
from pathlib import Path

import pandas as pd

def extract_timing(f):
    total_time = None
    ad_time = None
    with open(f, "r") as f:
        lines = f.readlines()
    for line in reversed(lines):
        if total_time is not None and ad_time is not None:
            break
        if "OverallAlgorithm" in line:
            total_time = float(line.split()[1])
        elif "Function Evaluations" in line:
            ad_time = float(line.split()[2])
    return total_time, ad_time

def main():
    json_path = Path(__file__).parent / "test_cases.json"

    with open(json_path, "r") as f:
        test_cases = json.load(f)

    # we do not care the first case5 result
    test_cases = test_cases[1:]

    solvers = ["ampl", "jump", "casadi", "poi"]

    # construct a dataframe, case as row, solver as column
    df = pd.DataFrame(columns=solvers, index=test_cases)

    for (case, solver) in product(test_cases, solvers):
        suffix = solver
        if solver == "poi":
            suffix = "poi_LLVM"
        # print(f"Get {case} with {solver}")
        logpath = Path(__file__).parent / "log" / f"{case}_{suffix}.log"
        with open(logpath, "r") as f:
            total_time, ad_time = extract_timing(logpath)
        
        content = f"{total_time:.1f}({ad_time:.1f})"
        df.at[case, solver] = content

    print(df)
    df.to_csv("result.csv")

if __name__ == "__main__":
    main()