from casadi_ import casadi_main
from poi import poi_main

import json
from pathlib import Path
import argparse

def main():
    json_path = Path(__file__).parent / "test_cases.json"

    with open(json_path, "r") as f:
        test_cases = json.load(f)

    # select solver method
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", type=str, default="poi")
    parser.add_argument("--method", type=str, default="polar")
    args = parser.parse_args()

    solver = args.solver
    method = args.method
    logdir = Path(__file__).parent / "log"

    if solver == "casadi":
        f = casadi_main
    elif solver == "poi":
        f = poi_main
    else:
        raise ValueError(f"Unknown solver: {solver}")
    
    if method == "rect":
        logdir = Path(__file__).parent / "log_rect"

    for case in test_cases:
        print(f"Running {case}")
        total_time, logpath = f(logdir, case, method)
        with open(logpath, "a") as file:
            file.write(f"Total time: {total_time}\n")

if __name__ == "__main__":
    main()
