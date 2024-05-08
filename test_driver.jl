include("jump.jl")

import JSON3

function test_driver()
    json_path = joinpath(@__DIR__, "test_cases.json")
    test_cases = open(json_path, "r") do f
        JSON3.read(json_path)
    end

    # get system args to select solver
    args = ARGS
    if length(args) == 0
        println("No solver selected. Defaulting to JuMP.")
        solver = "jump"
    else
        solver = args[1]
    end

    if solver == "jump"
        println("Using JuMP solver.")
        f = jump_main
    elseif solver == "ampl"
        println("Using AMPL solver.")
        f = ampl_main
    else
        println("Invalid solver selected. Defaulting to JuMP.")
        exit(1)
    end
    
    for case in test_cases
        println("Running test case: $case")
        total_time, logpath = f(case)

        # open(logpath, "a") do f
        #     write(f, "Total time: $total_time\n")
        # end
    end
end

test_driver()
