import JSON3

function gravity_main(filename, flag, logdir)
    logpath = joinpath(logdir, "$(filename)_gravity.log")

    mpath = joinpath(@__DIR__, "case", "$(filename).m")

    t0 = time()
    file = joinpath(@__DIR__, "acopf.exe")
    redirect_stdio(; stdout = logpath) do
        run(`$file $mpath $flag`)
    end
    t1 = time()

    return t1 - t0, logpath
end

function test_driver()
    json_path = joinpath(@__DIR__, "test_cases.json")
    test_cases = open(json_path, "r") do f
        JSON3.read(json_path)
    end

    # get system args
    args = ARGS

    logdir = joinpath(@__DIR__, "log")
    flag = ""

    if args[end] == "rect"
        logdir = joinpath(@__DIR__, "log_rect")
        flag = "ACRECT"
    end

    for case in test_cases
        println("Running test case: $case")
        total_time, logpath = gravity_main(case, flag, logdir)

        open(logpath, "a") do f
            write(f, "Total time: $total_time\n")
        end
    end
end

test_driver()
