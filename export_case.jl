include("parse_case.jl")
import JSON3

function main(case_file, json_file)
    data = parse_case_data(case_file)
    open(json_file, "w") do io
        JSON3.write(io, data)
    end
end

# iterate all .m files in 'case' folder and write '.json' file in 'json' folder
case_dir = joinpath(@__DIR__, "case")
json_dir = joinpath(@__DIR__, "json")
for file in readdir(case_dir)
    if endswith(file, ".m")
        case_file = joinpath(case_dir, file)
        json_file = joinpath(json_dir, replace(file, ".m" => ".json"))
        main(case_file, json_file)
        @info "exported $json_file"
    end
end
