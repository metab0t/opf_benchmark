Install PyOptInterface from the `nlp` branch

Firstly, generate json format of case files
```
mkdir json
mkdir log
julia export_case.jl
```

Run Python benchmark
```
python test_driver.py --solver poi
python test_driver.py --solver casadi
```

Run Julia benchmark
```
julia --project=. test_driver.jl ampl
julia --project=. test_driver.jl jump
```

Show Result
```
python analyze_result.py
```