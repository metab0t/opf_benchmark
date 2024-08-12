Install PyOptInterface:

```
pip install pyoptinterface[nlp]
pip install casadi
```

Download executable of Gravity from [the link](https://github.com/metab0t/Gravity/releases), extract `acopf.exe` to this directory.

Firstly, generate json format of case files

```
mkdir json
mkdir log
julia export_case.jl
```

Run Python benchmark

```
python test_driver.py --solver poi (--method rect)
python test_driver.py --solver casadi (--method rect)
```

Run Julia benchmark

```
julia --project=. test_driver.jl exa
julia --project=. test_driver.jl ampl
julia --project=. test_driver.jl jump
julia --project=. test_driver.jl jump_symbolicad
```

Run Gravity benchmark
```
julia --project=. test_gravity_driver.jl
```

Show Result

```
python analyze_result.py
```
