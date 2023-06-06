# NLopt.jl

[![Build Status](https://github.com/JuliaOpt/NLopt.jl/workflows/CI/badge.svg?branch=master)](https://github.com/JuliaOpt/NLopt.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaOpt/NLopt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaOpt/NLopt.jl)

[NLopt.jl](https://github.com/JuliaOpt/NLopt.jl) is a wrapper for the
[NLopt](https://nlopt.readthedocs.io/en/latest/) library.

## License

`NLopt.jl` is licensed under the [MIT License](https://github.com/JuliaOpt/NLopt.jl/blob/master/LICENSE.md).

The underlying solver, [stevengj/nlopt](https://github.com/stevengj/nlopt), is
licensed under the [LGPL v3.0 license](https://github.com/stevengj/nlopt/blob/master/COPYING).

## Installation

Install `NLopt.jl` using the Julia package manager:
```julia
import Pkg
Pkg.add("NLopt")
```

In addition to installing the `NLopt.jl` package, this will also download and
install the NLopt binaries. You do not need to install NLopt separately.

## Use with JuMP

You can use NLopt with JuMP as follows:
```julia
using JuMP, NLopt
model = Model(NLopt.Optimizer)
set_attribute(model, "algorithm", :LD_MMA)
```

## Options

The `algorithm` attribute is required. The value must be one of the supported
[NLopt algorithms](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/).

## Documentation

For more details, see the [NLopt.jl README](https://github.com/JuliaOpt/NLopt.jl/blob/master/README.md)
or the [NLopt documentation](https://nlopt.readthedocs.io/en/latest/).
