# ReverseDiffSparse

[![Build Status](https://travis-ci.org/mlubin/ReverseDiffSparse.jl.svg?branch=master)](https://travis-ci.org/mlubin/ReverseDiffSparse.jl) [![codecov.io](https://codecov.io/github/mlubin/ReverseDiffSparse.jl/coverage.svg?branch=master)](https://codecov.io/github/mlubin/ReverseDiffSparse.jl?branch=master)

Reverse-mode automatic differentiation for closed-form scalar algebraic expressions, producing gradients and Hessians. Efficient sparse Hessian computation is implemented by using state-of-the-art graph coloring approaches, written in pure Julia (see ``src/coloring.jl``). This package is primarily used by [JuMP](https://github.com/JuliaOpt/JuMP.jl), and we recommend using JuMP as the easiest interface to access it.
