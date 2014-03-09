# ReverseDiffSparse

[![Build Status](https://travis-ci.org/mlubin/ReverseDiffSparse.jl.png)](https://travis-ci.org/mlubin/ReverseDiffSparse.jl)

Reverse-mode automatic differentiation for closed-form scalar algebraic expressions, producing gradients and Hessians. Efficient sparse Hessian computation is implemented by using state-of-the-art graph coloring approaches, implemented in pure Julia (see ``src/coloring.jl``). This package is primarily used by [JuMP](https://github.com/JuliaOpt/JuMP.jl), but it works stand-alone as well.

Expressions are built around ``Placeholder`` objects (instead of symbols). These placeholders carry variable indices. The ``@processNLExpr`` macro is used to build expression trees. Examples:

```
using ReverseDiffSparse

x = placeholders(3)

ex = @processNLExpr sin(x[1])+sin(x[2])+sin(x[3])
# or, using the specialized sum{} syntax:
# ex = @processNLExpr sum{ sin(x[i]), i = 1:3 }
fg = genfgrad_simple(ex) # generates gradient evaluation function

out = zeros(3)
fval = fg([1.0,2.0,3.0], out) # evaluates the gradient at x = [1,2,3]
# fval = sin(1.0)+sin(2.0)+sin(3.0)
# out now stores [cos(1.0),cos(2.0),cos(3.0)]
```

For an example of Hessian computation, see ``test/test_hessian.jl``. More documentation is forthcoming. This package is experimental and still under development. 

