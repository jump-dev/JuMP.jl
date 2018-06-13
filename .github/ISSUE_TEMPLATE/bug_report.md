---
name: Bug report
about: Help us track down bugs in JuMP

---

Welcome to JuMP! The JuMP developers use the GitHub issue tracker for bug reports and feature requests only. Please read the following before posting a new issue:

- If you have a question or are unsure if the behavior you're experiencing is a bug, please search or post to our Discourse site: https://discourse.julialang.org/c/domain/opt. Questions posted to Discourse have broader visibility and are likely to be answered more quickly than issues filed here.

- If you're experiencing a bug that is solver-specific (e.g., it only happens when you use Gurobi), the issue is best raised at the relevant solver-wrapper repository (e.g., Gurobi.jl). If you are unsure, you should raise the problem on Discourse first.

- If you are reasonably confident your issue is a bug in JuMP, this is the right place. You can find a bug report template below. Be sure to include as much relevant information as possible, including a minimal reproducible example. See https://help.github.com/articles/basic-writing-and-formatting-syntax/ for background on how to format text and code on GitHub issues. Make sure to delete the introductory text (i.e., everything above `## Bug Report:` as well).

Thanks for contributing to JuMP!

## Bug Report: [descriptive name]

I have:
- [ ] searched discourse for related questions;
- [ ] searched Github for existing issues;
- [ ] included a minimum working example; and
- [ ] tried the example on different solvers.

### Description

Give a short description of the bug.

[When I try to add a variable with a lower bound of zero, I get a MethodError saying `no method matching getname(::Int64)`.]

### Minimum Working Example

Replace the following code with your minimum working example. (Note, this example will throw an error and is not a bug. It is expected behavior.)
```julia
using JuMP
m = Model()
@variable(m, 0 <= x)
```

Copy any output (including the full error message) here:
```julia
julia> using JuMP, Gurobi

julia> m = Model()
Feasibility problem with:
 * 0 linear constraints
 * 0 variables
Solver is default solver

julia> @variable(m, 0 <= x)
ERROR: MethodError: no method matching getname(::Int64)
Closest candidates are:
  getname(::Expr) at C:\Users\Oscar\.julia\v0.6\JuMP\src\macros.jl:238
  getname(::Void) at C:\Users\Oscar\.julia\v0.6\JuMP\src\macros.jl:235
  getname(::Symbol) at C:\Users\Oscar\.julia\v0.6\JuMP\src\macros.jl:234
  ...
```

### Version Info

Replace the following with your output from `versioninfo()`.
```
julia> versioninfo()
Julia Version 0.6.3
Commit d55cadc350* (2018-05-28 20:20 UTC)
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
  CPU: Intel(R) Core(TM) i7-3537U CPU @ 2.00GHz
  WORD_SIZE: 64
  BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Sandybridge)
  LAPACK: libopenblas64_
  LIBM: libopenlibm
  LLVM: libLLVM-3.9.1 (ORCJIT, ivybridge)
```
