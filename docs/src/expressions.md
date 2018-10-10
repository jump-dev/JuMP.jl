```@meta
DocTestSetup = quote
    using JuMP
end
```

Expressions
===========

DRAFT: JuMP has multiple types of expressions: affine, quadratic, and nonlinear.
Just talk about affine and quadratic here (see other section for nonlinear).
Describe the basic data structures and ways to construct expressions
(these are: constructors, operators, `add_to_expression!`, and macros).


Example code with doc tests:

```jldoctest
m = Model()
@variable(m, x)
@variable(m, y)
AffExpr(-1.0, x => 2.0, y => 1.0)

# output

2 x + y - 1
```

Objective functions
-------------------

TODO: Describe how JuMP expressions relate to MOI functions. How to set, query,
and modify an objective function.
