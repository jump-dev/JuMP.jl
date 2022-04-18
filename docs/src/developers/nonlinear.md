# [Nonlinear](@id nonlinear_developers)

The `JuMP.Nonlinear` submodule contains data structures and functions for
working with a nonlinear program in the form of an expression tree. This page
explains the API and describes the rationale behind its design.

## Standard form

[Nonlinear programs (NLPs)](https://en.wikipedia.org/wiki/Nonlinear_programming)
are a class of optimization problems in which some of the constraints or the
objective function are nonlinear:
```math
\begin{align}
    \min_{x \in \mathbb{R}^n} & f_0(x) \\
    \;\;\text{s.t.} & l_j \le f_j(x) \le u_j & j = 1 \ldots m
\end{align}
```
There may be additional constraints, as well as things like variable bounds
and integrality restrictions, but we do not consider them here because they are
best dealt with by other components of JuMP and MathOptInterface.

## API overview

The core element of the `Nonlinear` submodule is
[`Nonlinear.NonlinearData`](@ref):
```jldoctest nonlinear_developer
julia> import JuMP: Nonlinear

julia> data = Nonlinear.NonlinearData()
NonlinearData with available features:
  * :ExprGraph
```
[`Nonlinear.NonlinearData`](@ref) is a mutable struct that stores all of the
nonlinear information added to the model.

### Decision variables

Decision variables are represented by [`MOI.VariableIndex`](@ref)s. The user is
responsible for creating these.

### [Expressions](@id Nonlinear_Expressions)

The input data-structure is a Julia `Expr`. The input expressions can
incorporate [`MOI.VariableIndex`](@ref)es, but these must be interpolated into
the expression with `$`:
```jldoctest nonlinear_developer
julia> import JuMP: MOI

julia> x = MOI.VariableIndex(1)
MathOptInterface.VariableIndex(1)

julia> input = :(1 + sin($x)^2)
:(1 + sin(MathOptInterface.VariableIndex(1)) ^ 2)
```
There are a number of restrictions on the input `Expr`:
 * It cannot contain macros
 * It cannot contain broadcasting
 * It cannot contain splatting (except in limited situations)
 * It cannot contain linear algebra, such as matrix-vector products
 * It cannot contain generator expressions, including `sum(i for i in S)`

Given an input expression, add an expression using
[`Nonlinear.add_expression`](@ref):
```jldoctest nonlinear_developer
julia> expr = Nonlinear.add_expression(data, input)
JuMP.Nonlinear.ExpressionIndex(1)
```
The return value, `expr`, is a [`Nonlinear.ExpressionIndex`](@ref) that can
then be interpolated into other input expressions.

### [Parameters](@id Nonlinear_Parameters)

In addition to constant literals like `1` or `1.23`, you can create parameters.
Parameter are constants that you can change before passing the expression to the
solver. Create a parameter using [`Nonlinear.add_parameter`](@ref), which
accepts a default value:
```jldoctest nonlinear_developer
julia> p = Nonlinear.add_parameter(data, 1.23)
JuMP.Nonlinear.ParameterIndex(1)
```
The return value, `p`, is a [`Nonlinear.ParameterIndex`](@ref) that can then be
interpolated into other input expressions.

Update a parameter as follows:
```jldoctest nonlinear_developer
julia> data[p]
1.23

julia> data[p] = 4.56
4.56

julia> data[p]
4.56
```

### [Objectives](@id Nonlinear_Objectives)

Set a nonlinear objective using [`Nonlinear.set_objective`](@ref):
```jldoctest nonlinear_developer
julia> Nonlinear.set_objective(data, :($p + $expr + $x))
```

### [Constraints](@id Nonlinear_Constraints)

Add a constraint using [`Nonlinear.add_constraint`](@ref):
```jldoctest nonlinear_developer
julia> c = Nonlinear.add_constraint(data, :(1 + sqrt($x) <= 2.0))
JuMP.Nonlinear.ConstraintIndex(1)
```
The return value, `c`, is a [`Nonlinear.ConstraintIndex`](@ref) that is a unique
identifier for the constraint. Interval constraints are also supported:
```jldoctest nonlinear_developer
julia> c2 = Nonlinear.add_constraint(data, :(-1.0 <= 1 + sqrt($x) <= 2.0))
JuMP.Nonlinear.ConstraintIndex(2)
```

Delete a constraint using [`Nonlinear.delete`](@ref):
```jldoctest nonlinear_developer
julia> Nonlinear.delete(data, c2)
```

### User-defined operators

By default, `Nonlinear` supports a wide range of univariate and multivariate
operators. However, you can also define your own operators by _registering_
them.

#### Univariate operators

Register a univariate user-defined operator using
[`Nonlinear.register_operator`](@ref):
```jldoctest nonlinear_developer
julia> f(x) = 1 + sin(x)^2
f (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_f, 1, f)
```
Now, you can use `:my_f` in expressions:
```jldoctest nonlinear_developer
julia> new_expr = Nonlinear.add_expression(data, :(my_f($x + 1)))
JuMP.Nonlinear.ExpressionIndex(2)
```
By default, `Nonlinear` will compute first- and second-derivatives of the
registered operator using `ForwardDiff.jl`. Override this by passing functions
which compute the respective derivative:
```jldoctest nonlinear_developer
julia> f′(x) = 2 * sin(x) * cos(x)
f′ (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_f2, 1, f, f′)
```
or
```jldoctest nonlinear_developer
julia> f′′(x) = 2 * (cos(x)^2 - sin(x)^2)
f′′ (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_f3, 1, f, f′, f′′)
```

#### Multivariate operators

Register a multivariate user-defined operator using
[`Nonlinear.register_operator`](@ref):
```jldoctest nonlinear_developer
julia> g(x...) = x[1]^2 + x[1] * x[2] + x[2]^2
g (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_g, 2, g)
```
Now, you can use `:my_f` in expressions:
```jldoctest nonlinear_developer
julia> new_expr = Nonlinear.add_expression(data, :(my_g($x + 1, $x)))
JuMP.Nonlinear.ExpressionIndex(3)
```
By default, `Nonlinear` will compute the gradient of the registered
operator using `ForwardDiff.jl`. (Hessian information is not supported.)
Over-ride this by passing a function to compute the gradient:
```jldoctest nonlinear_developer
julia> function ∇g(ret, x...)
           ret[1] = 2 * x[1] + x[2]
           ret[2] = x[1] + 2 * x[2]
           return
       end
∇g (generic function with 1 method)

julia> Nonlinear.register_operator(data, :my_g2, 2, g, ∇g)
```

### [MathOptInterface](@id Nonlinear_MOI_interface)

`Nonlinear` implements the MathOptInterface API to allow solvers to query the
function and derivative information of our nonlinear model `data`. However,
before we can call [`MOI.initialize`](@ref), we need to set an
[`Nonlinear.AbstractAutomaticDifferentiation`](@ref).

There are two to choose from within JuMP, although other packages may add more
options by sub-typing [`Nonlinear.AbstractAutomaticDifferentiation`](@ref):
 * [`Nonlinear.Default`](@ref)
 * [`Nonlinear.SparseReverseMode`](@ref).

If we set [`Nonlinear.Default`](@ref), then we get access to `:ExprGraph`:
```jldoctest nonlinear_developer
julia> Nonlinear.set_differentiation_backend(data, Nonlinear.Default(), [x])

julia> data
NonlinearData with available features:
  * :ExprGraph
```
!!! note
    [`Nonlinear.set_differentiation_backend`](@ref) requires an ordered list of
    the variables that are included in the model. This order corresponds to the
    the order of the primal decision vector `x` which is passed to the various
    functions in MOI's nonlinear API.

The `:ExprGraph` feature means we can call [`MOI.objective_expr`](@ref) and
[`MOI.constraint_expr`](@ref) to retrieve the expression graph of the problem.
However, we cannot call gradient terms such as
[`MOI.eval_objective_gradient`](@ref) because [`Nonlinear.Default`](@ref) does
not know how to differentiate a nonlinear expression.

If, instead, we set [`Nonlinear.SparseReverseMode`](@ref), then we get access to
`:Grad`, the gradient of the objective function, `:Jac`, the jacobian matrix of
the constraints, `:JacVec`, the ability to compute Jacobian-vector products, and
`:ExprGraph`.
```jldoctest nonlinear_developer
julia> Nonlinear.set_differentiation_backend(
           data,
           Nonlinear.SparseReverseMode(),
           [x],
       )

julia> data
NonlinearData with available features:
  * :Grad
  * :Jac
  * :JacVec
  * :ExprGraph
```

However, before calling anything, we need to call [`MOI.initialize`](@ref):
```jldoctest nonlinear_developer
julia> MOI.initialize(data, [:Grad, :Jac, :JacVec, :ExprGraph])
```

Now we can call methods like [`MOI.eval_objective`](@ref):
```jldoctest nonlinear_developer
julia> x = [1.0]
1-element Vector{Float64}:
 1.0

julia> MOI.eval_objective(data, x)
7.268073418273571
```
and [`MOI.eval_objective_gradient`](@ref):
```jldoctest nonlinear_developer
julia> grad = [NaN]
1-element Vector{Float64}:
 NaN

julia> MOI.eval_objective_gradient(data, grad, x)

julia> grad
1-element Vector{Float64}:
 1.909297426825682
 ```

## Expression-graph representation

[`Nonlinear.NonlinearData`](@ref) stores nonlinear expressions in
[`Nonlinear.NonlinearExpression`](@ref)s. This section explains the design of
the expression graph datastructure in [`Nonlinear.NonlinearExpression`](@ref).

Given a nonlinear function like `f(x) = sin(x)^2 + x`, the first step is to
convert it into [Polish prefix notation](https://en.wikipedia.org/wiki/Polish_notation):
```
f(x, y) = (+ (^ (sin x) 2) x)
```
This format identifies each operator (function), as well as a list of arguments.
Operators can be univariate, like `sin`, or multivariate, like `+`.

A common way of representing Polish prefix notation in code is as follows:
```jldoctest expr_graph
julia> import JuMP: MOI

julia> x = MOI.VariableIndex(1);

julia> struct ExprNode
           op::Symbol
           children::Vector{Union{ExprNode,Float64,MOI.VariableIndex}}
       end

julia> expr = ExprNode(:+, [ExprNode(:^, [ExprNode(:sin, [x]), 2.0]), x]);
```

This datastructure follows our Polish prefix notation very closely, and we can
easily identify the arguments to an operator. However, it has a significant
draw-back: each node in the graph requires a `Vector`, which is heap-allocated
and tracked by Julia's garbage collector (GC). For large JuMP models, we can
expect to have millions of nodes in the expression graph, so this overhead
quickly becomes prohibiative for computation.

An alternative is to record the expression as a linear tape:
```jldoctest expr_graph
julia> expr = Any[:+, 2, :^, 2, :sin, 1, x, 2.0, x]
9-element Vector{Any}:
  :+
 2
  :^
 2
  :sin
 1
  MathOptInterface.VariableIndex(1)
 2.0
  MathOptInterface.VariableIndex(1)
```
The `Int` after each operator `Symbol` specifies the number of arguments.

This data-structure is a single vector, which resolves our problem with the GC,
but each element is the abstract type, `Any`, and so any operations on it will
lead to slower dynamic dispatch. It's also hard to identify the children of each
operation without reading the entire tape.

To summarize, representing expression graphs in Julia has the following
challenges:
 * Nodes in the expression graph should not contain a heap-allocated object
 * All data-structures should be concretely typed
 * It should be easy to identify the children of a node

### Sketch of the design in Nonlinear

`Nonlinear` overcomes these problems by decomposing the datastructure into a
number of different concrete-typed vectors.

First, we create vectors of the supported uni- and multivariate operators.
```jldoctest expr_graph
julia> const UNIVARIATE_OPERATORS = [:sin];

julia> const MULTIVARIATE_OPERATORS = [:+, :^];
```
In practice, there are many more supported operations than the ones listed here.

Second, we create an enum to represent the different types of nodes present in
the expression graph:
```jldoctest expr_graph
julia> @enum(
           NodeType,
           NODE_CALL_MULTIVARIATE,
           NODE_CALL_UNIVARIATE,
           NODE_VARIABLE,
           NODE_VALUE,
       )
```
In practice, there are node types other than the ones listed here.

Third, we create two concretely-typed structs as follows:
```jldoctest expr_graph
julia> struct Node
           type::NodeType
           parent::Int
           index::Int
       end

julia> struct NonlinearExpression
           nodes::Vector{Node}
           values::Vector{Float64}
       end
```

For each node `node` in the `.nodes` field, if `node.type` is:
 * `NODE_CALL_MULTIVARIATE`, we look up
   `MULTIVARIATE_OPERATORS[node.index]` to retrieve the operator
 * `NODE_CALL_UNIVARIATE`, we look up
   `UNIVARIATE_OPERATORS[node.index]` to retrieve the operator
 * `NODE_VARIABLE`, we create `MOI.VariableIndex(node.index)`
 * `NODE_VALUE`, we look up `values[node.index]`
The `.parent` field of each node is the integer index of the parent node in
`.nodes`. For the first node, the parent is `-1` by convention.

Therefore, we can represent our function as:
```jldoctest expr_graph
julia> expr = NonlinearExpression(
           [
               Node(NODE_CALL_MULTIVARIATE, 1, -1),
               Node(NODE_CALL_MULTIVARIATE, 2, 1),
               Node(NODE_CALL_UNIVARIATE, 1, 2),
               Node(NODE_VARIABLE, 1, 3),
               Node(NODE_VALUE, 1, 2),
               Node(NODE_VARIABLE, 1, 1),
           ],
           [2.0],
       );
```

This is less readable than the other options, but does this datastructure meet
our design goals?

Instead of a heap-allocated object for each node, we only have two `Vector`s for
each expression, `nodes` and `values`, as well as two constant vectors for the
`OPERATORS`. In addition, all fields are concretely typed, and there are no
`Union` or `Any` tyypes.

For our third goal, it is not easy to identify the children of a node, but it is
easy to identify the _parent_ of any node. Therefore, we can use
[`Nonlinear.adjacency_matrix`](@ref) to compute a sparse matrix that maps
children to their parents.

The tape is also ordered topologically, so that a reverse pass of the nodes
evaluates all children nodes before their parent.

### The design in practice

In practice, `Node` and `NonlinearExpression` are exactly [`Nonlinear.Node`](@ref)
and [`Nonlinear.NonlinearExpression`](@ref). However, [`Nonlinear.NodeType`](@ref)
has more terms to account for comparison operators such as `:>=` and `:<=`,
logic operators such as `:&&` and `:||`, nonlinear parameters, and nested
subexpressions.

Moreover, instead of storing the operators as global constants, they are stored
in [`Nonlinear.OperatorRegistry`](@ref), and it also stores a vector of logic
operators and a vector of comparison operators. In addition to
[`Nonlinear.DEFAULT_UNIVARIATE_OPERATORS`](@ref) and
[`Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS`](@ref), you can register
user-defined functions using [`Nonlinear.register_operator`](@ref).

[`Nonlinear.NonlinearData`](@ref) is a struct that stores the
[`Nonlinear.OperatorRegistry`](@ref), as well as a list of parameters and
subexpressions in the model.

## ReverseAD

`Nonlinear.ReverseAD` is a submodule for computing derivatives of the problem
inside [`Nonlinear.NonlinearData`](@ref) using sparse reverse-mode automatic
differentiation (AD).

This section does not attempt to explain how sparse reverse-mode AD works, but
instead explains why JuMP contains it's own implementation, and highlights
notable differences from similar packages.

!!! warning
    You should not interact with `ReverseAD` directly. Instead, you should
    create a [`Nonlinear.NonlinearData`](@ref) object, call
    [`Nonlinear.set_differentiation_backend`](@ref) with
    [`Nonlinear.SparseReverseMode`](@ref), and then query the MOI API methods.

### Why another AD package?

The JuliaDiff organization maintains a [list of packages](https://juliadiff.org)
for doing AD in Julia. At last count, there were at least ten packages–not
including `ReverseAD`–for reverse-mode AD in Julia. Given this multitude, why
does JuMP maintain another implementation instead of re-using existing tooling?

Here are four reasons:

 * **Scale and Sparsity:** the types of functions that JuMP computes derivatives
   of have two key characteristics: they can be very large scale (10^5 or more
   functions across 10^5 or more variables) and they are very sparse. For large
   problems, it is common for the hessian to have `O(n)` non-zero elements
   instead of `O(n^2)` if it was dense. To the best of our knowledge,
   `ReverseAD` is the only reverse-mode AD system in Julia that handles sparsity
   by default. The lack of sparsity support is _the_ main reason why we do no
   use a generic package.
 * **Limited scope:** most other AD packages accept arbitrary Julia functions as
   input and then trace an expression graph using operator overloading. This
   means they must deal (or detect and ignore) with control flow, I/O, and other
   vagaries of Julia. In contrast, `ReverseAD` only accepts functions in the
   form of [`Nonlinear.NonlinearExpression`](@ref), which greatly limits the
   range of syntax that it must deal with. By reducing the scope of what we
   accept as input to functions relevant for mathematical optimization, we can
   provide a simpler implementation with various performance optimizations.
 * **Historical:** `ReverseAD` started life as [ReverseDiffSparse.jl](https://github.com/mlubin/ReverseDiffSparse.jl),
   development of which begain in early 2014(!). This was well before the other
   packages started development. Because we had a well-tested, working AD in
   JuMP, there was less motivation to contribute to and explore other AD
   packages. The lack of historical interaction also meant that other packages
   were not optimized for the types of problems that JuMP is built for (i.e.,
   large-scale sparse problems).
 * **Technical debt** Prior to the introduction of `Nonlinear`, JuMP's nonlinear
   implementation was a confusing mix of functions and types spread across the
   code base and in the private `_Derivatives` submodule. This made it hard to
   swap the AD system for another. The main motivation for refactoring JuMP to
   create the `Nonlinear` submodule was to abstract the interface between JuMP
   and the AD system, allowing us to swap-in and test new AD systems in the
   future.
