```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Extensions](@id extensions_manual)

JuMP provides a variety of ways to extend the basic modeling functionality.

!!! tip
    This documentation in this section is still a work-in-progress. The best
    place to look for ideas and help when writing a new JuMP extension are
    existing JuMP extensions. Examples include:
     * [BilevelJuMP.jl](https://github.com/joaquimg/BilevelJuMP.jl)
     * [Coluna.jl](https://github.com/atoptima/Coluna.jl)
     * [InfiniteOpt.jl](https://github.com/infiniteopt/InfiniteOpt.jl)
     * [Plasmo.jl](https://github.com/zavalab/Plasmo.jl)
     * [PolyJuMP.jl](https://github.com/jump-dev/PolyJuMP.jl)
     * [SDDP.jl](https://github.com/odow/SDDP.jl)
     * [StochasticPrograms.jl](https://github.com/martinbiel/StochasticPrograms.jl)
     * [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl)
     * [vOptGeneric.jl](https://github.com/vOptSolver/vOptGeneric.jl)

## Compatibility

When writing JuMP extensions, you should carefully consider the compatibility
guarantees that JuMP makes. In particular:

 * All functions, structs, and constants which do not begin with an underscore
   (`_`) are public. These are always safe to use, and they should all have
   corresponding documentation.
 * All identifiers beginning with an underscore (`_`) are private. These are
   not safe to use, because they may break in any JuMP release, including
   patch releases.
 * Unless explicitly mentioned in the documentation, all fields of a struct are
   private. These are not safe to use, because they may break in any JuMP
   release, including patch releases. An example of a field which is safe to use
   is the `model.ext` extension dictionary, which is documented in
   [The extension dictionary](@ref).

In general, we strongly encourage you to use only the public API of JuMP. If you
are missing a feature, please open a GitHub issue.

However, if you _do_ use the private API (for example, because your feature
request has not been implemented yet), then you must carefully restrict the
versions of JuMP that your package is compatible with in the `Project.toml`
file. The easiest way to do this is via the [hyphen specifiers](https://pkgdocs.julialang.org/v1/compatibility/).
For example, if your package supports all JuMP versions between v1.0.0 and
v1.1.1, do:
```
JuMP = "1.0.0 - 1.1.1"
```
Then, whenever JuMP releases a new version, you should check if your package is
still compatible and update the bound accordingly.

## Define a new set

To define a new set for JuMP, subtype `MOI.AbstractScalarSet` or
`MOI.AbstractVectorSet` and implement `Base.copy` for the set.

```jldoctest define_new_set
julia> struct NewMOIVectorSet <: MOI.AbstractVectorSet
           dimension::Int
       end

julia> Base.copy(x::NewMOIVectorSet) = x

julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @constraint(model, x in NewMOIVectorSet(2))
[x[1], x[2]] ∈ NewMOIVectorSet(2)
```

However, for vector-sets, this requires the user to specify the dimension
argument to their set, even though we could infer it from the length of `x`!

You can make a more user-friendly set by subtyping [`AbstractVectorSet`](@ref)
and implementing [`moi_set`](@ref).

```jldoctest define_new_set
julia> struct NewVectorSet <: JuMP.AbstractVectorSet end

julia> JuMP.moi_set(::NewVectorSet, dim::Int) = NewMOIVectorSet(dim)

julia> @constraint(model, x in NewVectorSet())
[x[1], x[2]] ∈ NewMOIVectorSet(2)
```

## [Extend `@variable`](@id extend_variable_macro)

Just as `Bin` and `Int` create binary and integer variables, you can extend
the [`@variable`](@ref) macro to create new types of variables. Here is an
explanation by example, where we create a `AddTwice` type, that creates a tuple
of two JuMP variables instead of a single variable.

First, create a new struct. This can be anything. Our struct holds a
[`VariableInfo`](@ref) object that stores bound information, and whether the
variable is binary or integer.
```jldoctest new_variable
julia> struct AddTwice
           info::JuMP.VariableInfo
       end
```

Second, implement [`build_variable`](@ref), which takes `::Type{AddTwice}` as
an argument, and returns an instance of `AddTwice`. Note that you can also
receive keyword arguments.
```jldoctest new_variable
julia> function JuMP.build_variable(
           _err::Function,
           info::JuMP.VariableInfo,
           ::Type{AddTwice};
           kwargs...
       )
           println("Can also use $kwargs here.")
           return AddTwice(info)
       end
```

Third, implement [`add_variable`](@ref), which takes the instance of `AddTwice`
from the previous step, and returns something. Typically, you will want to call
[`add_variable`](@ref) here. For example, our `AddTwice` call is going to add
two JuMP variables.
```jldoctest new_variable
julia> function JuMP.add_variable(
           model::JuMP.Model,
           duplicate::AddTwice,
           name::String,
       )
           a = JuMP.add_variable(
               model,
               JuMP.ScalarVariable(duplicate.info),
               "$(name)_a",
            )
           b = JuMP.add_variable(
               model,
               JuMP.ScalarVariable(duplicate.info),
               "$(name)_b",
            )
           return (a, b)
       end
```

Now `AddTwice` can be passed to [`@variable`](@ref) similar to `Bin` or `Int`,
or through the `variable_type` keyword. However, now it adds two variables
instead of one.

```jldoctest new_variable
julia> model = Model();

julia> @variable(model, x[i=1:2], variable_type = AddTwice, kw = i)
Can also use Base.Pairs(:kw => 1) here.
Can also use Base.Pairs(:kw => 2) here.
2-element Vector{Tuple{VariableRef, VariableRef}}:
 (x[1]_a, x[1]_b)
 (x[2]_a, x[2]_b)

julia> num_variables(model)
4

julia> first(x[1])
x[1]_a

julia> last(x[2])
x[2]_b
```

## Extend [`@constraint`](@ref)

The [`@constraint`](@ref) macro has three steps that can be intercepted and
extended: parse time, build time, and add time.

### Parse

To extend the [`@constraint`](@ref) macro at parse time, implement one of the
following methods:

 * [`parse_constraint_head`](@ref)
 * [`parse_constraint_call`](@ref)

!!! warning
    Extending the constraint macro at parse time is an advanced operation and
    has the potential to interfere with existing JuMP syntax. Please discuss
    with the [developer chatroom](https://gitter.im/JuliaOpt/jump-dev) before
    publishing any code that implements these methods.

[`parse_constraint_head`](@ref) should be implemented to intercept an expression
based on the `.head` field of `Base.Expr`. For example:
```jldoctest
julia> using JuMP

julia> const MutableArithmetics = JuMP._MA;

julia> model = Model(); @variable(model, x);

julia> function JuMP.parse_constraint_head(
           error_fn::Function,
           ::Val{:≔},
           lhs,
           rhs,
       )
           println("Rewriting ≔ as ==")
           new_lhs, parse_code = MutableArithmetics.rewrite(lhs)
           build_code = :(
               build_constraint($(error_fn), $(new_lhs), MOI.EqualTo($(rhs)))
           )
           return false, parse_code, build_code
       end

julia> @constraint(model, x + x ≔ 1.0)
Rewriting ≔ as ==
2 x = 1
```

[`parse_constraint_call`](@ref) should be implemented to intercept an expression
of the form `Expr(:call, op, args...)`. For example:
```jldoctest
julia> using JuMP

julia> const MutableArithmetics = JuMP._MA;

julia> model = Model(); @variable(model, x);

julia> function JuMP.parse_constraint_call(
           error_fn::Function,
           is_vectorized::Bool,
           ::Val{:my_equal_to},
           lhs,
           rhs,
       )
           println("Rewriting my_equal_to to ==")
           new_lhs, parse_code = MutableArithmetics.rewrite(lhs)
           build_code = if is_vectorized
               :(build_constraint($(error_fn), $(new_lhs), MOI.EqualTo($(rhs)))
           )
           else
               :(build_constraint.($(error_fn), $(new_lhs), MOI.EqualTo($(rhs))))
           end
           return parse_code, build_code
       end

julia> @constraint(model, my_equal_to(x + x, 1.0))
Rewriting my_equal_to to ==
2 x = 1
```

!!! tip
    When parsing a constraint you can recurse into sub-constraint (for example, the
    `{expr}` in `z --> {x <= 1}`) by calling [`parse_constraint`](@ref).

To prevent JuMP from promoting the set to the same value type as the model, use
[`SkipModelConvertScalarSetWrapper`](@ref).

### Build

To extend the [`@constraint`](@ref) macro at build time, implement a new
[`build_constraint`](@ref) method.

This may mean implementing a method for a specific function or set created at
parse time, or it may mean implementing a method which handles additional
positional arguments.

[`build_constraint`](@ref) must return an [`AbstractConstraint`](@ref), which
can either be an [`AbstractConstraint`](@ref) already supported by JuMP, for example, `ScalarConstraint` or `VectorConstraint`, or a custom
[`AbstractConstraint`](@ref) with a corresponding [`add_constraint`](@ref)
method (see [Add](@ref extension_add_constraint)).

!!! tip
    The easiest way to extend [`@constraint`](@ref) is via an additional
    positional argument to [`build_constraint`](@ref).

Here is an example of adding extra arguments to [`build_constraint`](@ref):
```jldoctest
julia> model = Model(); @variable(model, x);

julia> struct MyConstrType end

julia> function JuMP.build_constraint(
            error_fn::Function,
            f::JuMP.GenericAffExpr,
            set::MOI.EqualTo,
            extra::Type{MyConstrType};
            d = 0,
       )
            new_set = MOI.LessThan(set.value + d)
            return JuMP.build_constraint(error_fn, f, new_set)
       end

julia> @constraint(model, my_con, x == 0, MyConstrType, d = 2)
my_con : x ≤ 2
```

!!! note
    Only a single positional argument can be given to a particular constraint.
    Extensions that seek to pass multiple arguments (for example, `Foo` and `Bar`)
    should combine them into one argument type (for example, `FooBar`).

### [Add](@id extension_add_constraint)

[`build_constraint`](@ref) returns an [`AbstractConstraint`](@ref) object. To
extend [`@constraint`](@ref) at add time, define a subtype of
[`AbstractConstraint`](@ref), implement [`build_constraint`](@ref) to return an
instance of the new type, and then implement [`add_constraint`](@ref).

Here is an example:
```jldoctest
julia> model = Model(); @variable(model, x);

julia> struct MyTag
           name::String
       end

julia> struct MyConstraint{S} <: AbstractConstraint
           name::String
           f::AffExpr
           s::S
       end

julia> function JuMP.build_constraint(
            error_fn::Function,
            f::AffExpr,
            set::MOI.AbstractScalarSet,
            extra::MyTag,
       )
            return MyConstraint(extra.name, f, set)
       end

julia> function JuMP.add_constraint(
            model::Model,
            con::MyConstraint,
            name::String,
       )
            return add_constraint(
                model,
                ScalarConstraint(con.f, con.s),
                "$(con.name)[$(name)]",
            )
       end

julia> @constraint(model, my_con, 2x <= 1, MyTag("my_prefix"))
my_prefix[my_con] : 2 x - 1 ≤ 0
```

## The extension dictionary

Every JuMP model has a field `.ext::Dict{Symbol,Any}` that can be used by
extensions. This is useful if your extensions to [`@variable`](@ref) and
[`@constraint`](@ref) need to store information between calls.

The most common way to initialize a model with information in the `.ext`
dictionary is to provide a new constructor:
```jldoctest
julia> function MyModel()
           model = Model()
           model.ext[:MyModel] = 1
           return model
       end
MyModel (generic function with 1 method)

julia> model = MyModel()
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none

julia> model.ext
Dict{Symbol, Any} with 1 entry:
  :MyModel => 1
```

If you define extension data, implement [`copy_extension_data`](@ref)
to support [`copy_model`](@ref).

## Defining new JuMP models

If extending individual calls to [`@variable`](@ref) and [`@constraint`](@ref)
is not sufficient, it is possible to implement a new model via a subtype of
[`AbstractModel`](@ref). You can also define new [`AbstractVariableRef`](@ref)s
to create different types of JuMP variables.

!!! warning
    Extending JuMP in this manner is an advanced operation. We strongly
    encourage you to consider how you can use the methods mentioned in the
    previous sections to achieve your aims instead of defining new model and
    variable types. Consult the [developer chatroom](https://gitter.im/JuliaOpt/jump-dev)
    _before_ starting work on this.

If you define new types, you will need to implement a considerable number of
methods, and doing so will require a detailed understanding of the JuMP
internals. Therefore, the list of methods to implement is currently
undocumented.

The easiest way to extend JuMP by defining a new model type is to follow an
existing example. A simple example to follow is the [JuMPExtension module](https://github.com/jump-dev/JuMP.jl/blob/master/test/JuMPExtension.jl)
in the JuMP test suite. The best example of an external JuMP extension that
implements an [`AbstractModel`](@ref) is [InfiniteOpt.jl](https://github.com/infiniteopt/InfiniteOpt.jl).

### Testing JuMP extensions

The JuMP test suite contains a large number of tests for JuMP extensions. You
can run these tests by copying the MIT-licensed [`Kokako.jl`](https://github.com/jump-dev/JuMP.jl/blob/master/test/Kokako.jl)
file in the JuMP tests into your `/test` folder, and then adding this snippet to
your `/test/runtests.jl` file:

```julia
using MyJuMPExtension
import JuMP
include("Kokako.jl")
const MODULES_TO_TEST = Kokako.include_modules_to_test(JuMP)
Kokako.run_tests(
    MODULES_TO_TEST,
    MyJuMPExtension.MyModel,
    MyJuMPExtension.MyVariableRef;
    test_prefix = "test_extension_",
)
```

## Set an `optimize!` hook

Some extensions require modification to the problem after the user has finished
constructing the problem, but before `optimize!` is called. For these
situations, JuMP provides [`set_optimize_hook`](@ref), which lets you intercept
the [`optimize!`](@ref) call.

Here's a simple example of adding an optimize hook that extends [`optimize!`](@ref)
to take a keyword argument `silent`:
```jldoctest
julia> using JuMP, HiGHS

julia> model = Model(HiGHS.Optimizer);

julia> @variable(model, x >= 1.5, Int);

julia> @objective(model, Min, x);

julia> function silent_hook(model; silent::Bool)
           if silent
               set_silent(model)
           else
               unset_silent(model)
           end
           ## Make sure you set ignore_optimize_hook = true, or we'll
           ## recursively enter the optimize hook!
           return optimize!(model; ignore_optimize_hook = true)
       end
silent_hook (generic function with 1 method)

julia> set_optimize_hook(model, silent_hook)
silent_hook (generic function with 1 method)

julia> optimize!(model; silent = true)

julia> optimize!(model; silent = false)
Coefficient ranges:
  Cost   [1e+00, 1e+00]
  Bound  [2e+00, 2e+00]
Assessing feasibility of MIP using primal feasibility and integrality tolerance of       1e-06
Solution has               num          max          sum
Col     infeasibilities      0            0            0
Integer infeasibilities      0            0            0
Row     infeasibilities      0            0            0
Row     residuals            0            0            0
Presolving model
0 rows, 0 cols, 0 nonzeros  0s
0 rows, 0 cols, 0 nonzeros  0s
Presolve: Optimal

Solving report
  Status            Optimal
  Primal bound      2
  Dual bound        2
  Gap               0% (tolerance: 0.01%)
  Solution status   feasible
                    2 (objective)
                    0 (bound viol.)
                    0 (int. viol.)
                    0 (row viol.)
  Timing            0.00 (total)
                    0.00 (presolve)
                    0.00 (postsolve)
  Nodes             0
  LP iterations     0 (total)
                    0 (strong br.)
                    0 (separation)
                    0 (heuristics)
```

## Creating new container types

JuMP macros (for example, [`@variable`](@ref)) accept a `container` keyword
argument to force the type of container that is chosen. By default, JuMP
supports `container = Array`, `container = DenseAxisArray`,
`container = SparseAxisArray` and `container = Auto`. You can extend support to
user-defined types by implementing [`Containers.container`](@ref).

For example, here is a container that reverses the order of the indices:
```jldoctest extend_containers
julia> struct Foo end

julia> function Containers.container(f::Function, indices, ::Type{Foo})
           return reverse([f(i...) for i in indices])
       end

julia> model = Model();

julia> @variable(model, x[1:3], container = Foo)
3-element Vector{VariableRef}:
 x[3]
 x[2]
 x[1]

julia> x[1]
x[3]

julia> @variable(model, y[1:3, 1:2], container = Foo)
3×2 Matrix{VariableRef}:
 y[3,2]  y[3,1]
 y[2,2]  y[2,1]
 y[1,2]  y[1,1]

julia> y[1, 1]
y[3,2]

julia> @variable(model, z[i=1:3; isodd(i)], container = Foo)
2-element Vector{VariableRef}:
 z[3]
 z[1]

julia> z[2]
z[1]
```

!!! warning
    If you are a general user, you should not need to create a new container
    type. Instead, consider following [User-defined containers](@ref) and create
    a new container using standard Julia syntax. For example:
    ```jldoctest
    julia> model = Model();

    julia> @variable(model, x[1:3])
    3-element Vector{VariableRef}:
     x[1]
     x[2]
     x[3]

    julia> y = reverse(x)
    3-element Vector{VariableRef}:
     x[3]
     x[2]
     x[1]
    ```

## Performance tips for extensions

The function-in-set design of MathOptInterface causes type stability issues in
Julia if you try to iterate over all of the constraints in a model. The easiest
way to fix this is to use a function barrier.

For example, instead of:
```julia
function all_names_slow(model)
    names = Set{String}()
    for ci in all_constraints(model)
        push!(names, name(ci))
    end
    return names
end
```
use:
```julia
function _function_barrier(names, model, ::Type{F}, ::Type{S}) where {F,S}
    for ci in all_constraints(model, F, S)
        push!(names, name(ci))
    end
    return
end

function all_names_fast(model)
    names = Set{String}()
    for (F, S) in list_of_constraint_types(model)
        _function_barrier(names, model, F, S)
    end
    return names
end
```

!!! note
    It is important to explicitly type the `F` and `S` arguments. If you leave
    them untyped, for example, `function _function_barrier(names, model, F, S)`,
    Julia will not specialize the function calls and performance will not be
    improved.
