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
     * [InfiniteOpt.jl](https://github.com/pulsipher/InfiniteOpt.jl)
     * [Plasmo.jl](https://github.com/zavalab/Plasmo.jl)
     * [PolyJuMP.jl](https://github.com/jump-dev/PolyJuMP.jl)
     * [SDDP.jl](https://github.com/odow/SDDP.jl)
     * [StochasticPrograms.jl](https://github.com/martinbiel/StochasticPrograms.jl)
     * [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl)
     * [vOptGeneric.jl](https://github.com/vOptSolver/vOptGeneric.jl)

## Define a new set

To define a new set for JuMP, subtype `MOI.AbstractScalarSet` or
`MOI.AbstractVectorSet` and implement `Base.copy` for the set. That's it!

```jldoctest define_new_set
struct _NewVectorSet <: MOI.AbstractVectorSet
    dimension::Int
end
Base.copy(x::_NewVectorSet) = x

model = Model()
@variable(model, x[1:2])
@constraint(model, x in _NewVectorSet(2))

# output

[x[1], x[2]] ∈ _NewVectorSet(2)
```

However, for vector-sets, this requires the user to specify the dimension
argument to their set, even though we could infer it from the length of `x`!

You can make a more user-friendly set by subtyping [`AbstractVectorSet`](@ref)
and implementing [`moi_set`](@ref).

```jldoctest define_new_set
struct NewVectorSet <: JuMP.AbstractVectorSet end
JuMP.moi_set(::NewVectorSet, dim::Int) = _NewVectorSet(dim)

model = Model()
@variable(model, x[1:2])
@constraint(model, x in NewVectorSet())

# output

[x[1], x[2]] ∈ _NewVectorSet(2)
```

## Extend [`@variable`](@ref)

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
               name * "_a",
            )
           b = JuMP.add_variable(
               model,
               JuMP.ScalarVariable(duplicate.info),
               name * "_b",
            )
           return (a, b)
       end
```

Now `AddTwice` can be passed to [`@variable`](@ref) just like `Bin` or `Int`.
However, now it adds two variables instead of one!
```jldoctest new_variable
julia> model = Model();

julia> @variable(model, x[i=1:2], AddTwice, kw=i)
Can also use Base.Iterators.Pairs(:kw => 1) here.
Can also use Base.Iterators.Pairs(:kw => 2) here.
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
           _error::Function,
           ::Val{:(:=)},
           lhs,
           rhs,
       )
           println("Rewriting := as ==")
           new_lhs, parse_code = MutableArithmetics.rewrite(lhs)
           build_code = :(
               build_constraint($(_error), $(new_lhs), MOI.EqualTo($(rhs)))
           )
           return false, parse_code, build_code
       end

julia> @constraint(model, x + x := 1.0)
Rewriting := as ==
2 x = 1.0
```

[`parse_constraint_call`](@ref) should be implemented to intercept an expression
of the form `Expr(:call, op, args...)`. For example:
```jldoctest
julia> using JuMP

julia> const MutableArithmetics = JuMP._MA;

julia> model = Model(); @variable(model, x);

julia> function JuMP.parse_constraint_call(
           _error::Function,
           is_vectorized::Bool,
           ::Val{:my_equal_to},
           lhs,
           rhs,
       )
           println("Rewriting my_equal_to to ==")
           new_lhs, parse_code = MutableArithmetics.rewrite(lhs)
           build_code = if is_vectorized
               :(build_constraint($(_error), $(new_lhs), MOI.EqualTo($(rhs)))
           )
           else
               :(build_constraint.($(_error), $(new_lhs), MOI.EqualTo($(rhs))))
           end
           return parse_code, build_code
       end

julia> @constraint(model, my_equal_to(x + x, 1.0))
Rewriting my_equal_to to ==
2 x = 1.0
```

!!! tip
    When parsing a constraint you can recurse into sub-constraint (e.g., the
    `{expr}` in `z => {x <= 1}`) by calling [`parse_constraint`](@ref).

### Build

To extend the [`@constraint`](@ref) macro at build time, implement a new
[`build_constraint`](@ref) method.

This may mean implementing a method for a specific function or set created at
parse time, or it may mean implementing a method which handles additional
positional arguments.

[`build_constraint`](@ref) must return an [`AbstractConstraint`](@ref), which
can either be an [`AbstractConstraint`](@ref) already supported by JuMP, e.g., `ScalarConstraint` or `VectorConstraint`, or a custom
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
            _error::Function,
            f::JuMP.GenericAffExpr,
            set::MOI.EqualTo,
            extra::Type{MyConstrType};
            d = 0,
       )
            new_set = MOI.LessThan(set.value + d)
            return JuMP.build_constraint(_error, f, new_set)
       end

julia> @constraint(model, my_con, x == 0, MyConstrType, d = 2)
my_con : x ≤ 2.0
```

!!! note
    Only a single positional argument can be given to a particular constraint.
    Extensions that seek to pass multiple arguments (e.g., `Foo` and `Bar`)
    should combine them into one argument type (e.g., `FooBar`).

### [Add](@id extension_add_constraint)

Work in progress.

### Shapes

Shapes allow vector constraints, which are represented as flat vectors in MOI,
to retain a matrix shape at the JuMP level. There is a `shape` field in
`VectorConstraint` that can be set in [`build_constraint`](@ref) and that is
used to reshape the result computed in [`value`](@ref) and [`dual`](@ref).

## Extend [`@objective`](@ref)

Work in progress.

### Adding a bridge

Work in progress.

## Defining new JuMP models

Work in progress.
