```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Extensions](@id extensions_manual)

JuMP provides a variety of ways to extend the basic modeling functionality.

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
Can also use Base.Iterators.Pairs(:kw=>1) here.
Can also use Base.Iterators.Pairs(:kw=>2) here.
2-element Array{Tuple{VariableRef,VariableRef},1}:
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

The [`@constraint`](@ref) macro always calls the same three functions:
* `parse_constraint`: is called at parsing time, it parses the constraint
  expression and returns a [`build_constraint`](@ref) call expression;
* [`build_constraint`](@ref): given the functions and sets involved in the
  constraints, it returns a `AbstractConstraint`;
* [`add_constraint`](@ref): given the model, the `AbstractConstraint`
  constructed in [`build_constraint`](@ref) and the constraint name, it stores
  them in the model and returns a `ConstraintRef`.

Adding methods to these functions is the recommended way to extend the
[`@constraint`](@ref) macro.

### Adding `parse_constraint` methods

Work in progress.
### Adding `build_constraint` methods

There are typically two choices when creating a [`build_constraint`](@ref)
method, either return an `AbstractConstraint` already supported by the
model, i.e. `ScalarConstraint` or `VectorConstraint`, or a custom
`AbstractConstraint` with a corresponding [`add_constraint`](@ref) method (see
[Adding `add_constraint` methods](@ref)).

### Adding `add_constraint` methods

Work in progress.

### Adding an extra positional argument

We can also extend `@constraint` to handle additional positional arguments that 
effectively "tag" a particular constraint type and/or pass along additional 
information that we may want. For example, we can make a `MyConstrType` that 
modifies affine equalities:
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
Note that only a single positional argument can be given to a particular 
constraint. Extensions that seek to pass multiple arguments (e.g., `Foo` and 
`Bar`) should combine them into one argument type (e.g., `FooBar`). 

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
