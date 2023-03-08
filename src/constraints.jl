#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

"""
    ConstraintRef

Holds a reference to the model and the corresponding MOI.ConstraintIndex.
"""
struct ConstraintRef{M<:AbstractModel,C,Shape<:AbstractShape}
    model::M
    index::C
    shape::Shape
end

"""
    index(cr::ConstraintRef)::MOI.ConstraintIndex

Return the index of the constraint that corresponds to `cr` in the MOI backend.
"""
index(cr::ConstraintRef) = cr.index

"""
    struct ConstraintNotOwned{C <: ConstraintRef} <: Exception
        constraint_ref::C
    end

The constraint `constraint_ref` was used in a model different to
`owner_model(constraint_ref)`.
"""
struct ConstraintNotOwned{C<:ConstraintRef} <: Exception
    constraint_ref::C
end

"""
    owner_model(con_ref::ConstraintRef)

Returns the model to which `con_ref` belongs.
"""
owner_model(con_ref::ConstraintRef) = con_ref.model

"""
    check_belongs_to_model(con_ref::ConstraintRef, model::AbstractModel)

Throw `ConstraintNotOwned` if `owner_model(con_ref)` is not `model`.
"""
function check_belongs_to_model(con_ref::ConstraintRef, model::AbstractModel)
    if owner_model(con_ref) !== model
        throw(ConstraintNotOwned(con_ref))
    end
end

Base.broadcastable(con_ref::ConstraintRef) = Ref(con_ref)

"""
    dual_start_value(con_ref::ConstraintRef)

Return the dual start value (MOI attribute `ConstraintDualStart`) of the
constraint `con_ref`.

Note: If no dual start value has been set, `dual_start_value` will return
`nothing`.

See also [`set_dual_start_value`](@ref).
"""
function dual_start_value(
    con_ref::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex},
)
    return reshape_vector(_dual_start(con_ref), dual_shape(con_ref.shape))
end

# Returns the value of MOI.ConstraintDualStart in a type-stable way
function _dual_start(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractScalarFunction,
            <:MOI.AbstractScalarSet,
        },
    },
)::Union{Nothing,Float64}
    return MOI.get(owner_model(con_ref), MOI.ConstraintDualStart(), con_ref)
end
function _dual_start(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractVectorFunction,
            <:MOI.AbstractVectorSet,
        },
    },
)::Union{Nothing,Vector{Float64}}
    return MOI.get(owner_model(con_ref), MOI.ConstraintDualStart(), con_ref)
end

"""
    set_dual_start_value(con_ref::ConstraintRef, value)

Set the dual start value (MOI attribute `ConstraintDualStart`) of the constraint
`con_ref` to `value`. To remove a dual start value set it to `nothing`.

See also [`dual_start_value`](@ref).
"""
function set_dual_start_value(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractVectorFunction,
            <:MOI.AbstractVectorSet,
        },
    },
    value,
)
    vectorized_value = vectorize(value, dual_shape(con_ref.shape))
    MOI.set(
        owner_model(con_ref),
        MOI.ConstraintDualStart(),
        con_ref,
        vectorized_value,
    )
    return
end
function set_dual_start_value(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractVectorFunction,
            <:MOI.AbstractVectorSet,
        },
    },
    ::Nothing,
)
    MOI.set(owner_model(con_ref), MOI.ConstraintDualStart(), con_ref, nothing)
    return
end
function set_dual_start_value(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractScalarFunction,
            <:MOI.AbstractScalarSet,
        },
    },
    value,
)
    MOI.set(owner_model(con_ref), MOI.ConstraintDualStart(), con_ref, value)
    return
end

"""
    set_start_value(con_ref::ConstraintRef, value)

Set the primal start value ([`MOI.ConstraintPrimalStart`](@ref)) of the
constraint `con_ref` to `value`. To remove a primal start value set it to
`nothing`.

See also [`start_value`](@ref).
"""
function set_start_value(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractVectorFunction,
            <:MOI.AbstractVectorSet,
        },
    },
    value,
)
    MOI.set(
        owner_model(con_ref),
        MOI.ConstraintPrimalStart(),
        con_ref,
        vectorize(value, con_ref.shape),
    )
    return
end

function set_start_value(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractVectorFunction,
            <:MOI.AbstractVectorSet,
        },
    },
    ::Nothing,
)
    MOI.set(owner_model(con_ref), MOI.ConstraintPrimalStart(), con_ref, nothing)
    return
end

function set_start_value(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractScalarFunction,
            <:MOI.AbstractScalarSet,
        },
    },
    value,
)
    MOI.set(owner_model(con_ref), MOI.ConstraintPrimalStart(), con_ref, value)
    return
end

"""
    start_value(con_ref::ConstraintRef)

Return the primal start value ([`MOI.ConstraintPrimalStart`](@ref)) of the
constraint `con_ref`.

Note: If no primal start value has been set, `start_value` will return
`nothing`.

See also [`set_start_value`](@ref).
"""
function start_value(
    con_ref::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex},
)
    return reshape_vector(
        MOI.get(owner_model(con_ref), MOI.ConstraintPrimalStart(), con_ref),
        con_ref.shape,
    )
end

"""
    name(con_ref::ConstraintRef)

Get a constraint's name attribute.
"""
function name(
    con_ref::ConstraintRef{<:AbstractModel,C},
) where {C<:MOI.ConstraintIndex}
    model = owner_model(con_ref)
    if !MOI.supports(backend(model), MOI.ConstraintName(), C)
        return ""
    end
    return MOI.get(model, MOI.ConstraintName(), con_ref)::String
end

# The name of VariableIndex constraints is empty.
function name(
    ::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex{MOI.VariableIndex}},
)
    return ""
end

"""
    set_name(con_ref::ConstraintRef, s::AbstractString)

Set a constraint's name attribute.
"""
function set_name(
    con_ref::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex},
    s::String,
)
    return MOI.set(con_ref.model, MOI.ConstraintName(), con_ref, s)
end

"""
    constraint_by_name(model::AbstractModel,
                       name::String)::Union{ConstraintRef, Nothing}

Return the reference of the constraint with name attribute `name` or `Nothing`
if no constraint has this name attribute. Throws an error if several
constraints have `name` as their name attribute.

    constraint_by_name(model::AbstractModel,
                       name::String,
                       F::Type{<:Union{AbstractJuMPScalar,
                                       Vector{<:AbstractJuMPScalar},
                                       MOI.AbstactFunction}},
                       S::Type{<:MOI.AbstractSet})::Union{ConstraintRef, Nothing}

Similar to the method above, except that it throws an error if the constraint is
not an `F`-in-`S` contraint where `F` is either the JuMP or MOI type of the
function, and `S` is the MOI type of the set. This method is recommended if you
know the type of the function and set since its returned type can be inferred
while for the method above (i.e. without `F` and `S`), the exact return type of
the constraint index cannot be inferred.

```jldoctest objective_function; filter = r"Stacktrace:.*"s
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, con, x^2 == 1)
con : x² = 1.0

julia> constraint_by_name(model, "kon")

julia> constraint_by_name(model, "con")
con : x² = 1.0

julia> constraint_by_name(model, "con", AffExpr, MOI.EqualTo{Float64})

julia> constraint_by_name(model, "con", QuadExpr, MOI.EqualTo{Float64})
con : x² = 1.0
```
"""
function constraint_by_name end

function constraint_by_name(model::Model, name::String)
    index = MOI.get(backend(model), MOI.ConstraintIndex, name)
    if index isa Nothing
        return nothing
    else
        return constraint_ref_with_index(model, index)
    end
end

function constraint_by_name(
    model::Model,
    name::String,
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
)
    index = MOI.get(backend(model), MOI.ConstraintIndex{F,S}, name)
    if index isa Nothing
        return nothing
    else
        return constraint_ref_with_index(model, index)
    end
end
function constraint_by_name(
    model::Model,
    name::String,
    F::Type{<:Union{ScalarType,Vector{ScalarType}}},
    S::Type,
) where {ScalarType<:AbstractJuMPScalar}
    return constraint_by_name(model, name, moi_function_type(F), S)
end

"""
    constraint_ref_with_index(model::AbstractModel, index::MOI.ConstraintIndex)

Return a `ConstraintRef` of `model` corresponding to `index`.
"""
function constraint_ref_with_index(
    model::AbstractModel,
    index::MOI.ConstraintIndex{
        <:MOI.AbstractScalarFunction,
        <:MOI.AbstractScalarSet,
    },
)
    return ConstraintRef(model, index, ScalarShape())
end
function constraint_ref_with_index(
    model::AbstractModel,
    index::MOI.ConstraintIndex{
        <:MOI.AbstractVectorFunction,
        <:MOI.AbstractVectorSet,
    },
)
    return ConstraintRef(model, index, get(model.shapes, index, VectorShape()))
end

"""
    delete(model::Model, con_ref::ConstraintRef)

Delete the constraint associated with `constraint_ref` from the model `model`.

Note that `delete` does not unregister the name from the model, so adding a new
constraint of the same name will throw an error. Use [`unregister`](@ref) to
unregister the name after deletion as follows:
```julia
@constraint(model, c, 2x <= 1)
delete(model, c)
unregister(model, :c)
```

See also: [`unregister`](@ref)
"""
function delete(model::Model, con_ref::ConstraintRef)
    if model !== con_ref.model
        error(
            "The constraint reference you are trying to delete does not " *
            "belong to the model.",
        )
    end
    model.is_model_dirty = true
    return MOI.delete(backend(model), index(con_ref))
end

"""
    delete(model::Model, con_refs::Vector{<:ConstraintRef})

Delete the constraints associated with `con_refs` from the model `model`.
Solvers may implement specialized methods for deleting multiple constraints of
the same concrete type, i.e., when `isconcretetype(eltype(con_refs))`. These
may be more efficient than repeatedly calling the single constraint delete
method.

See also: [`unregister`](@ref)
"""
function delete(
    model::Model,
    con_refs::Vector{<:ConstraintRef{<:AbstractModel}},
)
    if any(c -> model !== c.model, con_refs)
        error("A constraint reference you are trying to delete does not" * "
            belong to the model.")
    end
    model.is_model_dirty = true
    MOI.delete(backend(model), index.(con_refs))
    return
end

"""
    is_valid(model::Model, con_ref::ConstraintRef{<:AbstractModel})

Return `true` if `constraint_ref` refers to a valid constraint in `model`.
"""
function is_valid(model::Model, con_ref::ConstraintRef{<:AbstractModel})
    return (
        model === con_ref.model && MOI.is_valid(backend(model), con_ref.index)
    )
end

#############################################################################
# AbstractConstraint
"""
    abstract type AbstractConstraint

An abstract base type for all constraint types. `AbstractConstraint`s store the
function and set directly, unlike [`ConstraintRef`](@ref)s that are merely
references to constraints stored in a model. `AbstractConstraint`s do not need
to be attached to a model.
"""
abstract type AbstractConstraint end

"""
    struct BridgeableConstraint{C, B} <: AbstractConstraint
        constraint::C
        bridge_type::B
    end

Constraint `constraint` that can be bridged by the bridge of type `bridge_type`.
Adding this constraint to a model is equivalent to
```julia
add_bridge(model, bridge_type)
add_constraint(model, constraint)
```

## Examples

Given a new scalar set type `CustomSet` with a bridge `CustomBridge` that can
bridge `F`-in-`CustomSet` constraints, when the user does
```julia
model = Model()
@variable(model, x)
@constraint(model, x + 1 in CustomSet())
optimize!(model)
```
with an optimizer that does not support `F`-in-`CustomSet` constraints, the
constraint will not be bridged unless he manually calls `add_bridge(model,
CustomBridge)`. In order to automatically add the `CustomBridge` to any model to
which an `F`-in-`CustomSet` is added, simply add the following method:
```julia
function JuMP.build_constraint(_error::Function, func::AbstractJuMPScalar,
                               set::CustomSet)
    constraint = ScalarConstraint(func, set)
    return JuMP.BridgeableConstraint(constraint, CustomBridge)
end
```

### Note

JuMP extensions should extend `JuMP.build_constraint` only if they also defined
`CustomSet`, for three
reasons:
1. It is problematic if multiple extensions overload the same JuMP method.
2. A missing method will not inform the users that they forgot to load the
   extension module defining the `build_constraint` method.
3. Defining a method where neither the function nor any of the argument types
   are defined in the package is called [*type piracy*](https://docs.julialang.org/en/v1/manual/style-guide/index.html#Avoid-type-piracy-1)
   and is discouraged in the Julia style guide.
"""
struct BridgeableConstraint{C,B} <: AbstractConstraint
    constraint::C
    bridge_type::B
end

function add_constraint(
    model::Model,
    con::BridgeableConstraint,
    name::String = "",
)
    add_bridge(model, con.bridge_type)
    return add_constraint(model, con.constraint, name)
end

"""
    jump_function(constraint::AbstractConstraint)

Return the function of the constraint `constraint` in the function-in-set form
as a `AbstractJuMPScalar` or `Vector{AbstractJuMPScalar}`.
"""
jump_function(constraint::AbstractConstraint) = constraint.func

"""
    moi_function(constraint::AbstractConstraint)

Return the function of the constraint `constraint` in the function-in-set form
as a `MathOptInterface.AbstractFunction`.
"""
function moi_function(constraint::AbstractConstraint)
    return moi_function(jump_function(constraint))
end

"""
    moi_set(constraint::AbstractConstraint)

Return the set of the constraint `constraint` in the function-in-set form as a
`MathOptInterface.AbstractSet`.

    moi_set(s::AbstractVectorSet, dim::Int)

Returns the MOI set of dimension `dim` corresponding to the JuMP set `s`.
"""
moi_set(constraint::AbstractConstraint) = constraint.set

"""
    constraint_object(con_ref::ConstraintRef)

Return the underlying constraint data for the constraint referenced by `ref`.
"""
function constraint_object end

"""
    struct ScalarConstraint

The data for a scalar constraint. The `func` field contains a JuMP object
representing the function and the `set` field contains the MOI set.
See also the [documentation](@ref Constraints) on JuMP's representation of
constraints for more background.
"""
struct ScalarConstraint{F<:AbstractJuMPScalar,S<:MOI.AbstractScalarSet} <:
       AbstractConstraint
    func::F
    set::S
end

reshape_set(set::MOI.AbstractScalarSet, ::ScalarShape) = set
shape(::ScalarConstraint) = ScalarShape()

function constraint_object(
    con_ref::ConstraintRef{
        <:AbstractModel,
        MOI.ConstraintIndex{FuncType,SetType},
    },
) where {FuncType<:MOI.AbstractScalarFunction,SetType<:MOI.AbstractScalarSet}
    model = con_ref.model
    f = MOI.get(model, MOI.ConstraintFunction(), con_ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), con_ref)::SetType
    return ScalarConstraint(jump_function(model, f), s)
end
function check_belongs_to_model(con::ScalarConstraint, model)
    return check_belongs_to_model(con.func, model)
end

"""
    struct VectorConstraint

The data for a vector constraint. The `func` field contains a JuMP object
representing the function and the `set` field contains the MOI set. The
`shape` field contains an [`AbstractShape`](@ref) matching the form in which
the constraint was constructed (e.g., by using matrices or flat vectors).
See also the [documentation](@ref Constraints) on JuMP's representation of
constraints.
"""
struct VectorConstraint{
    F<:AbstractJuMPScalar,
    S<:MOI.AbstractVectorSet,
    Shape<:AbstractShape,
} <: AbstractConstraint
    func::Vector{F}
    set::S
    shape::Shape
end
function VectorConstraint(
    func::Vector{<:AbstractJuMPScalar},
    set::MOI.AbstractVectorSet,
)
    if length(func) != MOI.dimension(set)
        throw(
            DimensionMismatch(
                "Dimension of the function $(length(func)) does not match the " *
                "dimension of the set $(set).",
            ),
        )
    end
    return VectorConstraint(func, set, VectorShape())
end

function VectorConstraint(
    func::AbstractVector{<:AbstractJuMPScalar},
    set::MOI.AbstractVectorSet,
)
    # collect() is not used here so that DenseAxisArray will work
    f = [func[idx] for idx in eachindex(func)]
    return VectorConstraint(f, set)
end

reshape_set(set::MOI.AbstractVectorSet, ::VectorShape) = set
shape(con::VectorConstraint) = con.shape
function constraint_object(
    con_ref::ConstraintRef{
        <:AbstractModel,
        MOI.ConstraintIndex{FuncType,SetType},
    },
) where {FuncType<:MOI.AbstractVectorFunction,SetType<:MOI.AbstractVectorSet}
    model = con_ref.model
    f = MOI.get(model, MOI.ConstraintFunction(), con_ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), con_ref)::SetType
    return VectorConstraint(jump_function(model, f), s, con_ref.shape)
end
function check_belongs_to_model(con::VectorConstraint, model)
    for func in con.func
        check_belongs_to_model(func, model)
    end
end

function _moi_add_constraint(
    model::MOI.ModelLike,
    f::F,
    s::S,
) where {F<:MOI.AbstractFunction,S<:MOI.AbstractSet}
    if !MOI.supports_constraint(model, F, S)
        error(
            "Constraints of type $(F)-in-$(S) are not supported by the " *
            "solver.\n\nIf you expected the solver to support your problem, " *
            "you may have an error in your formulation. Otherwise, consider " *
            "using a different solver.\n\nThe list of available solvers, " *
            "along with the problem types they support, is available at " *
            "https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers.",
        )
    end
    return MOI.add_constraint(model, f, s)
end

"""
    add_constraint(model::Model, con::AbstractConstraint, name::String="")

Add a constraint `con` to `Model model` and sets its name.
"""
function add_constraint(
    model::Model,
    con::AbstractConstraint,
    name::String = "",
)
    # The type of backend(model) is unknown so we directly redirect to another
    # function.
    check_belongs_to_model(con, model)
    func, set = moi_function(con), moi_set(con)
    cindex = _moi_add_constraint(
        backend(model),
        func,
        set,
    )::MOI.ConstraintIndex{typeof(func),typeof(set)}
    cshape = shape(con)
    if !(cshape isa ScalarShape) && !(cshape isa VectorShape)
        model.shapes[cindex] = cshape
    end
    con_ref = ConstraintRef(model, cindex, cshape)
    # Only set names if appropriate!
    if !(func isa MOI.VariableIndex) &&
       !isempty(name) &&
       MOI.supports(backend(model), MOI.ConstraintName(), typeof(cindex))
        set_name(con_ref, name)
    end
    model.is_model_dirty = true
    return con_ref
end

"""
    set_normalized_coefficient(con_ref::ConstraintRef, variable::VariableRef, value)

Set the coefficient of `variable` in the constraint `constraint` to `value`.

Note that prior to this step, JuMP will aggregate multiple terms containing the
same variable. For example, given a constraint `2x + 3x <= 2`,
`set_normalized_coefficient(con, x, 4)` will create the constraint `4x <= 2`.

```jldoctest; filter=r"≤|<="
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, con, 2x + 3x <= 2)
con : 5 x ≤ 2.0

julia> set_normalized_coefficient(con, x, 4)

julia> con
con : 4 x ≤ 2.0
```
"""
function set_normalized_coefficient(
    con_ref::ConstraintRef{<:AbstractModel,MOI.ConstraintIndex{F,S}},
    variable,
    value,
) where {
    S,
    T,
    F<:Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}},
}
    model = owner_model(con_ref)
    MOI.modify(
        backend(model),
        index(con_ref),
        MOI.ScalarCoefficientChange(index(variable), convert(T, value)),
    )
    model.is_model_dirty = true
    return
end

"""
    set_normalized_coefficients(
        con_ref::ConstraintRef,
        variable,
        new_coefficients::Vector{Tuple{Int64,T}},
    )

Set the coefficients of `variable` in the constraint `con_ref` to
`new_coefficients`, where each element in `new_coefficients` is a tuple which
maps the row to a new coefficient.

Note that prior to this step, during constraint creation, JuMP will aggregate
multiple terms containing the same variable.

```jldoctest; filter=r"≤|<="
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, con, [2x + 3x, 4x] in MOI.Nonnegatives(2))
con : [5 x, 4 x] ∈ MathOptInterface.Nonnegatives(2)

julia> set_normalized_coefficients(con, x, [(1, 2.0), (2, 5.0)])

julia> con
con : [2 x, 5 x] ∈ MathOptInterface.Nonnegatives(2)
```
"""
function set_normalized_coefficients(
    con_ref::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex{F}},
    variable,
    new_coefficients::Vector{Tuple{Int64,T}},
) where {T,F<:Union{MOI.VectorAffineFunction{T},MOI.VectorQuadraticFunction{T}}}
    model = owner_model(con_ref)
    MOI.modify(
        backend(model),
        index(con_ref),
        MOI.MultirowChange(index(variable), new_coefficients),
    )
    model.is_model_dirty = true
    return
end

"""
    normalized_coefficient(con_ref::ConstraintRef, variable::VariableRef)

Return the coefficient associated with `variable` in `constraint` after JuMP has
normalized the constraint into its standard form. See also
[`set_normalized_coefficient`](@ref).
"""
function normalized_coefficient(
    con_ref::ConstraintRef{<:AbstractModel,MOI.ConstraintIndex{F,S}},
    variable,
) where {
    S,
    T,
    F<:Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}},
}
    con = constraint_object(con_ref)
    return _affine_coefficient(con.func, variable)
end

"""
    set_normalized_rhs(con_ref::ConstraintRef, value)

Set the right-hand side term of `constraint` to `value`.

Note that prior to this step, JuMP will aggregate all constant terms onto the
right-hand side of the constraint. For example, given a constraint `2x + 1 <=
2`, `set_normalized_rhs(con, 4)` will create the constraint `2x <= 4`, not `2x +
1 <= 4`.

```jldoctest; filter=r"≤|<="
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, con, 2x + 1 <= 2)
con : 2 x <= 1.0

julia> set_normalized_rhs(con, 4)

julia> con
con : 2 x <= 4.0
```
"""
function set_normalized_rhs(
    con_ref::ConstraintRef{<:AbstractModel,MOI.ConstraintIndex{F,S}},
    value,
) where {
    T,
    S<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T}},
    F<:Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}},
}
    MOI.set(
        owner_model(con_ref),
        MOI.ConstraintSet(),
        con_ref,
        S(convert(T, value)),
    )
    return
end

"""
    normalized_rhs(con_ref::ConstraintRef)

Return the right-hand side term of `con_ref` after JuMP has converted the
constraint into its normalized form. See also [`set_normalized_rhs`](@ref).
"""
function normalized_rhs(
    con_ref::ConstraintRef{<:AbstractModel,MOI.ConstraintIndex{F,S}},
) where {
    T,
    S<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T}},
    F<:Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}},
}
    con = constraint_object(con_ref)
    return MOI.constant(con.set)
end

function _moi_add_to_function_constant(
    model::MOI.ModelLike,
    ci::MOI.ConstraintIndex{
        <:MOI.AbstractScalarFunction,
        <:MOI.AbstractScalarSet,
    },
    value,
)
    set = MOI.get(model, MOI.ConstraintSet(), ci)
    if !MOI.Utilities.supports_shift_constant(typeof(set))
        error(
            "Unable to add to function constant for constraint type " *
            "$(typeof(ci))",
        )
    end
    new_set = MOIU.shift_constant(set, convert(Float64, -value))
    return MOI.set(model, MOI.ConstraintSet(), ci, new_set)
end
function _moi_add_to_function_constant(
    model::MOI.ModelLike,
    ci::MOI.ConstraintIndex{
        <:Union{MOI.VectorAffineFunction,MOI.VectorQuadraticFunction},
        <:MOI.AbstractVectorSet,
    },
    value,
)
    func = MOI.get(model, MOI.ConstraintFunction(), ci)
    new_constant = value + MOI.constant(func)
    return MOI.modify(model, ci, MOI.VectorConstantChange(new_constant))
end

"""
    add_to_function_constant(constraint::ConstraintRef, value)

Add `value` to the function constant term.

Note that for scalar constraints, JuMP will aggregate all constant terms onto the
right-hand side of the constraint so instead of modifying the function, the set
will be translated by `-value`. For example, given a constraint `2x <=
3`, `add_to_function_constant(c, 4)` will modify it to `2x <= -1`.

## Examples

For scalar constraints, the set is translated by `-value`:
```jldoctest; filter=r"≤|<="
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, con, 0 <= 2x - 1 <= 2)
con : 2 x ∈ [1.0, 3.0]

julia> add_to_function_constant(con, 4)

julia> con
con : 2 x ∈ [-3.0, -1.0]
```

For vector constraints, the constant is added to the function:
```jldoctest; filter=r"≤|<="
julia> model = Model();

julia> @variable(model, x);

julia> @variable(model, y);

julia> @constraint(model, con, [x + y, x, y] in SecondOrderCone())
con : [x + y, x, y] ∈ MathOptInterface.SecondOrderCone(3)

julia> add_to_function_constant(con, [1, 2, 2])

julia> con
con : [x + y + 1, x + 2, y + 2] ∈ MathOptInterface.SecondOrderCone(3)
```

"""
function add_to_function_constant(
    constraint::ConstraintRef{<:AbstractModel},
    value,
)
    model = owner_model(constraint)
    # The type of `backend(model)` is not type-stable, so we use a function
    # barrier (`_moi_add_to_function_constant`) to improve performance.
    _moi_add_to_function_constant(backend(model), index(constraint), value)
    model.is_model_dirty = true
    return
end

"""
    value(con_ref::ConstraintRef; result::Int = 1)

Return the primal value of constraint `con_ref` associated with result index
`result` of the most-recent solution returned by the solver.

That is, if `con_ref` is the reference of a constraint `func`-in-`set`, it
returns the value of `func` evaluated at the value of the variables (given by
[`value(::VariableRef)`](@ref)).

Use [`has_values`](@ref) to check if a result exists before asking for values.

See also: [`result_count`](@ref).

## Note

For scalar constraints, the constant is moved to the `set` so it is not taken
into account in the primal value of the constraint. For instance, the constraint
`@constraint(model, 2x + 3y + 1 == 5)` is transformed into
`2x + 3y`-in-`MOI.EqualTo(4)` so the value returned by this function is the
evaluation of `2x + 3y`.
```
"""
function value(
    con_ref::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex};
    result::Int = 1,
)
    return reshape_vector(_constraint_primal(con_ref, result), con_ref.shape)
end

"""
    value(var_value::Function, con_ref::ConstraintRef)

Evaluate the primal value of the constraint `con_ref` using `var_value(v)`
as the value for each variable `v`.
"""
function value(
    var_value::Function,
    con_ref::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex},
)
    f = jump_function(constraint_object(con_ref))
    return reshape_vector(value.(var_value, f), con_ref.shape)
end

# Returns the value of MOI.ConstraintPrimal in a type-stable way
function _constraint_primal(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractScalarFunction,
            <:MOI.AbstractScalarSet,
        },
    },
    result::Int,
)::Float64
    return MOI.get(con_ref.model, MOI.ConstraintPrimal(result), con_ref)
end
function _constraint_primal(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractVectorFunction,
            <:MOI.AbstractVectorSet,
        },
    },
    result,
)::Vector{Float64}
    return MOI.get(con_ref.model, MOI.ConstraintPrimal(result), con_ref)
end

"""
    has_duals(model::Model; result::Int = 1)

Return `true` if the solver has a dual solution in result index `result`
available to query, otherwise return `false`.

See also [`dual`](@ref), [`shadow_price`](@ref), and [`result_count`](@ref).
"""
function has_duals(model::Model; result::Int = 1)
    return dual_status(model; result = result) != MOI.NO_SOLUTION
end

"""
    dual(con_ref::ConstraintRef; result::Int = 1)

Return the dual value of constraint `con_ref` associated with result index
`result` of the most-recent solution returned by the solver.

Use `has_dual` to check if a result exists before asking for values.

See also: [`result_count`](@ref), [`shadow_price`](@ref).
"""
function dual(
    con_ref::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex};
    result::Int = 1,
)
    return reshape_vector(
        _constraint_dual(con_ref, result),
        dual_shape(con_ref.shape),
    )
end

# Returns the value of MOI.ConstraintPrimal in a type-stable way
function _constraint_dual(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractScalarFunction,
            <:MOI.AbstractScalarSet,
        },
    },
    result::Int,
)::Float64
    return MOI.get(con_ref.model, MOI.ConstraintDual(result), con_ref)
end
function _constraint_dual(
    con_ref::ConstraintRef{
        <:AbstractModel,
        <:MOI.ConstraintIndex{
            <:MOI.AbstractVectorFunction,
            <:MOI.AbstractVectorSet,
        },
    },
    result::Int,
)::Vector{Float64}
    return MOI.get(con_ref.model, MOI.ConstraintDual(result), con_ref)
end

"""
    shadow_price(con_ref::ConstraintRef)

Return the change in the objective from an infinitesimal relaxation of the
constraint.

This value is computed from [`dual`](@ref) and can be queried only when
`has_duals` is `true` and the objective sense is `MIN_SENSE` or `MAX_SENSE`
(not `FEASIBILITY_SENSE`). For linear constraints, the shadow prices differ at
most in sign from the `dual` value depending on the objective sense.

See also [`reduced_cost`](@ref JuMP.reduced_cost).

## Notes

- The function simply translates signs from `dual` and does not validate
  the conditions needed to guarantee the sensitivity interpretation of the
  shadow price. The caller is responsible, e.g., for checking whether the solver
  converged to an optimal primal-dual pair or a proof of infeasibility.
- The computation is based on the current objective sense of the model. If this
  has changed since the last solve, the results will be incorrect.
- Relaxation of equality constraints (and hence the shadow price) is defined
  based on which sense of the equality constraint is active.
"""
function shadow_price(
    con_ref::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex},
)
    return error(
        "The shadow price is not defined or not implemented for this type " *
        "of constraint.",
    )
end

function _shadow_price_less_than(dual_value, sense::MOI.OptimizationSense)
    # When minimizing, the shadow price is nonpositive and when maximizing the
    # shadow price is nonnegative (because relaxing a constraint can only
    # improve the objective). By MOI convention, a feasible dual on a LessThan
    # set is nonpositive, so we flip the sign when maximizing.
    if sense == MAX_SENSE
        return -dual_value
    elseif sense == MIN_SENSE
        return dual_value
    else
        error(
            "The shadow price is not available because the objective sense " *
            "$sense is not minimization or maximization.",
        )
    end
end

function _shadow_price_greater_than(dual_value, sense::MOI.OptimizationSense)
    # By MOI convention, a feasible dual on a GreaterThan set is nonnegative,
    # so we flip the sign when minimizing. (See comment in the method above).
    if sense == MAX_SENSE
        return dual_value
    elseif sense == MIN_SENSE
        return -dual_value
    else
        error(
            "The shadow price is not available because the objective sense " *
            "$sense is not minimization or maximization.",
        )
    end
end

function shadow_price(
    con_ref::ConstraintRef{<:AbstractModel,MOI.ConstraintIndex{F,S}},
) where {S<:MOI.LessThan,F}
    model = con_ref.model
    if !has_duals(model)
        error(
            "The shadow price is not available because no dual result is " *
            "available.",
        )
    end
    return _shadow_price_less_than(dual(con_ref), objective_sense(model))
end

function shadow_price(
    con_ref::ConstraintRef{<:AbstractModel,MOI.ConstraintIndex{F,S}},
) where {S<:MOI.GreaterThan,F}
    model = con_ref.model
    if !has_duals(model)
        error(
            "The shadow price is not available because no dual result is " *
            "available.",
        )
    end
    return _shadow_price_greater_than(dual(con_ref), objective_sense(model))
end

function shadow_price(
    con_ref::ConstraintRef{<:AbstractModel,MOI.ConstraintIndex{F,S}},
) where {S<:MOI.EqualTo,F}
    model = con_ref.model
    if !has_duals(model)
        error(
            "The shadow price is not available because no dual result is " *
            "available.",
        )
    end
    sense = objective_sense(model)
    dual_val = dual(con_ref)
    if dual_val > 0
        # Treat the equality constraint as if it were a GreaterThan constraint.
        return _shadow_price_greater_than(dual_val, sense)
    else
        # Treat the equality constraint as if it were a LessThan constraint.
        return _shadow_price_less_than(dual_val, sense)
    end
end

function _error_if_not_concrete_type(t)
    if !isconcretetype(t)
        error("`$t` is not a concrete type. Did you miss a type parameter?")
    end
    return
end
# `isconcretetype(Vector{Integer})` is `true`
function _error_if_not_concrete_type(t::Type{Vector{ElT}}) where {ElT}
    return _error_if_not_concrete_type(ElT)
end

"""
    num_constraints(model::Model, function_type, set_type)::Int64

Return the number of constraints currently in the model where the function
has type `function_type` and the set has type `set_type`.

See also [`list_of_constraint_types`](@ref) and [`all_constraints`](@ref).

# Example
```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0, Bin);

julia> @variable(model, y);

julia> @constraint(model, y in MOI.GreaterThan(1.0));

julia> @constraint(model, y <= 1.0);

julia> @constraint(model, 2x <= 1);

julia> num_constraints(model, VariableRef, MOI.GreaterThan{Float64})
2

julia> num_constraints(model, VariableRef, MOI.ZeroOne)
1

julia> num_constraints(model, AffExpr, MOI.LessThan{Float64})
2
```
"""
function num_constraints(
    model::Model,
    function_type::Type{
        <:Union{AbstractJuMPScalar,Vector{<:AbstractJuMPScalar}},
    },
    set_type::Type{<:MOI.AbstractSet},
)::Int64
    _error_if_not_concrete_type(function_type)
    _error_if_not_concrete_type(set_type)
    # TODO: Support JuMP's set helpers like SecondOrderCone().
    f_type = moi_function_type(function_type)
    return MOI.get(model, MOI.NumberOfConstraints{f_type,set_type}())
end

"""
    all_constraints(model::Model, function_type, set_type)::Vector{<:ConstraintRef}

Return a list of all constraints currently in the model where the function
has type `function_type` and the set has type `set_type`. The constraints are
ordered by creation time.

See also [`list_of_constraint_types`](@ref) and [`num_constraints`](@ref).

# Example
```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0, Bin);

julia> @constraint(model, 2x <= 1);

julia> all_constraints(model, VariableRef, MOI.GreaterThan{Float64})
1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.VariableIndex,MathOptInterface.GreaterThan{Float64}},ScalarShape},1}:
 x ≥ 0.0

julia> all_constraints(model, VariableRef, MOI.ZeroOne)
1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.VariableIndex,MathOptInterface.ZeroOne},ScalarShape},1}:
 x binary

julia> all_constraints(model, AffExpr, MOI.LessThan{Float64})
1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},1}:
 2 x ≤ 1.0
```
"""
function all_constraints(
    model::Model,
    function_type::Type{
        <:Union{AbstractJuMPScalar,Vector{<:AbstractJuMPScalar}},
    },
    set_type::Type{<:MOI.AbstractSet},
)
    _error_if_not_concrete_type(function_type)
    _error_if_not_concrete_type(set_type)
    # TODO: Support JuMP's set helpers like SecondOrderCone().
    f_type = moi_function_type(function_type)
    if set_type <: MOI.AbstractScalarSet
        constraint_ref_type = ConstraintRef{
            Model,
            MOI.ConstraintIndex{f_type,set_type},
            ScalarShape,
        }
    else
        constraint_ref_type =
            ConstraintRef{Model,MOI.ConstraintIndex{f_type,set_type}}
    end
    result = constraint_ref_type[]
    for idx in MOI.get(model, MOI.ListOfConstraintIndices{f_type,set_type}())
        push!(result, constraint_ref_with_index(model, idx))
    end
    return result
end

# TODO: Support vector function types. This is blocked by not having the shape
# information available.

"""
    list_of_constraint_types(model::Model)::Vector{Tuple{Type,Type}}

Return a list of tuples of the form `(F, S)` where `F` is a JuMP function type
and `S` is an MOI set type such that `all_constraints(model, F, S)` returns
a nonempty list.

# Example
```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0, Bin);

julia> @constraint(model, 2x <= 1);

julia> list_of_constraint_types(model)
3-element Array{Tuple{Type,Type},1}:
 (GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64})
 (VariableRef, MathOptInterface.GreaterThan{Float64})
 (VariableRef, MathOptInterface.ZeroOne)
```

## Performance considerations

Iterating over the list of function and set types is a type-unstable operation.
Consider using a function barrier. See the [Performance tips for extensions](@ref)
section of the documentation for more details.
"""
function list_of_constraint_types(model::Model)::Vector{Tuple{Type,Type}}
    # We include an annotated return type here because Julia fails terribly at
    # inferring it, even though we annotate the type of the return vector.
    return Tuple{Type,Type}[
        (jump_function_type(model, F), S) for
        (F, S) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
    ]
end

"""
    num_constraints(model::Model; count_variable_in_set_constraints::Bool)

Return the number of constraints in `model`.

If `count_variable_in_set_constraints == true`, then `VariableRef` constraints
such as `VariableRef`-in-`Integer` are included. To count only the number of
structural constraints (e.g., the rows in the constraint matrix of a linear
program), pass `count_variable_in_set_constraints = false`.

## Examples

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0, Int);

julia> @constraint(model, 2x <= 1);

julia> num_constraints(model; count_variable_in_set_constraints = true)
3

julia> num_constraints(model; count_variable_in_set_constraints = false)
1
```
"""
function num_constraints(model::Model; count_variable_in_set_constraints::Bool)
    ret = num_nonlinear_constraints(model)
    for (F, S) in list_of_constraint_types(model)
        if F != VariableRef || count_variable_in_set_constraints
            ret += num_constraints(model, F, S)
        end
    end
    return ret
end

"""
    all_constraints(
        model::Model;
        include_variable_in_set_constraints::Bool,
    )::Vector{ConstraintRef}

Return a list of all constraints in `model`.

If `include_variable_in_set_constraints == true`, then `VariableRef` constraints
such as `VariableRef`-in-`Integer` are included. To return only the structural
constraints (e.g., the rows in the constraint matrix of a linear program), pass
`include_variable_in_set_constraints = false`.

## Examples

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0, Int);

julia> @constraint(model, 2x <= 1);

julia> @NLconstraint(model, x^2 <= 1);

julia> all_constraints(model; include_variable_in_set_constraints = true)
4-element Vector{ConstraintRef}:
 2 x ≤ 1.0
 x ≥ 0.0
 x integer
 x ^ 2.0 - 1.0 ≤ 0

julia> all_constraints(model; include_variable_in_set_constraints = false)
2-element Vector{ConstraintRef}:
 2 x ≤ 1.0
 x ^ 2.0 - 1.0 ≤ 0
```

## Performance considerations

Note that this function is type-unstable because it returns an abstractly typed
vector. If performance is a problem, consider using [`list_of_constraint_types`](@ref)
and a function barrier. See the [Performance tips for extensions](@ref) section
of the documentation for more details.
"""
function all_constraints(
    model::Model;
    include_variable_in_set_constraints::Bool,
)
    ret = ConstraintRef[]
    for (F, S) in list_of_constraint_types(model)
        if F != VariableRef || include_variable_in_set_constraints
            append!(ret, all_constraints(model, F, S))
        end
    end
    append!(ret, all_nonlinear_constraints(model))
    return ret
end

"""
    relax_with_penalty!(
        model::Model,
        [penalties::Dict{ConstraintRef,Float64}];
        [default::Union{Nothing,Real} = nothing,]
    )

Destructively modify the model in-place to create a penalized relaxation of the
constraints.

!!! warning
    This is a destructive routine that modifies the model in-place. If you don't
    want to modify the original model, use [`copy_model`](@ref) to create a copy
    before calling [`relax_with_penalty!`](@ref).

## Reformulation

See [`MOI.Utilities.ScalarPenaltyRelaxation`](@ref) for details of the
reformulation.

For each constraint `ci`, the penalty passed to
[`MOI.Utilities.ScalarPenaltyRelaxation`](@ref) is `get(penalties, ci, default)`.
If the value is `nothing`, because `ci` does not exist in `penalties` and
`default = nothing`, then the constraint is skipped.

## Return value

This function returns a `Dict{ConstraintRef,AffExpr}` that maps each constraint
index to the corresponding `y + z` as an [`AffExpr`](@ref). In an optimal
solution, query the value of these functions to compute the violation of each
constraint.

## Relax a subset of constraints

To relax a subset of constraints, pass a `penalties` dictionary and set
`default = nothing`.

## Examples

```jldoctest
julia> function new_model()
           model = Model()
           @variable(model, x)
           @objective(model, Max, 2x + 1)
           @constraint(model, c1, 2x - 1 <= -2)
           @constraint(model, c2, 3x >= 0)
           return model
       end
new_model (generic function with 1 method)

julia> model_1 = new_model();

julia> relax_with_penalty!(model_1; default = 2.0)
Dict{ConstraintRef{Model, C, ScalarShape} where C, AffExpr} with 2 entries:
  c1 : 2 x - _[3] ≤ -1.0 => _[3]
  c2 : 3 x + _[2] ≥ 0.0  => _[2]

julia> print(model_1)
Max 2 x - 2 _[2] - 2 _[3] + 1
Subject to
 c2 : 3 x + _[2] ≥ 0.0
 c1 : 2 x - _[3] ≤ -1.0
 _[2] ≥ 0.0
 _[3] ≥ 0.0

julia> model_2 = new_model();

julia> relax_with_penalty!(model_2, Dict(model_2[:c2] => 3.0))
Dict{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.GreaterThan{Float64}}, ScalarShape}, AffExpr} with 1 entry:
  c2 : 3 x + _[2] ≥ 0.0 => _[2]

julia> print(model_2)
Max 2 x - 3 _[2] + 1
Subject to
 c2 : 3 x + _[2] ≥ 0.0
 c1 : 2 x ≤ -1.0
 _[2] ≥ 0.0
```
"""
function relax_with_penalty!(
    model::Model,
    penalties::Dict;
    default::Union{Nothing,Real} = nothing,
)
    if default !== nothing
        default = Float64(default)
    end
    moi_penalties = Dict{MOI.ConstraintIndex,Float64}(
        index(k) => Float64(v) for (k, v) in penalties
    )
    map = MOI.modify(
        backend(model),
        MOI.Utilities.PenaltyRelaxation(moi_penalties; default = default),
    )
    return Dict(
        constraint_ref_with_index(model, k) => jump_function(model, v) for
        (k, v) in map
    )
end

function relax_with_penalty!(model::Model; default::Real = 1.0)
    return relax_with_penalty!(model, Dict(); default = default)
end
