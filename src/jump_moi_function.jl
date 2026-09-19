#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    moi_function_type(::Type{T}) where {T}

Given a JuMP object type `T`, return the MathOptInterface equivalent.

See also: [`jump_function_type`](@ref).

## Example

```jldoctest
julia> moi_function_type(AffExpr)
MathOptInterface.ScalarAffineFunction{Float64}
```
"""
function moi_function_type end

"""
    moi_function(model::GenericModel, x::AbstractJuMPScalar)
    moi_function(model::GenericModel, x::AbstractArray{<:AbstractJuMPScalar})

Check that `x` belongs to `model`, then return the MathOptInterface equivalent
of the function `x`.

This is equivalent to calling `check_belongs_to_model(x, model)` followed by
`moi_function(x)`, but some methods may fuse the ownership check and conversion
into a single function call, and some types do not support the single-argument
    [`moi_function`](@ref).

See also: [`jump_function`](@ref).

!!! compat
    The `model` argument was added in JuMP v1.31.  New functions should use the
    two-argument version.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> f = 2.0 * x + 1.0
2 x + 1

julia> moi_function(model, f)
1.0 + 2.0 MOI.VariableIndex(1)
```
"""
function moi_function end

# Default fallback for backwards compatibility. The first argument `model` was
# introduced in JuMP@1.31.0.
function moi_function(model::GenericModel, f)
    check_belongs_to_model(f, model)
    return moi_function(f)
end

# Plasmo combines variables from multiple models, so
# check_belongs_to_model(owner_model(f), f) may not work.
moi_function(model, f) = moi_function(f)

"""
    check_belongs_to_model(x::AbstractJuMPScalar, model::AbstractModel)

Throw [`VariableNotOwned`](@ref) if the [`owner_model`](@ref) of `x` is not
`model`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> check_belongs_to_model(x, model)

julia> model_2 = Model();

julia> check_belongs_to_model(x, model_2)
ERROR: VariableNotOwned{VariableRef}(x): the variable x cannot be used in this model because
it belongs to a different model.
[...]
```
"""
function check_belongs_to_model end

"""
    jump_function_type(model::AbstractModel, ::Type{T}) where {T}

Given an MathOptInterface object type `T`, return the JuMP equivalent.

See also: [`moi_function_type`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> jump_function_type(model, MOI.ScalarAffineFunction{Float64})
AffExpr (alias for GenericAffExpr{Float64, GenericVariableRef{Float64}})
```
"""
function jump_function_type end

"""
    jump_function(model::AbstractModel, x::MOI.AbstractFunction)

Given an MathOptInterface object `x`, return the JuMP equivalent.

See also: [`moi_function`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> f = 2.0 * index(x) + 1.0
1.0 + 2.0 MOI.VariableIndex(1)

julia> jump_function(model, f)
2 x + 1
```
"""
function jump_function end

# MOI.VariableIndex

moi_function_type(::Type{<:AbstractVariableRef}) = MOI.VariableIndex

moi_function(variable::AbstractVariableRef) = index(variable)

function check_belongs_to_model(v::AbstractVariableRef, model::AbstractModel)
    if owner_model(v) !== model
        throw(VariableNotOwned(v))
    end
    return
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VariableIndex},
) where {T}
    return GenericVariableRef{T}
end

function jump_function(
    model::GenericModel{T},
    variable::MOI.VariableIndex,
) where {T}
    return GenericVariableRef{T}(model, variable)
end

# MOI.ScalarAffineFunction

function moi_function_type(::Type{<:GenericAffExpr{T}}) where {T}
    return MOI.ScalarAffineFunction{T}
end

moi_function(a::GenericAffExpr) = MOI.ScalarAffineFunction(a)

function check_belongs_to_model(a::GenericAffExpr, model::AbstractModel)
    for variable in keys(a.terms)
        check_belongs_to_model(variable, model)
    end
    return
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.ScalarAffineFunction{C}},
) where {C,T}
    S = promote_type(C, T)
    return GenericAffExpr{S,GenericVariableRef{T}}
end

function jump_function(
    model::GenericModel{T},
    f::MOI.ScalarAffineFunction{C},
) where {C,T}
    S = promote_type(C, T)
    return GenericAffExpr{S,GenericVariableRef{T}}(model, f)
end

# MOI.ScalarQuadraticFunction

function moi_function_type(::Type{<:GenericQuadExpr{T}}) where {T}
    return MOI.ScalarQuadraticFunction{T}
end

function moi_function(aff::GenericQuadExpr)
    return MOI.ScalarQuadraticFunction(aff)
end

function check_belongs_to_model(q::GenericQuadExpr, model::AbstractModel)
    check_belongs_to_model(q.aff, model)
    for variable_pair in keys(q.terms)
        check_belongs_to_model(variable_pair.a, model)
        check_belongs_to_model(variable_pair.b, model)
    end
    return
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.ScalarQuadraticFunction{C}},
) where {C,T}
    return GenericQuadExpr{promote_type(T, C),GenericVariableRef{T}}
end

function jump_function(
    model::GenericModel{T},
    f::MOI.ScalarQuadraticFunction{C},
) where {C,T}
    S = promote_type(T, C)
    return GenericQuadExpr{S,GenericVariableRef{T}}(model, f)
end

# MOI.ScalarNonlinearFunction

moi_function_type(::Type{<:GenericNonlinearExpr}) = MOI.ScalarNonlinearFunction

function moi_function(model::GenericModel, f::GenericNonlinearExpr{V}) where {V}
    if (cache = get(model.subexpressions, f, nothing)) !== nothing
        return cache
    end
    ret = MOI.ScalarNonlinearFunction(f.head, similar(f.args))
    stack = Tuple{MOI.ScalarNonlinearFunction,Int,GenericNonlinearExpr{V}}[]
    for i in length(f.args):-1:1
        if f.args[i] isa GenericNonlinearExpr{V}
            push!(stack, (ret, i, f.args[i]))
        else
            ret.args[i] = moi_function(model, f.args[i])
        end
    end
    while !isempty(stack)
        parent, i, arg = pop!(stack)
        if (cache = get(model.subexpressions, arg, nothing)) !== nothing
            parent.args[i] = cache
            continue
        end
        child = MOI.ScalarNonlinearFunction(arg.head, similar(arg.args))
        parent.args[i] = child
        for j in length(arg.args):-1:1
            if arg.args[j] isa GenericNonlinearExpr{V}
                push!(stack, (child, j, arg.args[j]))
            else
                child.args[j] = moi_function(model, arg.args[j])
            end
        end
        model.subexpressions[arg] = child
    end
    model.subexpressions[f] = ret
    return ret
end

# A backwards-compatible function to preserve behavior prior to #4032. As one
# example, this method was used by Plasmo.jl.
function moi_function(f::GenericNonlinearExpr{V}) where {V}
    model = owner_model(f)
    if model isa GenericModel
        # If `f` has a `GenericModel` as its owner, redirect to the two-arg
        # vesion so that we can cache common subexpressions.
        return moi_function(model, f)
    end
    # There are two reasons we might reach here:
    #  1. The function `f` has no `AbstractJuMPScalar` terms, like
    #     `NonlinearExpr(:+, Any[0.0])`. In this case, `model === nothing`, and
    #     we use the single-argument method below.
    #  2. The function `f` contains terms from a JuMP extension, for example
    #     InfiniteOpt. We call the two-argument version in case they have
    #     implemented it.
    ret = MOI.ScalarNonlinearFunction(f.head, similar(f.args))
    stack = Tuple{MOI.ScalarNonlinearFunction,Int,GenericNonlinearExpr{V}}[]
    for i in length(f.args):-1:1
        if f.args[i] isa GenericNonlinearExpr{V}
            push!(stack, (ret, i, f.args[i]))
        elseif model === nothing
            ret.args[i] = moi_function(f.args[i])
        else
            ret.args[i] = moi_function(model, f.args[i])
        end
    end
    while !isempty(stack)
        parent, i, arg = pop!(stack)
        child = MOI.ScalarNonlinearFunction(arg.head, similar(arg.args))
        parent.args[i] = child
        for j in length(arg.args):-1:1
            if arg.args[j] isa GenericNonlinearExpr{V}
                push!(stack, (child, j, arg.args[j]))
            elseif model === nothing
                child.args[j] = moi_function(arg.args[j])
            else
                child.args[j] = moi_function(model, arg.args[j])
            end
        end
    end
    return ret
end

function check_belongs_to_model(
    expr::GenericNonlinearExpr,
    model::AbstractModel,
)
    # TODO: Consider keeping an `IdDict` of visited expressions so that aliases
    # are checked only once. This traversal treats the expression as a tree, so
    # repeatedly aliased subexpressions can cause the work to grow
    # exponentially in the depth of the expression, even though the underlying
    # expression is a much smaller DAG. This is not urgent because JuMP's
    # internal conversion path checks ownership while converting and caches
    # aliases; this method is now used only when a user calls it directly.
    stack = Any[expr]
    while !isempty(stack)
        child = pop!(stack)
        if child isa GenericNonlinearExpr
            for arg in child.args
                push!(stack, arg)
            end
        elseif child isa AbstractJuMPScalar
            check_belongs_to_model(child, model)
        end
    end
    return
end

function jump_function_type(
    model::GenericModel,
    ::Type{<:MOI.ScalarNonlinearFunction},
)
    return GenericNonlinearExpr{variable_ref_type(typeof(model))}
end

function jump_function(model::GenericModel, f::MOI.ScalarNonlinearFunction)
    V = variable_ref_type(typeof(model))
    ret = GenericNonlinearExpr{V}(f.head, Any[])
    stack = Tuple{GenericNonlinearExpr,Any}[]
    for arg in reverse(f.args)
        push!(stack, (ret, arg))
    end
    while !isempty(stack)
        parent, arg = pop!(stack)
        if arg isa MOI.ScalarNonlinearFunction
            new_ret = GenericNonlinearExpr{V}(arg.head, Any[])
            push!(parent.args, new_ret)
            for child in reverse(arg.args)
                push!(stack, (new_ret, child))
            end
        else
            push!(parent.args, jump_function(model, arg))
        end
    end
    return ret
end

# MOI.VectorOfVariables

function moi_function_type(::Type{<:Vector{<:AbstractVariableRef}})
    return MOI.VectorOfVariables
end

function moi_function(variables::Vector{<:AbstractVariableRef})
    return MOI.VectorOfVariables(variables)
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VectorOfVariables},
) where {T}
    return Vector{GenericVariableRef{T}}
end

function jump_function(
    model::GenericModel{T},
    variables::MOI.VectorOfVariables,
) where {T}
    return GenericVariableRef{T}[
        GenericVariableRef{T}(model, v) for v in variables.variables
    ]
end

# MOI.VectorAffineFunction

function moi_function_type(::Type{<:Vector{<:GenericAffExpr{T}}}) where {T}
    return MOI.VectorAffineFunction{T}
end

moi_function(a::Vector{<:GenericAffExpr}) = MOI.VectorAffineFunction(a)

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VectorAffineFunction{C}},
) where {C,T}
    S = promote_type(C, T)
    return Vector{GenericAffExpr{S,GenericVariableRef{T}}}
end

function jump_function(
    model::GenericModel{T},
    f::MOI.VectorAffineFunction{C},
) where {T,C}
    S = promote_type(C, T)
    ret = GenericAffExpr{S,GenericVariableRef{T}}[]
    for scalar_f in MOIU.eachscalar(f)
        g = GenericAffExpr{S,GenericVariableRef{T}}(scalar_f.constant)
        for t in scalar_f.terms
            add_to_expression!(
                g,
                t.coefficient,
                GenericVariableRef(model, t.variable),
            )
        end
        push!(ret, g)
    end
    return ret
end

# MOI.VectorQuadraticFunction

function moi_function_type(::Type{<:Vector{<:GenericQuadExpr{T}}}) where {T}
    return MOI.VectorQuadraticFunction{T}
end

moi_function(a::Vector{<:GenericQuadExpr}) = MOI.VectorQuadraticFunction(a)

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VectorQuadraticFunction{C}},
) where {C,T}
    S = promote_type(T, C)
    return Vector{GenericQuadExpr{S,GenericVariableRef{T}}}
end

function jump_function(
    model::GenericModel{T},
    f::MOI.VectorQuadraticFunction{C},
) where {C,T}
    S = promote_type(T, C)
    return GenericQuadExpr{S,GenericVariableRef{T}}[
        GenericQuadExpr{S,GenericVariableRef{T}}(model, f) for
        f in MOIU.eachscalar(f)
    ]
end

# MOI.VectorNonlinearFunction

function moi_function_type(::Type{<:AbstractVector{<:GenericNonlinearExpr}})
    return MOI.VectorNonlinearFunction
end

function moi_function(f::AbstractVector{<:GenericNonlinearExpr})
    return MOI.VectorNonlinearFunction(f)
end

function moi_function(
    model::GenericModel,
    f::AbstractVector{<:GenericNonlinearExpr},
)
    return MOI.VectorNonlinearFunction([moi_function(model, row) for row in f])
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VectorNonlinearFunction},
) where {T}
    return Vector{GenericNonlinearExpr{GenericVariableRef{T}}}
end

function jump_function(
    model::GenericModel{T},
    f::MOI.VectorNonlinearFunction,
) where {T}
    return GenericNonlinearExpr{GenericVariableRef{T}}[
        jump_function(model, fi) for fi in MOI.Utilities.eachscalar(f)
    ]
end

# MOI.Nonlinear.Expression

function jump_function(model::GenericModel, expr::MOI.Nonlinear.Expression)
    V = variable_ref_type(typeof(model))
    nlp = nonlinear_model(model)::MOI.Nonlinear.Model
    parsed = Vector{Any}(undef, length(expr.nodes))
    adj = MOI.Nonlinear.adjacency_matrix(expr.nodes)
    rowvals = SparseArrays.rowvals(adj)
    for i in length(expr.nodes):-1:1
        node = expr.nodes[i]
        parsed[i] = if node.type == MOI.Nonlinear.NODE_CALL_UNIVARIATE
            GenericNonlinearExpr{V}(
                nlp.operators.univariate_operators[node.index],
                parsed[rowvals[SparseArrays.nzrange(adj, i)[1]]],
            )
        elseif node.type == MOI.Nonlinear.NODE_CALL_MULTIVARIATE
            GenericNonlinearExpr{V}(
                nlp.operators.multivariate_operators[node.index],
                Any[parsed[rowvals[j]] for j in SparseArrays.nzrange(adj, i)],
            )
        elseif node.type == MOI.Nonlinear.NODE_MOI_VARIABLE
            V(model, MOI.VariableIndex(node.index))
        elseif node.type == MOI.Nonlinear.NODE_VALUE
            expr.values[node.index]
        else
            # node.type == MOI.Nonlinear.NODE_COMPARISON
            # node.type == MOI.Nonlinear.NODE_LOGIC
            # node.type == MOI.Nonlinear.NODE_PARAMETER
            # node.type == MOI.Nonlinear.NODE_SUBEXPRESSION
            error(
                """
                Encountered an unsupported node type `$(node.type)` when converting \
                a nonlinear expression to a JuMP expression.

                This conversion is not currently supported. Use the MOI \
                representation directly or reformulate the expression.
                """,
            )
        end
    end
    return parsed[1]
end

# AbstractConstraint

"""
    moi_function(constraint::AbstractConstraint)

Return the function of the constraint `constraint` in the function-in-set form
as a `MathOptInterface.AbstractFunction`.
"""
function moi_function(constraint::AbstractConstraint)
    return moi_function(jump_function(constraint))
end

function moi_function(model, constraint::AbstractConstraint)
    return moi_function(model, jump_function(constraint))
end

function moi_function(model::GenericModel, constraint::AbstractConstraint)
    return moi_function(model, jump_function(constraint))
end

function check_belongs_to_model(con::AbstractConstraint, model::AbstractModel)
    check_belongs_to_model(jump_function(con), model)
    return
end

"""
    jump_function(constraint::AbstractConstraint)

Return the function of the constraint `constraint` in the function-in-set form
as a `AbstractJuMPScalar` or `Vector{AbstractJuMPScalar}`.
"""
jump_function(constraint::AbstractConstraint) = constraint.func

# Base.Number

moi_function(x::Number) = x

check_belongs_to_model(::Number, ::AbstractModel) = nothing

jump_function(::GenericModel{T}, x::Number) where {T} = convert(T, x)

# Base.AbstractArray

# `moi_function(::Array)` would be ambiguous with
# `moi_function(AbstractArray{<:AbstractVariableRef})`
moi_function(x::AbstractArray) = moi_function.(x)

function moi_function(x::AbstractArray{AbstractJuMPScalar})
    return error(
        """
        Unable to convert an array of type `::$(typeof(x))` to an equivalent function
        in MathOptInterface because the array has the abstract element type
        `AbstractJuMPScalar`.

        To fix this error, convert every element in the array to the same concrete
        element type.

        For example, instead of:
        ```julia
        model = Model();
        @variable(model, x);
        y = AbstractJuMPScalar[x, sin(x)]
        @objective(model, Min, y)
        ```
        do
        ```julia
        @objective(model, Min, convert.(NonlinearExpr, y))
        ```
        """,
    )
end

function check_belongs_to_model(f::AbstractArray, model::AbstractModel)
    for func in f
        check_belongs_to_model(func, model)
    end
    return
end
