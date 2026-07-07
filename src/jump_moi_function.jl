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
    moi_function([model::GenericModel,] x::AbstractJuMPScalar)
    moi_function([model::GenericModel,] x::AbstractArray{<:AbstractJuMPScalar})

Given a JuMP object `x`, return the MathOptInterface equivalent.

If a `model` is passed, this function may, as a performance optimization, cache
the mapping of `x` to the MOI equivalent in order to improve handling of common
subexpressions.

See also: [`jump_function`](@ref).

!!! compat
    Using the `model` argument requires or JuMP v1.31 later.

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

# A default fallback for backwards compatibility. The first argument `model` was
# introduced in JuMP@1.31.0.
moi_function(model, f) = moi_function(f)

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
    key = objectid(f)
    if haskey(model.subexpressions, key)
        return model.subexpressions[key]
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
        arg_key = objectid(arg)
        if haskey(model.subexpressions, arg_key)
            parent.args[i] = model.subexpressions[arg_key]
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
        model.subexpressions[arg_key] = child
    end
    model.subexpressions[key] = ret
    return ret
end

# A backwards-compatible function to preserve behavior prior to #4032. This was
# used by Plasmo.jl
function moi_function(f::GenericNonlinearExpr{V}) where {V}
    ret = MOI.ScalarNonlinearFunction(f.head, similar(f.args))
    stack = Tuple{MOI.ScalarNonlinearFunction,Int,GenericNonlinearExpr{V}}[]
    for i in length(f.args):-1:1
        if f.args[i] isa GenericNonlinearExpr{V}
            push!(stack, (ret, i, f.args[i]))
        else
            ret.args[i] = moi_function(f.args[i])
        end
    end
    while !isempty(stack)
        parent, i, arg = pop!(stack)
        child = MOI.ScalarNonlinearFunction(arg.head, similar(arg.args))
        parent.args[i] = child
        for j in length(arg.args):-1:1
            if arg.args[j] isa GenericNonlinearExpr{V}
                push!(stack, (child, j, arg.args[j]))
            else
                child.args[j] = moi_function(arg.args[j])
            end
        end
    end
    return ret
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

"""
    jump_function(constraint::AbstractConstraint)

Return the function of the constraint `constraint` in the function-in-set form
as a `AbstractJuMPScalar` or `Vector{AbstractJuMPScalar}`.
"""
jump_function(constraint::AbstractConstraint) = constraint.func

# Base.Number

moi_function(x::Number) = x

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
