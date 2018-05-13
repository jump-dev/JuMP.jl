#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/affexpr.jl
# Defines all types relating to affine expressions
# - GenericAffExpr              ∑ aᵢ xᵢ  +  c
#   - AffExpr                   Alias for (Float64, VariableRef)
#   - AffExprConstraint         AffExpr-in-set constraint
# Operator overloads in src/operators.jl
#############################################################################

# Utilities for OrderedDict
function add_or_set!(dict::OrderedDict{K,V}, k::K, v::V) where {K,V}
    # TODO: This unnecessarily requires two lookups for k.
    # TODO: Decide if we want to drop zeros here.
    dict[k] = get!(dict, k, zero(V)) + v
    return dict
end

function new_ordered_dict(::Type{K}, ::Type{V}, kv::AbstractArray{<:Pair}) where {K,V}
    dict = OrderedDict{K,V}()
    sizehint!(dict, length(kv))
    for pair in kv
        add_or_set!(dict, convert(K, pair.first), convert(V, pair.second))
    end
    return dict
end

function new_ordered_dict(::Type{K}, ::Type{V}, kv::Pair...) where {K,V}
    dict = OrderedDict{K,V}()
    sizehint!(dict, length(kv))
    for pair in kv
        add_or_set!(dict, convert(K, pair.first), convert(V, pair.second))
    end
    return dict
end



#############################################################################
# GenericAffExpr
# ∑ aᵢ xᵢ  +  c
mutable struct GenericAffExpr{CoefType,VarType} <: AbstractJuMPScalar
    constant::CoefType
    terms::OrderedDict{VarType,CoefType}
end

variablereftype(::GenericAffExpr{C, V}) where {C, V} = V

function GenericAffExpr(constant::V, kv::AbstractArray{Pair{K,V}}) where {K,V}
    result = GenericAffExpr{V,K}(constant, new_ordered_dict(K, V, kv))
end

function GenericAffExpr(constant::V, kv::Pair{K,V}...) where {K,V}
    return GenericAffExpr{V,K}(constant, new_ordered_dict(K, V, kv...))
end

function GenericAffExpr{V,K}(constant, kv::AbstractArray{<:Pair}) where {K,V}
    return GenericAffExpr{V,K}(convert(V, constant), new_ordered_dict(K, V, kv))
end

function GenericAffExpr{V,K}(constant, kv::Pair...) where {K,V}
    return GenericAffExpr{V,K}(convert(V, constant), new_ordered_dict(K, V, kv...))
end

Base.iszero(a::GenericAffExpr) = isempty(a.terms) && iszero(a.constant)
Base.zero(::Type{GenericAffExpr{C,V}}) where {C,V} = GenericAffExpr{C,V}(zero(C), OrderedDict{V,C}())
Base.one(::Type{GenericAffExpr{C,V}}) where {C,V}  = GenericAffExpr{C,V}(one(C), OrderedDict{V,C}())
Base.zero(a::GenericAffExpr) = zero(typeof(a))
Base.one( a::GenericAffExpr) =  one(typeof(a))
Base.copy(a::GenericAffExpr) = GenericAffExpr(copy(a.constant), copy(a.terms))

GenericAffExpr{C, V}() where {C, V} = zero(GenericAffExpr{C, V})

function map_coefficients!(a::GenericAffExpr, f::Function)
    # The iterator remains valid if existing elements are updated.
    for (coef, var) in linearterms(a)
        a.terms[var] = f(coef)
    end
    a.constant = f(a.constant)
    return a
end

function map_coefficients(a::GenericAffExpr, f::Function)
    return map_coefficients!(copy(a), f)
end

Base.sizehint!(a::GenericAffExpr, n::Int) = sizehint!(a.terms, n)

"""
    value(a::GenericAffExpr, map::Function)

Evaluate `a` given the value `map(v)` for each variable `v`.
"""
function value(a::GenericAffExpr{T, V}, map::Function) where {T, V}
    S = Base.promote_op(map, V)
    U = Base.promote_op(*, T, S)
    ret = U(a.constant)
    for (var, coef) in a.terms
        ret += coef * map(var)
    end
    ret
end

# Iterator protocol - iterates over tuples (aᵢ,xᵢ)
struct LinearTermIterator{GAE<:GenericAffExpr}
    aff::GAE
end

"""
    linearterms(aff::GenericAffExpr{C, V})

Provides an iterator over coefficient-variable tuples `(a_i::C, x_i::V)` in the
linear part of the affine expression.
"""
linearterms(aff::GenericAffExpr) = LinearTermIterator(aff)

reorder_iterator(p::Pair, state::Int) = ((p.second, p.first), state)

Base.start(lti::LinearTermIterator) = start(lti.aff.terms)
Base.done( lti::LinearTermIterator, state::Int) = done(lti.aff.terms, state)
Base.next( lti::LinearTermIterator, state::Int) = reorder_iterator(next(lti.aff.terms, state)...)
Base.length(lti::LinearTermIterator) = length(lti.aff.terms)

# TODO: Consider renaming. push! is an unusual name to update a term in a
# dictionary.
"""
    Base.push!{C,V}(aff::GenericAffExpr{C,V}, new_coef::C, new_var::V)

An efficient way to add to an affine expression in place. For example, to add
`5x` to an existing expression `aff`, use `push!(aff, 5.0, x)`. This is
*significantly* more efficient than `aff += 5.0*x`.
"""
function Base.push!(aff::GenericAffExpr{C,V}, new_coef::C, new_var::V) where {C,V}
    add_or_set!(aff.terms, new_var, new_coef)
    aff
end

# TODO: Consider renaming. append! is an unusual name to add two expressions.
"""
    Base.append!{C,V}(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V})

Efficiently append the terms of an affine expression to an existing affine
expression. For example, given `aff = 5.0*x` and `other = 7.0*y + 3.0*z`,
we can add to `aff` in place by using `append!(aff, other)`.This results in
`aff` equal to `5x + 7y + 3z`. `append!` is *significantly* more efficient than
`aff += other`.
"""
function Base.append!(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V}) where {C,V}
    merge!(+, aff.terms, other.terms)
    aff.constant += other.constant
    aff
end
# For consistency, allow appending constants and individual variables
Base.append!(aff::GenericAffExpr{C,C}, other::C) where {C} = error() # for ambiguity
function Base.append!(aff::GenericAffExpr{C,V}, other::C) where {C,V}
    aff.constant += other
    aff
end
function Base.append!(aff::GenericAffExpr{C,V}, other::Real) where {C,V}
    aff.constant += other
    aff
end
Base.append!(aff::GenericAffExpr{C,V}, other::V) where {C,V} = push!(aff,one(C),other)

function Base.isequal(aff::GenericAffExpr{C,V},other::GenericAffExpr{C,V}) where {C,V}
    return isequal(aff.constant, other.constant) && isequal(aff.terms, other.terms)
end

function Base.dropzeros(aff::GenericAffExpr)
    result = copy(aff)
    for (coef, var) in linearterms(aff)
        if iszero(coef)
            delete!(result.terms, var)
        end
    end
    return result
end

# Check if two AffExprs are equal after dropping zeros and disregarding the
# order. Mostly useful for testing.
function isequal_canonical(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V}) where {C,V}
    aff_nozeros = dropzeros(aff)
    other_nozeros = dropzeros(other)
    # Note: This depends on equality of OrderedDicts ignoring order.
    # This is the current behavior, but it seems questionable.
    return isequal(aff_nozeros.terms, other_nozeros.terms) && aff.constant == other.constant
end

Base.convert(::Type{GenericAffExpr{T,V}}, v::V)    where {T,V} = GenericAffExpr(zero(T), v => one(T))
Base.convert(::Type{GenericAffExpr{T,V}}, v::Real) where {T,V} = GenericAffExpr{T,V}(convert(T, v))

# Alias for (Float64, VariableRef), the specific GenericAffExpr used by JuMP
const AffExpr = GenericAffExpr{Float64,VariableRef}

# Check all coefficients are finite, i.e. not NaN, not Inf, not -Inf
function assert_isfinite(a::AffExpr)
    for (coef, var) in linearterms(a)
        isfinite(coef) || error("Invalid coefficient $coef on variable $var.")
    end
end

"""
    resultvalue(v::GenericAffExpr)

Evaluate an `GenericAffExpr` given the result returned by a solver.
Replaces `getvalue` for most use cases.
"""
resultvalue(a::GenericAffExpr) = value(a, resultvalue)

# Note: No validation is performed that the variables in the AffExpr belong to
# the same model.
function MOI.ScalarAffineFunction(a::AffExpr)
    vars = MOIVAR[]
    coefs = Float64[]
    sizehint!(vars, length(a.terms))
    sizehint!(coefs, length(a.terms))
    for (coef, var) in linearterms(a)
        push!(vars, index(var))
        push!(coefs, coef)
    end
    return MOI.ScalarAffineFunction(vars, coefs, a.constant)
end

function AffExpr(m::Model, f::MOI.ScalarAffineFunction)
    aff = AffExpr()
    for i in 1:length(f.variables)
        push!(aff, f.coefficients[i], VariableRef(m,f.variables[i]))
    end
    aff.constant = f.constant
    return aff
end

"""
    _fillvaf!(outputindex, variables, coefficients, offset::Int, oi::Int, aff::AffExpr)

Fills the vectors outputindex, variables, coefficients at indices starting at `offset+1` with the terms of `aff`.
The output index for all terms is `oi`.
"""
function _fillvaf!(outputindex, variables, coefficients, offset::Int, oi::Int, aff::AffExpr)
    i = 1
    for (coef, var) in linearterms(aff)
        outputindex[offset+i] = oi
        variables[offset+i] = index(var)
        coefficients[offset+i] = coef
        i += 1
    end
    offset + length(aff.terms)
end

function MOI.VectorAffineFunction(affs::Vector{AffExpr})
    len = sum(aff -> length(linearterms(aff)), affs)
    outputindex = Vector{Int}(len)
    variables = Vector{MOIVAR}(len)
    coefficients = Vector{Float64}(len)
    constant = Vector{Float64}(length(affs))
    offset = 0
    for (i, aff) in enumerate(affs)
        constant[i] = aff.constant
        offset = _fillvaf!(outputindex, variables, coefficients, offset, i, aff)
    end
    MOI.VectorAffineFunction(outputindex, variables, coefficients, constant)
end

function setobjective(m::Model, sense::Symbol, a::AffExpr)
    if sense == :Min
        moisense = MOI.MinSense
    else
        @assert sense == :Max
        moisense = MOI.MaxSense
    end
    MOI.set!(m.moibackend, MOI.ObjectiveSense(), moisense)
    MOI.set!(m.moibackend, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(a))
    nothing
end

"""
    objectivefunction(m::Model, ::Type{AffExpr})

Return an `AffExpr` object representing the objective function.
Error if the objective is not linear.
"""
function objectivefunction(m::Model, ::Type{AffExpr})
    f = MOI.get(m.moibackend, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())::MOI.ScalarAffineFunction
    return AffExpr(m, f)
end


# Copy an affine expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(a::GenericAffExpr, new_model::Model)
    result = zero(a)
    for (coef, var) in linearterms(a)
        push!(result, coef, copy(v, new_model))
    end
    result.constant = a.constant
    return result
end

# TODO GenericAffExprConstraint

struct AffExprConstraint{S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::AffExpr
    set::S
end

moi_function_and_set(c::AffExprConstraint) = (MOI.ScalarAffineFunction(c.func), c.set)

# TODO: Find somewhere to put this error message.
#addconstraint(m::Model, c::Array{AffExprConstraint}) =
#    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")

struct VectorAffExprConstraint{S <: MOI.AbstractVectorSet} <: AbstractConstraint
    func::Vector{AffExpr}
    set::S
end

moi_function_and_set(c::VectorAffExprConstraint) = (MOI.VectorAffineFunction(c.func), c.set)

function constraintobject(cr::ConstraintRef{Model}, ::Type{AffExpr}, ::Type{SetType}) where {SetType <: MOI.AbstractScalarSet}
    f = MOI.get(cr.m, MOI.ConstraintFunction(), cr)::MOI.ScalarAffineFunction
    s = MOI.get(cr.m, MOI.ConstraintSet(), cr)::SetType
    return AffExprConstraint(AffExpr(cr.m, f), s)
end

function constraintobject(cr::ConstraintRef{Model}, ::Type{Vector{AffExpr}}, ::Type{SetType}) where {SetType <: MOI.AbstractVectorSet}
    m = cr.m
    f = MOI.get(m, MOI.ConstraintFunction(), cr)::MOI.VectorAffineFunction
    s = MOI.get(m, MOI.ConstraintSet(), cr)::SetType
    return VectorAffExprConstraint(map(f -> AffExpr(m, f), MOIU.eachscalar(f)), s)
end
