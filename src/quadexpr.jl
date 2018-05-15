#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/quadexpr.jl
# Defines all types relating to expressions with a quadratic and affine part
# - GenericQuadExpr             ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c
#   - QuadExpr                  Alias for (Float64, VariableRef)
# - QuadExprConstraint       ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c  in set
# Operator overloads in src/operators.jl
#############################################################################


struct UnorderedPair{T}
    a::T
    b::T
end

Base.hash(p::UnorderedPair, h::UInt) = hash(hash(p.a) + hash(p.b), h)
function Base.isequal(p1::UnorderedPair, p2::UnorderedPair)
    return (p1.a == p2.a && p1.b == p2.b) || (p1.a == p2.b && p1.b == p2.a)
end

# GenericQuadExpr
# ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c
mutable struct GenericQuadExpr{CoefType,VarType} <: AbstractJuMPScalar
    aff::GenericAffExpr{CoefType,VarType}
    terms::OrderedDict{UnorderedPair{VarType}, CoefType}
end

function GenericQuadExpr(aff::GenericAffExpr{V,K}, kv::AbstractArray{Pair{UnorderedPair{K},V}}) where {K,V}
    return GenericQuadExpr{V,K}(aff, new_ordered_dict(UnorderedPair{K}, V, kv))
end

function GenericQuadExpr(aff::GenericAffExpr{V,K}, kv::Pair{UnorderedPair{K},V}...) where {K,V}
    return GenericQuadExpr{V,K}(aff, new_ordered_dict(UnorderedPair{K}, V, kv...))
end

function GenericAffExpr{V,K}(aff::GenericAffExpr{V,K}, kv::AbstractArray{<:Pair}) where {K,V}
    return GenericQuadExpr{V,K}(aff, new_ordered_dict(UnorderedPair{K}, V, kv))
end

function GenericQuadExpr{V,K}(aff::GenericAffExpr{V,K}, kv::Pair...) where {K,V}
    return GenericQuadExpr{V,K}(aff, new_ordered_dict(UnorderedPair{K}, V, kv...))
end

Base.iszero(q::GenericQuadExpr) = isempty(q.terms) && iszero(q.aff)
function Base.zero(::Type{GenericQuadExpr{C,V}}) where {C,V}
    return GenericQuadExpr(zero(GenericAffExpr{C,V}), OrderedDict{UnorderedPair{V}, C}())
end
function Base.one(::Type{GenericQuadExpr{C,V}}) where {C,V}
    return GenericQuadExpr(one(GenericAffExpr{C,V}), OrderedDict{UnorderedPair{V}, C}())
end
Base.zero(q::GenericQuadExpr) = zero(typeof(q))
Base.one(q::GenericQuadExpr)  = one(typeof(q))
Base.copy(q::GenericQuadExpr) = GenericQuadExpr(copy(q.aff), copy(q.terms))

function map_coefficients_inplace!(f::Function, q::GenericQuadExpr)
    # The iterator remains valid if existing elements are updated.
    for (key, value) in q.terms
        q.terms[key] = f(value)
    end
    map_coefficients_inplace!(f, q.aff)
    return q
end

function map_coefficients(f::Function, q::GenericQuadExpr)
    return map_coefficients_inplace!(f, copy(q))
end

"""
    linearterms(aff::GenericQuadExpr{C, V})

Provides an iterator over tuples `(coefficient::C, variable::V)` in the
linear part of the quadratic expression.
"""
linearterms(quad::GenericQuadExpr) = LinearTermIterator(quad.aff)

struct QuadTermIterator{GQE<:GenericQuadExpr}
    quad::GQE
end

"""
    quadterms(quad::GenericQuadExpr{C, V})

Provides an iterator over tuples `(coefficient::C, var_1::V, var_2::V)` in the
quadratic part of the quadratic expression.
"""
quadterms(quad::GenericQuadExpr) = QuadTermIterator(quad)

function reorder_iterator(p::Pair{UnorderedPair{V},C}, state::Int) where {C,V}
    return ((p.second, p.first.a, p.first.b), state)
end

Base.start(qti::QuadTermIterator) = start(qti.quad.terms)
Base.done(qti::QuadTermIterator, state::Int) = done(qti.quad.terms, state)
Base.next(qti::QuadTermIterator, state::Int) = reorder_iterator(next(qti.quad.terms, state)...)
Base.length(qti::QuadTermIterator) = length(qti.quad.terms)

function add_to_expression!(quad::GenericQuadExpr{C,V}, new_coef::C, new_var1::V, new_var2::V) where {C,V}
    # Node: OrderedDict updates the *key* as well. That is, if there was a
    # previous value for UnorderedPair(new_var2, new_var1), it's key will now be
    # UnorderedPair(new_var1, new_var2) (because these are defined as equal).
    key = UnorderedPair(new_var1, new_var2)
    add_or_set!(quad.terms, key, new_coef)
    quad
end

function add_to_expression!(quad::GenericQuadExpr{C, V}, new_coef::C, new_var::V) where {C,V}
    add_to_expression!(quad.aff, new_coef, new_var)
    quad
end

function add_to_expression!(q::GenericQuadExpr{T,S}, other::GenericQuadExpr{T,S}) where {T,S}
    merge!(+, q.terms, other.terms)
    add_to_expression!(q.aff, other.aff)
    q
end

function assert_isfinite(q::GenericQuadExpr)
    assert_isfinite(q.aff)
    for (coef, var1, var2) in quadterms(q)
        isfinite(coef) || error("Invalid coefficient $coef on quadratic term $var1*$var2.")
    end
end

function Base.isequal(q::GenericQuadExpr{T,S}, other::GenericQuadExpr{T,S}) where {T,S}
    return isequal(q.aff,other.aff) && isequal(q.terms, other.terms)
end

function Base.dropzeros(quad::GenericQuadExpr)
    quad_terms = copy(quad.terms)
    for (key, value) in quad.terms
        if iszero(value)
            delete!(quad_terms, key)
        end
    end
    return GenericQuadExpr(dropzeros(quad.aff), quad_terms)
end

# Check if two QuadExprs are equal regardless of the order, and after dropping zeros.
# Mostly useful for testing.
function isequal_canonical(quad::GenericQuadExpr{CoefType,VarType}, other::GenericQuadExpr{CoefType,VarType}) where {CoefType,VarType}
    quad_nozeros = dropzeros(quad)
    other_nozeros = dropzeros(other)
    return isequal(quad_nozeros, other_nozeros)
end

# Alias for (Float64, VariableRef)
const QuadExpr = GenericQuadExpr{Float64,VariableRef}
function Base.convert(::Type{GenericQuadExpr{C, V}}, v::Union{Real,AbstractVariableRef,GenericAffExpr}) where {C, V}
    return GenericQuadExpr(convert(GenericAffExpr{C, V}, v))
end
GenericQuadExpr{C, V}() where {C, V} = zero(GenericQuadExpr{C, V})

function MOI.ScalarQuadraticFunction(q::QuadExpr)
    assert_isfinite(q)
    qterms = map(t -> MOI.ScalarQuadraticTerm(t[2] == t[3] ? 2t[1] : t[1],
                                              index(t[2]),
                                              index(t[3])), quadterms(q))
    moi_aff = MOI.ScalarAffineFunction(q.aff)
    return MOI.ScalarQuadraticFunction(moi_aff.terms,
                                       qterms, moi_aff.constant)
end

function QuadExpr(m::Model, f::MOI.ScalarQuadraticFunction)
    quad = QuadExpr(AffExpr(m, MOI.ScalarAffineFunction(f.affine_terms,
                                                        f.constant)))
    for t in f.quadratic_terms
        v1 = t.variable_index_1
        v2 = t.variable_index_2
        coef = t.coefficient
        if v1 == v2
            coef /= 2
        end
        add_to_expression!(quad, coef, VariableRef(m, v1), VariableRef(m, v2))
    end
    return quad
end

function setobjective(m::Model, sense::Symbol, a::QuadExpr)
    if sense == :Min
        moisense = MOI.MinSense
    else
        @assert sense == :Max
        moisense = MOI.MaxSense
    end
    MOI.set!(m.moibackend, MOI.ObjectiveSense(), moisense)
    MOI.set!(m.moibackend, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), MOI.ScalarQuadraticFunction(a))
    nothing
end

"""
    objectivefunction(m::Model, ::Type{QuadExpr})

Return a `QuadExpr` object representing the objective function.
Error if the objective is not quadratic.
"""
function objectivefunction(m::Model, ::Type{QuadExpr})
    f = MOI.get(m.moibackend, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}())::MOI.ScalarQuadraticFunction
    return QuadExpr(m, f)
end

# Copy a quadratic expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(q::GenericQuadExpr, new_model::Model)
    GenericQuadExpr(copy(q.qvars1, new_model), copy(q.qvars2, new_model),
                copy(q.qcoeffs), copy(q.aff, new_model))
end

# TODO: resultvalue for QuadExpr

##########################################################################
# TODO: GenericQuadExprConstraint

struct QuadExprConstraint{V <: AbstractVariableRef, S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::GenericQuadExpr{Float64, V}
    set::S
end

moi_function_and_set(c::QuadExprConstraint) = (MOI.ScalarQuadraticFunction(c.func), c.set)

function constraintobject(cr::ConstraintRef{Model}, ::Type{QuadExpr}, ::Type{SetType}) where {SetType <: MOI.AbstractScalarSet}
    m = cr.m
    f = MOI.get(m.moibackend, MOI.ConstraintFunction(), index(cr))::MOI.ScalarQuadraticFunction
    s = MOI.get(m.moibackend, MOI.ConstraintSet(), index(cr))::SetType
    return QuadExprConstraint(QuadExpr(m, f), s)
end

# TODO: VectorQuadExprConstraint
