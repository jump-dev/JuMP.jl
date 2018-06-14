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
#   - AffExprConstraint         AffExpr-in-set constraint
# Operator overloads in src/operators.jl
#############################################################################

# Utilities for OrderedDict
function add_or_set!(dict::OrderedDict{K,V}, k::K, v::V) where {K,V}
    # TODO: This unnecessarily requires two lookups for k.
    # TODO: Decide if we want to drop zeros here after understanding the
    # performance implications.
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
    return GenericAffExpr{V,K}(constant, new_ordered_dict(K, V, kv))
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

function map_coefficients_inplace!(f::Function, a::GenericAffExpr)
    # The iterator remains valid if existing elements are updated.
    for (coef, var) in linearterms(a)
        a.terms[var] = f(coef)
    end
    a.constant = f(a.constant)
    return a
end

function map_coefficients(f::Function, a::GenericAffExpr)
    return map_coefficients_inplace!(f, copy(a))
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

"""
    add_to_expression!(expression, terms...)

Updates `expression` *in place* to `expression + (*)(terms...)`. This is
typically much more efficient than `expression += (*)(terms...)`. For example,
`add_to_expression!(expression, a, b)` produces the same result as `expression
+= a*b`, and `add_to_expression!(expression, a)` produces the same result as
`expression += a`.

Only a few methods are defined, mostly for internal use, and only for the cases
when (1) they can be implemented efficiently and (2) `expression` is capable of
storing the result. For example, `add_to_expression!(::AffExpr, ::VariableRef,
::VariableRef)` is not defined because a `GenericAffExpr` cannot store the
product of two variables.
"""
function add_to_expression! end

# TODO: add deprecations for Base.push! and Base.append!

function add_to_expression!(aff::GenericAffExpr{C,V}, new_coef::C, new_var::V) where {C,V}
    add_or_set!(aff.terms, new_var, new_coef)
    aff
end

function add_to_expression!(aff::GenericAffExpr{C,V}, new_var::V) where {C,V}
    add_or_set!(aff.terms, new_var, one(C))
    aff
end

function add_to_expression!(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V}) where {C,V}
    merge!(+, aff.terms, other.terms)
    aff.constant += other.constant
    aff
end

function add_to_expression!(aff::GenericAffExpr{C,V}, other::C) where {C,V}
    aff.constant += other
    aff
end
function add_to_expression!(aff::GenericAffExpr{C,V}, other::Real) where {C,V}
    aff.constant += other
    aff
end

function Base.isequal(aff::GenericAffExpr{C,V},other::GenericAffExpr{C,V}) where {C,V}
    return isequal(aff.constant, other.constant) && isequal(aff.terms, other.terms)
end

Base.hash(aff::GenericAffExpr, h::UInt) = hash(aff.constant, hash(aff.terms, h))

function Base.dropzeros(aff::GenericAffExpr)
    result = copy(aff)
    for (coef, var) in linearterms(aff)
        if iszero(coef)
            delete!(result.terms, var)
        end
    end
    if iszero(result.constant)
        # This is to work around isequal(0.0, -0.0) == false.
        result.constant = zero(typeof(result.constant))
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
    return isequal(aff_nozeros, other_nozeros)
end

Base.convert(::Type{GenericAffExpr{T,V}}, v::V)    where {T,V} = GenericAffExpr(zero(T), v => one(T))
Base.convert(::Type{GenericAffExpr{T,V}}, v::Real) where {T,V} = GenericAffExpr{T,V}(convert(T, v))

# Check all coefficients are finite, i.e. not NaN, not Inf, not -Inf
function assert_isfinite(a::GenericAffExpr)
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

# Note: No validation is performed that the variables in the GenericAffExpr
# belong to the same model.
function MOI.ScalarAffineFunction(a::GenericAffExpr)
    assert_isfinite(a)
    terms = map(t -> MOI.ScalarAffineTerm(t[1], index(t[2])), linearterms(a))
    return MOI.ScalarAffineFunction(terms, a.constant)
end

function GenericAffExpr{T, VariableRef{MT}}(m::MT, f::MOI.ScalarAffineFunction) where {T, MT}
    aff = GenericAffExpr{T, VariableRef{MT}}()
    for t in f.terms
        add_to_expression!(aff, t.coefficient, VariableRef(m, t.variable_index))
    end
    aff.constant = f.constant
    return aff
end

"""
    _fillvaf!(terms, offset::Int, oi::Int, aff::GenericAffExpr)

Fills the vectors terms at indices starting at `offset+1` with the terms of `aff`.
The output index for all terms is `oi`.
"""
function _fillvaf!(terms, offset::Int, oi::Int, aff::GenericAffExpr)
    i = 1
    for (coef, var) in linearterms(aff)
        terms[offset+i] = MOI.VectorAffineTerm(Int64(oi), MOI.ScalarAffineTerm(coef, index(var)))
        i += 1
    end
    offset + length(linearterms(aff))
end

function MOI.VectorAffineFunction(affs::Vector{<:GenericAffExpr})
    len = sum(aff -> length(linearterms(aff)), affs)
    terms = Vector{MOI.VectorAffineTerm{Float64}}(len)
    constant = Vector{Float64}(length(affs))
    offset = 0
    for (i, aff) in enumerate(affs)
        constant[i] = aff.constant
        offset = _fillvaf!(terms, offset, i, aff)
    end
    MOI.VectorAffineFunction(terms, constant)
end

function setobjective(m::Model, sense::Symbol, a::GenericAffExpr)
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
    objectivefunction(m::Model, ::Type{<:GenericAffExpr})

Return an `GenericAffExpr` object representing the objective function.
Error if the objective is not linear.
"""
function objectivefunction(m::Model, AffExprType::Type{<:GenericAffExpr})
    f = MOI.get(m.moibackend, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())::MOI.ScalarAffineFunction
    return AffExprType(m, f)
end


# Copy an affine expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(a::GenericAffExpr, new_model::Model)
    result = zero(a)
    for (coef, var) in linearterms(a)
        add_to_expression!(result, coef, copy(var, new_model))
    end
    result.constant = a.constant
    return result
end

# TODO GenericAffExprConstraint

struct AffExprConstraint{V <: AbstractVariableRef, S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::GenericAffExpr{Float64, V}
    set::S
end

moi_function_and_set(c::AffExprConstraint) = (MOI.ScalarAffineFunction(c.func), c.set)

# TODO: Find somewhere to put this error message.
#addconstraint(m::Model, c::Array{AffExprConstraint}) =
#    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")

struct VectorAffExprConstraint{V <: AbstractVariableRef, S <: MOI.AbstractVectorSet} <: AbstractConstraint
    func::Vector{GenericAffExpr{Float64, V}}
    set::S
end

moi_function_and_set(c::VectorAffExprConstraint) = (MOI.VectorAffineFunction(c.func), c.set)

function constraintobject(cr::ConstraintRef{<:Model}, AffExprType::Type{<:GenericAffExpr}, ::Type{SetType}) where {SetType <: MOI.AbstractScalarSet}
    f = MOI.get(cr.m, MOI.ConstraintFunction(), cr)::MOI.ScalarAffineFunction
    s = MOI.get(cr.m, MOI.ConstraintSet(), cr)::SetType
    return AffExprConstraint(AffExprType(cr.m, f), s)
end

function constraintobject(cr::ConstraintRef{<:Model}, ::Type{Vector{AffExprType}}, ::Type{SetType}) where {AffExprType <: GenericAffExpr, SetType <: MOI.AbstractVectorSet}
    m = cr.m
    f = MOI.get(m, MOI.ConstraintFunction(), cr)::MOI.VectorAffineFunction
    s = MOI.get(m, MOI.ConstraintSet(), cr)::SetType
    return VectorAffExprConstraint(map(f -> AffExprType(m, f), MOIU.eachscalar(f)), s)
end
