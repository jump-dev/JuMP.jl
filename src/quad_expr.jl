#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/quad_expr.jl
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
    return GenericQuadExpr{V,K}(aff, _new_ordered_dict(UnorderedPair{K}, V, kv))
end

function GenericQuadExpr(aff::GenericAffExpr{V,K}, kv::Pair{UnorderedPair{K},V}...) where {K,V}
    return GenericQuadExpr{V,K}(aff, _new_ordered_dict(UnorderedPair{K}, V, kv...))
end

function GenericAffExpr{V,K}(aff::GenericAffExpr{V,K}, kv::AbstractArray{<:Pair}) where {K,V}
    return GenericQuadExpr{V,K}(aff, _new_ordered_dict(UnorderedPair{K}, V, kv))
end

function GenericQuadExpr{V,K}(aff::GenericAffExpr{V,K}, kv::Pair...) where {K,V}
    return GenericQuadExpr{V,K}(aff, _new_ordered_dict(UnorderedPair{K}, V, kv...))
end

function Base.iszero(expr::GenericQuadExpr)
    return iszero(expr.aff) && all(iszero, values(expr.terms))
end
function Base.zero(::Type{GenericQuadExpr{C,V}}) where {C,V}
    return GenericQuadExpr(zero(GenericAffExpr{C,V}), OrderedDict{UnorderedPair{V}, C}())
end
function Base.one(::Type{GenericQuadExpr{C,V}}) where {C,V}
    return GenericQuadExpr(one(GenericAffExpr{C,V}), OrderedDict{UnorderedPair{V}, C}())
end
Base.zero(q::GenericQuadExpr) = zero(typeof(q))
Base.one(q::GenericQuadExpr)  = one(typeof(q))
Base.copy(q::GenericQuadExpr) = GenericQuadExpr(copy(q.aff), copy(q.terms))
Base.broadcastable(q::GenericQuadExpr) = Ref(q)

"""
    drop_zeros!(expr::GenericQuadExpr)

Remove terms in the quadratic expression with `0` coefficients.
"""
function drop_zeros!(expr::GenericQuadExpr)
    drop_zeros!(expr.aff)
    for (key, coef) in expr.terms
        if iszero(coef)
            delete!(expr.terms, key)
        end
    end
    return
end

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

function _affine_coefficient(f::GenericQuadExpr{C, V}, variable::V) where {C, V}
    return _affine_coefficient(f.aff, variable)
end

"""
    constant(aff::GenericQuadExpr{C, V})::C

Return the constant of the quadratic expression.
"""
constant(quad::GenericQuadExpr) = constant(quad.aff)

"""
    linear_terms(quad::GenericQuadExpr{C, V})

Provides an iterator over tuples `(coefficient::C, variable::V)` in the
linear part of the quadratic expression.
"""
linear_terms(quad::GenericQuadExpr) = LinearTermIterator(quad.aff)

struct QuadTermIterator{GQE<:GenericQuadExpr}
    quad::GQE
end

"""
    quad_terms(quad::GenericQuadExpr{C, V})

Provides an iterator over tuples `(coefficient::C, var_1::V, var_2::V)` in the
quadratic part of the quadratic expression.
"""
quad_terms(quad::GenericQuadExpr) = QuadTermIterator(quad)

function _reorder_and_flatten(p::Pair{<:UnorderedPair})
    return (p.second, p.first.a, p.first.b)
end
function Base.iterate(qti::QuadTermIterator)
    ret = iterate(qti.quad.terms)
    if ret === nothing
        return nothing
    else
        return _reorder_and_flatten(ret[1]), ret[2]
    end
end
function Base.iterate(qti::QuadTermIterator, state)
    ret = iterate(qti.quad.terms, state)
    if ret === nothing
        return nothing
    else
        return _reorder_and_flatten(ret[1]), ret[2]
    end
end
Base.length(qti::QuadTermIterator) = length(qti.quad.terms)
function Base.eltype(qti::QuadTermIterator{GenericQuadExpr{C, V}}
                    ) where {C, V}
    return Tuple{C, V, V}
end

# With one factor.

function add_to_expression!(quad::GenericQuadExpr{C}, other::C) where C
    return add_to_expression!(quad.aff, other)
end

function add_to_expression!(q::GenericQuadExpr{T,S},
                            other::GenericAffExpr{T,S}) where {T,S}
    add_to_expression!(q.aff, other)
    return q
end

function add_to_expression!(q::GenericQuadExpr{T,S},
                            other::GenericQuadExpr{T,S}) where {T,S}
    merge!(+, q.terms, other.terms)
    add_to_expression!(q.aff, other.aff)
    return q
end

# With two factors.

function add_to_expression!(quad::GenericQuadExpr{C, V},
                            new_coef::Real,
                            new_var::V) where {C,V}
    add_to_expression!(quad.aff, new_coef, new_var)
    return quad
end

function add_to_expression!(quad::GenericQuadExpr{C, V},
                            new_var::Union{V, GenericAffExpr{C, V}},
                            new_coef::Real) where {C,V}
    return add_to_expression!(quad, new_coef, new_var)
end

function add_to_expression!(quad::GenericQuadExpr{C},
                            new_coef::Real,
                            new_aff::GenericAffExpr{C}) where {C}
    add_to_expression!(quad.aff, new_coef, new_aff)
    return quad
end

function add_to_expression!(quad::GenericQuadExpr{C, V}, coef::Real,
                            other::GenericQuadExpr{C, V}) where {C, V}
    for (key, term_coef) in other.terms
        _add_or_set!(quad.terms, key, coef * term_coef)
    end
    return add_to_expression!(quad, coef, other.aff)
end

function add_to_expression!(quad::GenericQuadExpr{C, V},
                            other::GenericQuadExpr{C, V},
                            coef::Real) where {C, V}
    return add_to_expression!(quad, coef, other)
end

function add_to_expression!(quad::GenericQuadExpr{C},
                            var_1::AbstractVariableRef,
                            var_2::AbstractVariableRef) where {C}
    return add_to_expression!(quad, one(C), var_1, var_2)
end

function add_to_expression!(quad::GenericQuadExpr{C,V},
                            var::V,
                            aff::GenericAffExpr{C,V}) where {C,V}
    for (coef, term_var) in linear_terms(aff)
        key = UnorderedPair(var, term_var)
        _add_or_set!(quad.terms, key, coef)
    end
    return add_to_expression!(quad, var, aff.constant)
end

function add_to_expression!(quad::GenericQuadExpr{C,V},
                            aff::GenericAffExpr{C,V},
                            var::V) where {C,V}
    return add_to_expression!(quad, var, aff)
end

function add_to_expression!(quad::GenericQuadExpr{C,V},
                            lhs::GenericAffExpr{C,V},
                            rhs::GenericAffExpr{C,V}) where {C,V}
    lhs_length = length(linear_terms(lhs))
    rhs_length = length(linear_terms(rhs))

    # Quadratic terms
    for (lhscoef, lhsvar) in linear_terms(lhs)
        for (rhscoef, rhsvar) in linear_terms(rhs)
            add_to_expression!(quad, lhscoef*rhscoef, lhsvar, rhsvar)
        end
    end

    # Try to preallocate space for aff
    cur = length(linear_terms(quad))
    if !iszero(lhs.constant) && !iszero(rhs.constant)
        sizehint!(quad.aff, cur + lhs_length + rhs_length)
    elseif !iszero(lhs.constant)
        sizehint!(quad.aff, cur + rhs_length)
    elseif !iszero(rhs.constant)
        sizehint!(quad.aff, cur + lhs_length)
    end

    # [LHS constant] * [RHS linear terms]
    if !iszero(lhs.constant)
        c = lhs.constant
        for (rhscoef, rhsvar) in linear_terms(rhs)
            add_to_expression!(quad.aff, c*rhscoef, rhsvar)
        end
    end

    # [RHS constant] * [LHS linear terms]
    if !iszero(rhs.constant)
        c = rhs.constant
        for (lhscoef, lhsvar) in linear_terms(lhs)
            add_to_expression!(quad.aff, c*lhscoef, lhsvar)
        end
    end

    quad.aff.constant += lhs.constant * rhs.constant

    return quad
end

# With three factors.

function add_to_expression!(quad::GenericQuadExpr{C,V}, new_coef::C,
                            new_var1::V, new_var2::V) where {C,V}
    # Node: OrderedDict updates the *key* as well. That is, if there was a
    # previous value for UnorderedPair(new_var2, new_var1), it's key will now be
    # UnorderedPair(new_var1, new_var2) (because these are defined as equal).
    key = UnorderedPair(new_var1, new_var2)
    _add_or_set!(quad.terms, key, new_coef)
    return quad
end


function _assert_isfinite(q::GenericQuadExpr)
    _assert_isfinite(q.aff)
    for (coef, var1, var2) in quad_terms(q)
        isfinite(coef) || error("Invalid coefficient $coef on quadratic term $var1*$var2.")
    end
end

function Base.isequal(q::GenericQuadExpr{T,S}, other::GenericQuadExpr{T,S}) where {T,S}
    return isequal(q.aff,other.aff) && isequal(q.terms, other.terms)
end

Base.hash(quad::GenericQuadExpr, h::UInt) = hash(quad.aff, hash(quad.terms, h))

function SparseArrays.dropzeros(quad::GenericQuadExpr)
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
# Used in `JuMP._mul!`.
function Base.convert(::Type{T}, quad::GenericQuadExpr{T}) where T
    if !isempty(quad.terms)
        throw(InexactError(:convert, T, quad))
    end
    return convert(T, quad.aff)
end

function check_belongs_to_model(q::GenericQuadExpr, model::AbstractModel)
    check_belongs_to_model(q.aff, model)
    for variable_pair in keys(q.terms)
        check_belongs_to_model(variable_pair.a, model)
        check_belongs_to_model(variable_pair.b, model)
    end
end

"""
    _moi_quadratic_term(t::Tuple)

Return the MOI.ScalarQuadraticTerm for the quadratic term `t`, element of the
[`quad_terms`](@ref) iterator. Note that the `VariableRef`s are transformed
into `MOI.VariableIndex`s hence the owner model information is lost.
"""
function _moi_quadratic_term(t::Tuple)
    return MOI.ScalarQuadraticTerm(t[2] == t[3] ? 2t[1] : t[1], index(t[2]),
                                   index(t[3]))
end
function MOI.ScalarQuadraticFunction(q::QuadExpr)
    _assert_isfinite(q)
    qterms = MOI.ScalarQuadraticTerm{Float64}[_moi_quadratic_term(t)
                                              for t in quad_terms(q)]
    moi_aff = MOI.ScalarAffineFunction(q.aff)
    return MOI.ScalarQuadraticFunction(moi_aff.terms,
                                       qterms, moi_aff.constant)
end
function moi_function(aff::GenericQuadExpr)
    return MOI.ScalarQuadraticFunction(aff)
end
function moi_function_type(::Type{<:GenericQuadExpr{T}}) where T
    return MOI.ScalarQuadraticFunction{T}
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
function jump_function_type(::Model,
                            ::Type{MOI.ScalarQuadraticFunction{T}}) where T
    return GenericQuadExpr{T, VariableRef}
end
function jump_function(model::Model, f::MOI.ScalarQuadraticFunction{T}) where T
    return GenericQuadExpr{T, VariableRef}(model, f)
end
function jump_function_type(::Model,
                            ::Type{MOI.VectorQuadraticFunction{T}}) where T
    return Vector{GenericQuadExpr{T, VariableRef}}
end
function jump_function(model::Model, f::MOI.VectorQuadraticFunction{T}) where T
    return GenericQuadExpr{T, VariableRef}[
        GenericQuadExpr{T, VariableRef}(model, f) for f in MOIU.eachscalar(f)]
end

"""
    _fill_vqf!(terms::Vector{<:MOI.VectorQuadraticTerm}, offset::Int, oi::Int,
               quad::AbstractJuMPScalar)

Fills the vectors terms at indices starting at `offset+1` with the quadratic
terms of `quad`. The output index for all terms is `oi`. Return the index of the
last term added.
"""
function _fill_vqf!(terms::Vector{<:MOI.VectorQuadraticTerm}, offset::Int,
                    oi::Int, aff::AbstractJuMPScalar)
    i = 1
    for term in quad_terms(aff)
        terms[offset + i] = MOI.VectorQuadraticTerm(Int64(oi),
                                                    _moi_quadratic_term(term))
        i += 1
    end
    return offset + length(quad_terms(aff))
end

function MOI.VectorQuadraticFunction(quads::Vector{QuadExpr})
    num_quadratic_terms = sum(quad -> length(quad_terms(quad)), quads)
    quadratic_terms = Vector{MOI.VectorQuadraticTerm{Float64}}(undef,
                                                            num_quadratic_terms)
    num_lin_terms = sum(quad -> length(linear_terms(quad)), quads)
    lin_terms = Vector{MOI.VectorAffineTerm{Float64}}(undef, num_lin_terms)
    constants = Vector{Float64}(undef, length(quads))
    quad_offset = 0
    lin_offset = 0
    for (i, quad) in enumerate(quads)
        quad_offset = _fill_vqf!(quadratic_terms, quad_offset, i, quad)
        lin_offset = _fill_vaf!(lin_terms, lin_offset, i, quad)
        constants[i] = constant(quad)
    end
    MOI.VectorQuadraticFunction(lin_terms, quadratic_terms, constants)
end
moi_function(a::Vector{<:GenericQuadExpr}) = MOI.VectorQuadraticFunction(a)
function moi_function_type(::Type{<:Vector{<:GenericQuadExpr{T}}}) where {T}
    return MOI.VectorQuadraticFunction{T}
end


# Copy a quadratic expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(q::GenericQuadExpr, new_model::Model)
    GenericQuadExpr(copy(q.qvars1, new_model), copy(q.qvars2, new_model),
                    copy(q.qcoeffs), copy(q.aff, new_model))
end

# Requires that value_func(::VarType) is defined.
function value(ex::GenericQuadExpr{CoefType, VarType},
               value_func::Function) where {CoefType, VarType}
    RetType = Base.promote_op(
        (ctype, vtype) -> ctype * value_func(vtype) * value_func(vtype),
        CoefType, VarType)
    ret = convert(RetType, value(ex.aff, value_func))
    for (vars, coef) in ex.terms
        ret += coef * value_func(vars.a) * value_func(vars.b)
    end
    return ret
end

JuMP.value(ex::JuMP.GenericQuadExpr) = value(ex, JuMP.value)
