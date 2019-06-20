#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/aff_expr.jl
# Defines all types relating to affine expressions
# - GenericAffExpr              ∑ aᵢ xᵢ  +  c
#   - AffExpr                   Alias for (Float64, VariableRef)
#   - AffExprConstraint         AffExpr-in-set constraint
# Operator overloads in src/operators.jl
#############################################################################

# Utilities for OrderedDict
function _add_or_set!(dict::OrderedDict{K,V}, k::K, v::V) where {K,V}
    # Adding zero terms to this dictionary leads to unacceptable performance
    # degradations. See, e.g., https://github.com/JuliaOpt/JuMP.jl/issues/1946.
    if iszero(v)
        return dict  # No-op.
    end
    # TODO: This unnecessarily requires two lookups for k.
    dict[k] = get!(dict, k, zero(V)) + v
    return dict
end

function _new_ordered_dict(::Type{K}, ::Type{V}, kv::AbstractArray{<:Pair}) where {K,V}
    dict = OrderedDict{K,V}()
    sizehint!(dict, length(kv))
    for pair in kv
        _add_or_set!(dict, convert(K, pair.first), convert(V, pair.second))
    end
    return dict
end

function _new_ordered_dict(::Type{K}, ::Type{V}, kv::Pair...) where {K,V}
    dict = OrderedDict{K,V}()
    sizehint!(dict, length(kv))
    for pair in kv
        _add_or_set!(dict, convert(K, pair.first), convert(V, pair.second))
    end
    return dict
end
# Shortcut for one and two arguments to avoid creating an empty dict and add
# elements one by one with `JuMP._add_or_set!`
function _new_ordered_dict(::Type{K}, ::Type{V}, kv::Pair) where {K, V}
    return OrderedDict{K, V}(kv)
end
function _new_ordered_dict(::Type{K}, ::Type{V}, kv1::Pair, kv2::Pair) where {K, V}
    if isequal(kv1.first, kv2.first)
        return OrderedDict{K, V}(kv1.first => kv1.second + kv2.second)
    else
        return OrderedDict{K, V}(kv1, kv2)
    end
end



#############################################################################
# GenericAffExpr
# ∑ aᵢ xᵢ  +  c
mutable struct GenericAffExpr{CoefType,VarType} <: AbstractJuMPScalar
    constant::CoefType
    terms::OrderedDict{VarType,CoefType}
end

variable_ref_type(::GenericAffExpr{C, V}) where {C, V} = V

function GenericAffExpr(constant::V, kv::AbstractArray{Pair{K,V}}) where {K,V}
    return GenericAffExpr{V,K}(constant, _new_ordered_dict(K, V, kv))
end

function GenericAffExpr(constant::V, kv::Pair{K,V}...) where {K,V}
    return GenericAffExpr{V,K}(constant, _new_ordered_dict(K, V, kv...))
end

function GenericAffExpr{V,K}(constant, kv::AbstractArray{<:Pair}) where {K,V}
    return GenericAffExpr{V,K}(convert(V, constant), _new_ordered_dict(K, V, kv))
end

function GenericAffExpr{V,K}(constant, kv::Pair...) where {K,V}
    return GenericAffExpr{V,K}(convert(V, constant), _new_ordered_dict(K, V, kv...))
end

function Base.iszero(expr::GenericAffExpr)
    return iszero(expr.constant) && all(iszero, values(expr.terms))
end
Base.zero(::Type{GenericAffExpr{C,V}}) where {C,V} = GenericAffExpr{C,V}(zero(C), OrderedDict{V,C}())
Base.one(::Type{GenericAffExpr{C,V}}) where {C,V}  = GenericAffExpr{C,V}(one(C), OrderedDict{V,C}())
Base.zero(a::GenericAffExpr) = zero(typeof(a))
Base.one( a::GenericAffExpr) =  one(typeof(a))
Base.copy(a::GenericAffExpr) = GenericAffExpr(copy(a.constant), copy(a.terms))
Base.broadcastable(a::GenericAffExpr) = Ref(a)

"""
    drop_zeros!(expr::GenericAffExpr)

Remove terms in the affine expression with `0` coefficients.
"""
function drop_zeros!(expr::GenericAffExpr)
    for (key, coef) in expr.terms
        if iszero(coef)
            delete!(expr.terms, key)
        end
    end
    return
end

GenericAffExpr{C, V}() where {C, V} = zero(GenericAffExpr{C, V})

function _affine_coefficient(f::GenericAffExpr{C, V}, variable::V) where {C, V}
    return get(f.terms, variable, zero(C))
end

function map_coefficients_inplace!(f::Function, a::GenericAffExpr)
    # The iterator remains valid if existing elements are updated.
    for (coef, var) in linear_terms(a)
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
    value(ex::GenericAffExpr, var_value::Function)

Evaluate `ex` using `var_value(v)` as the value for each variable `v`.
"""
function value(ex::GenericAffExpr{T, V}, var_value::Function) where {T, V}
    S = Base.promote_op(var_value, V)
    U = Base.promote_op(*, T, S)
    ret = convert(U, ex.constant)
    for (var, coef) in ex.terms
        ret += coef * var_value(var)
    end
    ret
end

"""
    constant(aff::GenericAffExpr{C, V})::C

Return the constant of the affine expression.
"""
constant(aff::GenericAffExpr) = aff.constant

# Iterator protocol - iterates over tuples (aᵢ,xᵢ)
struct LinearTermIterator{GAE<:GenericAffExpr}
    aff::GAE
end

"""
    linear_terms(aff::GenericAffExpr{C, V})

Provides an iterator over coefficient-variable tuples `(a_i::C, x_i::V)` in the
linear part of the affine expression.
"""
linear_terms(aff::GenericAffExpr) = LinearTermIterator(aff)


_reverse_pair_to_tuple(p::Pair) = (p.second, p.first)
function Base.iterate(lti::LinearTermIterator)
    ret = iterate(lti.aff.terms)
    if ret === nothing
        return nothing
    else
        return _reverse_pair_to_tuple(ret[1]), ret[2]
    end
end
function Base.iterate(lti::LinearTermIterator, state)
    ret = iterate(lti.aff.terms, state)
    if ret === nothing
        return nothing
    else
        return _reverse_pair_to_tuple(ret[1]), ret[2]
    end
end
Base.length(lti::LinearTermIterator) = length(lti.aff.terms)
function Base.eltype(lti::LinearTermIterator{GenericAffExpr{C, V}}
                    ) where {C, V}
    return Tuple{C, V}
end

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

# With one factor.

function add_to_expression!(aff::GenericAffExpr{C,V},
                            other::Real) where {C,V}
    aff.constant += other
    return aff
end

function add_to_expression!(aff::GenericAffExpr{C,V}, new_var::V) where {C,V}
    _add_or_set!(aff.terms, new_var, one(C))
    return aff
end

function add_to_expression!(aff::GenericAffExpr{C,V},
                            other::GenericAffExpr{C,V}) where {C,V}
    # Note: merge!() doesn't appear to call sizehint!(). Is this important?
    merge!(+, aff.terms, other.terms)
    aff.constant += other.constant
    return aff
end

# With two factors.

function add_to_expression!(aff::GenericAffExpr{C,V}, new_coef::Real,
                            new_var::V) where {C,V}
    _add_or_set!(aff.terms, new_var, convert(C, new_coef))
    return aff
end

function add_to_expression!(aff::GenericAffExpr{C,V}, new_var::V,
                            new_coef::Real) where {C,V}
    return add_to_expression!(aff, new_coef, new_var)
end

function add_to_expression!(aff::GenericAffExpr{C,V}, coef::Real,
                            other::GenericAffExpr{C,V}) where {C,V}
    sizehint!(aff, length(linear_terms(aff)) + length(linear_terms(other)))
    for (term_coef, var) in linear_terms(other)
        _add_or_set!(aff.terms, var, coef * term_coef)
    end
    aff.constant += coef * other.constant
    return aff
end

function add_to_expression!(aff::GenericAffExpr{C,V},
                            other::GenericAffExpr{C,V},
                            coef::Real) where {C,V}
    return add_to_expression!(aff, coef, other)
end


function Base.isequal(aff::GenericAffExpr{C,V},
                      other::GenericAffExpr{C,V}) where {C,V}
    return isequal(aff.constant, other.constant) &&
        isequal(aff.terms, other.terms)
end

Base.hash(aff::GenericAffExpr, h::UInt) = hash(aff.constant, hash(aff.terms, h))

function SparseArrays.dropzeros(aff::GenericAffExpr)
    result = copy(aff)
    for (coef, var) in linear_terms(aff)
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
# Used in `JuMP._mul!`.
function Base.convert(::Type{T}, aff::GenericAffExpr{T}) where T
    if !isempty(aff.terms)
        throw(InexactError(:convert, T, aff))
    end
    return convert(T, aff.constant)
end

# Alias for (Float64, VariableRef), the specific GenericAffExpr used by JuMP
const AffExpr = GenericAffExpr{Float64,VariableRef}

# Check all coefficients are finite, i.e. not NaN, not Inf, not -Inf
function _assert_isfinite(a::AffExpr)
    for (coef, var) in linear_terms(a)
        isfinite(coef) || error("Invalid coefficient $coef on variable $var.")
    end
end

"""
    value(v::GenericAffExpr)

Evaluate an `GenericAffExpr` given the result returned by a solver.
Replaces `getvalue` for most use cases.
"""
value(a::GenericAffExpr) = value(a, value)

function check_belongs_to_model(a::GenericAffExpr, model::AbstractModel)
    for variable in keys(a.terms)
        check_belongs_to_model(variable, model)
    end
end

# Note: No validation is performed that the variables in the AffExpr belong to
# the same model. The verification is done in `check_belongs_to_model` which
# should be called before calling `MOI.ScalarAffineFunction`.
function MOI.ScalarAffineFunction(a::AffExpr)
    _assert_isfinite(a)
    terms = MOI.ScalarAffineTerm{Float64}[MOI.ScalarAffineTerm(t[1],
                                                               index(t[2]))
                                          for t in linear_terms(a)]
    return MOI.ScalarAffineFunction(terms, a.constant)
end
moi_function(a::GenericAffExpr) = MOI.ScalarAffineFunction(a)
function moi_function_type(::Type{<:GenericAffExpr{T}}) where T
    return MOI.ScalarAffineFunction{T}
end


function AffExpr(m::Model, f::MOI.ScalarAffineFunction)
    aff = AffExpr()
    for t in f.terms
        add_to_expression!(aff, t.coefficient, VariableRef(m, t.variable_index))
    end
    aff.constant = f.constant
    return aff
end
function jump_function_type(::Model,
                            ::Type{MOI.ScalarAffineFunction{T}}) where T
    return GenericAffExpr{T, VariableRef}
end
function jump_function(model::Model, f::MOI.ScalarAffineFunction{T}) where T
    return GenericAffExpr{T, VariableRef}(model, f)
end
function jump_function_type(::Model,
                            ::Type{MOI.VectorAffineFunction{T}}) where T
    return Vector{GenericAffExpr{T, VariableRef}}
end
function jump_function(model::Model, f::MOI.VectorAffineFunction{T}) where T
    return GenericAffExpr{T, VariableRef}[
        GenericAffExpr{T, VariableRef}(model, f) for f in MOIU.eachscalar(f)]
end

"""
    _fill_vaf!(terms::Vector{<:MOI.VectorAffineTerm}, offset::Int, oi::Int,
               aff::AbstractJuMPScalar)

Fills the vectors terms at indices starting at `offset+1` with the affine terms
of `aff`. The output index for all terms is `oi`. Return the index of the last
term added.
"""
function _fill_vaf!(terms::Vector{<:MOI.VectorAffineTerm}, offset::Int, oi::Int,
                    aff::AbstractJuMPScalar)
    i = 1
    for (coef, var) in linear_terms(aff)
        terms[offset+i] = MOI.VectorAffineTerm(Int64(oi), MOI.ScalarAffineTerm(coef, index(var)))
        i += 1
    end
    return offset + length(linear_terms(aff))
end

function MOI.VectorAffineFunction(affs::Vector{AffExpr})
    len = sum(aff -> length(linear_terms(aff)), affs)
    terms = Vector{MOI.VectorAffineTerm{Float64}}(undef, len)
    constant = Vector{Float64}(undef, length(affs))
    offset = 0
    for (i, aff) in enumerate(affs)
        constant[i] = aff.constant
        offset = _fill_vaf!(terms, offset, i, aff)
    end
    MOI.VectorAffineFunction(terms, constant)
end
moi_function(a::Vector{<:GenericAffExpr}) = MOI.VectorAffineFunction(a)
function moi_function_type(::Type{<:Vector{<:GenericAffExpr{T}}}) where {T}
    return MOI.VectorAffineFunction{T}
end

# Copy an affine expression to a new model by converting all the
# variables to the new model's variables
function Base.copy(a::GenericAffExpr, new_model::AbstractModel)
    result = zero(a)
    for (coef, var) in linear_terms(a)
        add_to_expression!(result, coef, copy(var, new_model))
    end
    result.constant = a.constant
    return result
end

# TODO: Find somewhere to put this error message.
#add_constraint(m::Model, c::Array{AffExprConstraint}) =
#    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")
