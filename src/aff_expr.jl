#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
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
    # degradations. See, e.g., https://github.com/jump-dev/JuMP.jl/issues/1946.
    if iszero(v)
        return dict  # No-op.
    end
    # TODO: This unnecessarily requires two lookups for k.
    dict[k] = get!(dict, k, zero(V)) + v
    return dict
end

function _new_ordered_dict(
    ::Type{K},
    ::Type{V},
    kv::AbstractArray{<:Pair},
) where {K,V}
    dict = OrderedDict{K,V}()
    sizehint!(dict, length(kv))
    for pair in kv
        _add_or_set!(dict, convert(K, pair.first), convert(V, pair.second))
    end
    return dict
end

function _new_ordered_dict(
    ::Type{K},
    ::Type{V},
    kv::Vararg{Pair,N},
) where {K,V,N}
    dict = OrderedDict{K,V}()
    sizehint!(dict, length(kv))
    for pair in kv
        _add_or_set!(dict, convert(K, pair.first), convert(V, pair.second))
    end
    return dict
end
# Shortcut for one and two arguments to avoid creating an empty dict and add
# elements one by one with `JuMP._add_or_set!`
function _new_ordered_dict(::Type{K}, ::Type{V}, kv::Pair) where {K,V}
    return OrderedDict{K,V}(kv)
end
function _new_ordered_dict(
    ::Type{K},
    ::Type{V},
    kv1::Pair,
    kv2::Pair,
) where {K,V}
    if isequal(kv1.first, kv2.first)
        return OrderedDict{K,V}(kv1.first => kv1.second + kv2.second)
    else
        return OrderedDict{K,V}(kv1, kv2)
    end
end

# As `!isbits(VariableRef)`, creating a pair allocates, with this API, we avoid
# this allocation.
function _build_aff_expr(constant::V, coef::V, var::K) where {V,K}
    terms = OrderedDict{K,V}()
    terms[var] = coef
    return GenericAffExpr{V,K}(constant, terms)
end
function _build_aff_expr(
    constant::V,
    coef1::V,
    var1::K,
    coef2::V,
    var2::K,
) where {V,K}
    if isequal(var1, var2)
        return _build_aff_expr(constant, coef1 + coef2, var1)
    end
    terms = OrderedDict{K,V}()
    terms[var1] = coef1
    terms[var2] = coef2
    return GenericAffExpr{V,K}(constant, terms)
end

#############################################################################

"""
    mutable struct GenericAffExpr{CoefType,VarType} <: AbstractJuMPScalar
        constant::CoefType
        terms::OrderedDict{VarType,CoefType}
    end

An expression type representing an affine expression of the form:
``\\sum a_i x_i + c``.

## Fields

 * `.constant`: the constant `c` in the expression.
 * `.terms`: an `OrderedDict`, with keys of `VarType` and values of `CoefType`
   describing the sparse vector `a`.
"""
mutable struct GenericAffExpr{CoefType,VarType} <: AbstractJuMPScalar
    constant::CoefType
    terms::OrderedDict{VarType,CoefType}
end

"""
    variable_ref_type(::GenericAffExpr{C, V}) where {C, V}

A helper function used internally by JuMP and some JuMP extensions. Returns the
variable type `V` from a [`GenericAffExpr`](@ref)
"""
variable_ref_type(::GenericAffExpr{C,V}) where {C,V} = V

"""
    GenericAffExpr(constant::V, kv::AbstractArray{Pair{K,V}}) where {K,V}

Create a [`GenericAffExpr`](@ref) by passing a constant and a vector of pairs.

## Examples

```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> GenericAffExpr(1.0, [x => 1.0])
x + 1
"""
function GenericAffExpr(constant::V, kv::AbstractArray{Pair{K,V}}) where {K,V}
    return GenericAffExpr{V,K}(constant, _new_ordered_dict(K, V, kv))
end

"""
    GenericAffExpr(constant::V, kv::Vararg{Pair{K,V},N}) where {K,V,N}

Create a [`GenericAffExpr`](@Ref) by passing a constant and pairs of additional
arguments.

## Examples

```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> GenericAffExpr(1.0, x => 1.0)
x + 1
"""
function GenericAffExpr(constant::V, kv::Vararg{Pair{K,V},N}) where {K,V,N}
    return GenericAffExpr{V,K}(constant, _new_ordered_dict(K, V, kv...))
end

function GenericAffExpr{V,K}(constant, kv::AbstractArray{<:Pair}) where {K,V}
    return GenericAffExpr{V,K}(
        convert(V, constant),
        _new_ordered_dict(K, V, kv),
    )
end

function GenericAffExpr{V,K}(constant, kv::Vararg{Pair,N}) where {K,V,N}
    return GenericAffExpr{V,K}(
        convert(V, constant),
        _new_ordered_dict(K, V, kv...),
    )
end

function Base.iszero(expr::GenericAffExpr)
    return iszero(expr.constant) && all(iszero, values(expr.terms))
end
function Base.zero(::Type{GenericAffExpr{C,V}}) where {C,V}
    return GenericAffExpr{C,V}(zero(C), OrderedDict{V,C}())
end
function Base.one(::Type{GenericAffExpr{C,V}}) where {C,V}
    return GenericAffExpr{C,V}(one(C), OrderedDict{V,C}())
end
Base.zero(a::GenericAffExpr) = zero(typeof(a))
Base.one(a::GenericAffExpr) = one(typeof(a))
Base.copy(a::GenericAffExpr) = GenericAffExpr(copy(a.constant), copy(a.terms))
Base.broadcastable(a::GenericAffExpr) = Ref(a)

"""
    coefficient(a::GenericAffExpr{C,V}, v::V) where {C,V}

Return the coefficient associated with variable `v` in the affine expression `a`.
"""
coefficient(a::GenericAffExpr{C,V}, v::V) where {C,V} = get(a.terms, v, zero(C))
coefficient(a::GenericAffExpr{C,V}, v1::V, v2::V) where {C,V} = zero(C)

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

GenericAffExpr{C,V}() where {C,V} = zero(GenericAffExpr{C,V})

function _affine_coefficient(f::GenericAffExpr{C,V}, variable::V) where {C,V}
    return get(f.terms, variable, zero(C))
end

"""
    map_coefficients_inplace!(f::Function, a::GenericAffExpr)

Apply `f` to the coefficients and constant term of an [`GenericAffExpr`](@ref)
`a` and update them in-place.

See also: [`map_coefficients`](@ref)

## Example

```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> a = GenericAffExpr(1.0, x => 1.0)
x + 1

julia> map_coefficients_inplace!(c -> 2 * c, a)
2 x + 2

julia> a
2 x + 2
```
"""
function map_coefficients_inplace!(f::Function, a::GenericAffExpr)
    # The iterator remains valid if existing elements are updated.
    for (coef, var) in linear_terms(a)
        a.terms[var] = f(coef)
    end
    a.constant = f(a.constant)
    return a
end

"""
    map_coefficients(f::Function, a::GenericAffExpr)

Apply `f` to the coefficients and constant term of an [`GenericAffExpr`](@ref)
`a` and return a new expression.

See also: [`map_coefficients_inplace!`](@ref)

## Example

```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> a = GenericAffExpr(1.0, x => 1.0)
x + 1

julia> map_coefficients(c -> 2 * c, a)
2 x + 2

julia> a
x + 1
```
"""
function map_coefficients(f::Function, a::GenericAffExpr)
    return map_coefficients_inplace!(f, copy(a))
end

Base.sizehint!(a::GenericAffExpr, n::Int) = sizehint!(a.terms, n)

"""
    value(ex::GenericAffExpr, var_value::Function)

Evaluate `ex` using `var_value(v)` as the value for each variable `v`.
"""
function value(ex::GenericAffExpr{T,V}, var_value::Function) where {T,V}
    S = Base.promote_op(var_value, V)
    U = Base.promote_op(*, T, S)
    ret = convert(U, ex.constant)
    for (var, coef) in ex.terms
        ret += coef * var_value(var)
    end
    return ret
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
function Base.eltype(lti::LinearTermIterator{GenericAffExpr{C,V}}) where {C,V}
    return Tuple{C,V}
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

# With one factor.

function add_to_expression!(aff::GenericAffExpr, other::_Constant)
    aff.constant += _constant_to_number(other)
    return aff
end

function add_to_expression!(aff::GenericAffExpr{C,V}, new_var::V) where {C,V}
    _add_or_set!(aff.terms, new_var, one(C))
    return aff
end

function add_to_expression!(
    aff::GenericAffExpr{C,V},
    other::GenericAffExpr{C,V},
) where {C,V}
    # Note: merge!() doesn't appear to call sizehint!(). Is this important?
    merge!(+, aff.terms, other.terms)
    aff.constant += other.constant
    return aff
end

# With two factors.

function add_to_expression!(
    aff::GenericAffExpr{C,V},
    new_coef::_Constant,
    new_var::V,
) where {C,V}
    _add_or_set!(aff.terms, new_var, convert(C, _constant_to_number(new_coef)))
    return aff
end

function add_to_expression!(
    aff::GenericAffExpr{C,V},
    new_var::V,
    new_coef::_Constant,
) where {C,V}
    return add_to_expression!(aff, new_coef, new_var)
end

function add_to_expression!(
    aff::GenericAffExpr{C,V},
    coef::_Constant,
    other::GenericAffExpr{C,V},
) where {C,V}
    sizehint!(aff, length(linear_terms(aff)) + length(linear_terms(other)))
    for (term_coef, var) in linear_terms(other)
        _add_or_set!(aff.terms, var, coef * term_coef)
    end
    aff.constant += coef * other.constant
    return aff
end

function add_to_expression!(
    aff::GenericAffExpr{C,V},
    other::GenericAffExpr{C,V},
    coef::_Constant,
) where {C,V}
    return add_to_expression!(aff, coef, other)
end

function Base.isequal(
    aff::GenericAffExpr{C,V},
    other::GenericAffExpr{C,V},
) where {C,V}
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

"""
    isequal_canonical(
        aff::GenericAffExpr{C,V},
        other::GenericAffExpr{C,V}
    ) where {C,V}

Return `true` if `aff` is equal to `other` after dropping zeros and disregarding
the order. Mainly useful for testing.
"""
function isequal_canonical(
    aff::GenericAffExpr{C,V},
    other::GenericAffExpr{C,V},
) where {C,V}
    aff_nozeros = dropzeros(aff)
    other_nozeros = dropzeros(other)
    # Note: This depends on equality of OrderedDicts ignoring order.
    # This is the current behavior, but it seems questionable.
    return isequal(aff_nozeros, other_nozeros)
end

function Base.convert(::Type{GenericAffExpr{T,V}}, v::V) where {T,V}
    return GenericAffExpr(zero(T), v => one(T))
end
function Base.convert(::Type{GenericAffExpr{T,V}}, v::_Constant) where {T,V}
    return GenericAffExpr{T,V}(convert(T, _constant_to_number(v)))
end

"""
    AffExpr

Alias for `GenericAffExpr{Float64,VariableRef}`, the specific
[`GenericAffExpr`](@ref) used by JuMP.
"""
const AffExpr = GenericAffExpr{Float64,VariableRef}

# Check all coefficients are finite, i.e. not NaN, not Inf, not -Inf
function _assert_isfinite(a::AffExpr)
    for (coef, var) in linear_terms(a)
        isfinite(coef) || error("Invalid coefficient $coef on variable $var.")
    end
    if isnan(a.constant)
        error(
            "Expression contains an invalid NaN constant. This could be " *
            "produced by `Inf - Inf`.",
        )
    end
end

"""
    value(v::GenericAffExpr; result::Int = 1)

Return the value of the `GenericAffExpr` `v` associated with result index
`result` of the most-recent solution returned by the solver.

Replaces `getvalue` for most use cases.

See also: [`result_count`](@ref).
"""
function value(a::GenericAffExpr; result::Int = 1)
    return value(a, (x) -> value(x; result = result))
end

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
    terms = MOI.ScalarAffineTerm{Float64}[
        MOI.ScalarAffineTerm(t[1], index(t[2])) for t in linear_terms(a)
    ]
    return MOI.ScalarAffineFunction(terms, a.constant)
end

"""
    moi_function(x)

Given a JuMP object `x`, return the MathOptInterface equivalent.

See also: [`jump_function`](@ref).
"""
function moi_function end

"""
    moi_function_type(::Type{T}) where {T}

Given a JuMP object type `T`, return the MathOptInterface equivalent.

See also: [`jump_function_type`](@ref).
"""
function moi_function_type end

"""
    jump_function(x)

Given an MathOptInterface object `x`, return the JuMP equivalent.

See also: [`moi_function`](@ref).
"""
function jump_function end

"""
    jump_function_type(::Type{T}) where {T}

Given an MathOptInterface object type `T`, return the JuMP equivalent.

See also: [`moi_function_type`](@ref).
"""
function jump_function_type end

moi_function(a::GenericAffExpr) = MOI.ScalarAffineFunction(a)
function moi_function_type(::Type{<:GenericAffExpr{T}}) where {T}
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
function jump_function_type(
    ::Model,
    ::Type{MOI.ScalarAffineFunction{T}},
) where {T}
    return GenericAffExpr{T,VariableRef}
end
function jump_function(model::Model, f::MOI.ScalarAffineFunction{T}) where {T}
    return GenericAffExpr{T,VariableRef}(model, f)
end
function jump_function_type(
    ::Model,
    ::Type{MOI.VectorAffineFunction{T}},
) where {T}
    return Vector{GenericAffExpr{T,VariableRef}}
end
function jump_function(model::Model, f::MOI.VectorAffineFunction{T}) where {T}
    return GenericAffExpr{T,VariableRef}[
        GenericAffExpr{T,VariableRef}(model, f) for f in MOIU.eachscalar(f)
    ]
end

"""
    _fill_vaf!(terms::Vector{<:MOI.VectorAffineTerm}, offset::Int, oi::Int,
               aff::AbstractJuMPScalar)

Fills the vectors terms at indices starting at `offset+1` with the affine terms
of `aff`. The output index for all terms is `oi`. Return the index of the last
term added.
"""
function _fill_vaf!(
    terms::Vector{<:MOI.VectorAffineTerm},
    offset::Int,
    oi::Int,
    aff::AbstractJuMPScalar,
)
    i = 1
    for (coef, var) in linear_terms(aff)
        terms[offset+i] = MOI.VectorAffineTerm(
            Int64(oi),
            MOI.ScalarAffineTerm(coef, index(var)),
        )
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
    return MOI.VectorAffineFunction(terms, constant)
end
moi_function(a::Vector{<:GenericAffExpr}) = MOI.VectorAffineFunction(a)
function moi_function_type(::Type{<:Vector{<:GenericAffExpr{T}}}) where {T}
    return MOI.VectorAffineFunction{T}
end

# TODO: Find somewhere to put this error message.
#add_constraint(m::Model, c::Array{AffExprConstraint}) =
#    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")
