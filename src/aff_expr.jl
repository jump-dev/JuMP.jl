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

function _add_or_set!(dict::OrderedDict{K,V}, k::K, v::V) where {K,V}
    # Adding zero terms to this dictionary leads to unacceptable performance
    # degradations. See, for example, https://github.com/jump-dev/JuMP.jl/issues/1946.
    if iszero(v)
        return dict  # No-op.
    end
    index = OrderedCollections.ht_keyindex2(dict, k)
    if index <= 0  # Key does not exist. We pay the penalty of a second lookup.
        setindex!(dict, v, k)
    else
        dict.vals[index] += v
        dict.keys[index] = k
    end
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

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> expr = x[2] + 3.0 * x[1] + 4.0
x[2] + 3 x[1] + 4

julia> expr.constant
4.0

julia> expr.terms
OrderedCollections.OrderedDict{VariableRef, Float64} with 2 entries:
  x[2] => 1.0
  x[1] => 3.0
```
"""
mutable struct GenericAffExpr{CoefType,VarType} <: AbstractJuMPScalar
    constant::CoefType
    terms::OrderedDict{VarType,CoefType}
end

variable_ref_type(::Type{GenericAffExpr{C,V}}) where {C,V} = V

function owner_model(x::GenericAffExpr)
    if !isempty(x.terms)
        return owner_model(first(keys(x.terms)))
    end
    return nothing
end

"""
    GenericAffExpr(constant::V, kv::AbstractArray{Pair{K,V}}) where {K,V}

Create a [`GenericAffExpr`](@ref) by passing a constant and a vector of pairs.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> GenericAffExpr(1.0, [x => 1.0])
x + 1
```
"""
function GenericAffExpr(constant::V, kv::AbstractArray{Pair{K,V}}) where {K,V}
    return GenericAffExpr{V,K}(constant, _new_ordered_dict(K, V, kv))
end

"""
    GenericAffExpr(constant::V, kv::Vararg{Pair{K,V},N}) where {K,V,N}

Create a [`GenericAffExpr`](@ref) by passing a constant and pairs of additional
arguments.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> GenericAffExpr(1.0, x => 1.0)
x + 1
```
"""
function GenericAffExpr(
    constant::V,
    kv1::Pair{K,V},
    tail::Vararg{Pair{K,V},N},
) where {K,V,N}
    return GenericAffExpr{V,K}(constant, _new_ordered_dict(K, V, kv1, tail...))
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

Base.zero(a::GenericAffExpr) = zero(typeof(a))

function Base.one(::Type{GenericAffExpr{C,V}}) where {C,V}
    return GenericAffExpr{C,V}(one(C), OrderedDict{V,C}())
end

function Base.oneunit(::Type{GenericAffExpr{C,V}}) where {C,V}
    return GenericAffExpr{C,V}(oneunit(C), OrderedDict{V,C}())
end

Base.one(a::GenericAffExpr) = one(typeof(a))

Base.oneunit(a::GenericAffExpr) = oneunit(typeof(a))

Base.copy(a::GenericAffExpr) = GenericAffExpr(copy(a.constant), copy(a.terms))

Base.broadcastable(a::GenericAffExpr) = Ref(a)

Base.conj(a::GenericAffExpr{<:Real}) = a
Base.real(a::GenericAffExpr{<:Real}) = a
Base.imag(a::GenericAffExpr{<:Real}) = zero(a)
Base.abs2(a::GenericAffExpr{<:Real}) = a^2
Base.isreal(x::GenericAffExpr{<:Real}) = true

Base.conj(a::GenericAffExpr{<:Complex}) = map_coefficients(conj, a)

function _map_coefs(f::Function, a::GenericAffExpr{Complex{T},V}) where {T,V}
    output = convert(GenericAffExpr{T,V}, f(a.constant))
    for (coef, var) in linear_terms(a)
        new_coef = f(coef)
        if !iszero(new_coef)
            output.terms[var] = new_coef
        end
    end
    return output
end

Base.real(a::GenericAffExpr{<:Complex}) = _map_coefs(real, a)
Base.imag(a::GenericAffExpr{<:Complex}) = _map_coefs(imag, a)
function Base.abs2(a::GenericAffExpr{<:Complex})
    imag_a = imag(a)
    return add_to_expression!(real(a)^2, imag_a, imag_a)
end

function Base.isreal(x::GenericAffExpr{<:Complex})
    return isreal(x.constant) && all(isreal, values(x.terms))
end

# Needed for cases when Julia uses `x == 0` instead of `iszero(x)` (for example,
# in the stdlib).
Base.:(==)(x::GenericAffExpr, y::Number) = isempty(x.terms) && x.constant == y

"""
    coefficient(a::GenericAffExpr{C,V}, v::V) where {C,V}

Return the coefficient associated with variable `v` in the affine expression `a`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> expr = 2.0 * x + 1.0;

julia> coefficient(expr, x)
2.0
```
"""
coefficient(a::GenericAffExpr{C,V}, v::V) where {C,V} = get(a.terms, v, zero(C))

coefficient(::GenericAffExpr{C,V}, ::V, ::V) where {C,V} = zero(C)

"""
    drop_zeros!(expr::GenericAffExpr)

Remove terms in the affine expression with `0` coefficients.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> expr = x[1] + x[2];

julia> add_to_expression!(expr, -1.0, x[1])
0 x[1] + x[2]

julia> drop_zeros!(expr)

julia> expr
x[2]
```
"""
function drop_zeros!(expr::GenericAffExpr)
    _drop_zeros!(expr.terms)
    return
end

GenericAffExpr{C,V}() where {C,V} = zero(GenericAffExpr{C,V})

"""
    map_coefficients_inplace!(f::Function, a::GenericAffExpr)

Apply `f` to the coefficients and constant term of an [`GenericAffExpr`](@ref)
`a` and update them in-place.

See also: [`map_coefficients`](@ref)

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

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

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> a = GenericAffExpr(1.0, x => 1.0)
x + 1

julia> map_coefficients(c -> 2 * c, a)
2 x + 2

julia> a
x + 1
```
"""
function map_coefficients(f::Function, a::GenericAffExpr)
    # `map_coefficients(f, a.constant)` infers the coefficient type
    # which is then picked up in the method signature of `_map_aff`
    # and then used to build the `OrderedDict`.
    return _map_aff(f, f(a.constant), a)
end
function _map_aff(f, constant::C, a::GenericAffExpr{T,V}) where {C,T,V}
    terms = OrderedDict{V,C}()
    for (coef, var) in linear_terms(a)
        terms[var] = f(coef)
    end
    return GenericAffExpr(constant, terms)
end

Base.sizehint!(a::GenericAffExpr, n::Int) = sizehint!(a.terms, n)

"""
    value(var_value::Function, ex::GenericAffExpr)

Evaluate `ex` using `var_value(v)` as the value for each variable `v`.
"""
function value(var_value::Function, ex::GenericAffExpr{T,V}) where {T,V}
    S = Base.promote_op(var_value, V)
    U = Base.promote_op(*, T, S)
    ret = convert(U, ex.constant)
    for (var, coef) in ex.terms
        ret += coef * var_value(var)
    end
    return ret
end

"""
    constant(aff::GenericAffExpr{C,V})::C

Return the constant of the affine expression.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> aff = 2.0 * x + 3.0;

julia> constant(aff)
3.0
```
"""
constant(aff::GenericAffExpr) = aff.constant

"""
    LinearTermIterator{GAE<:GenericAffExpr}

A struct that implements the `iterate` protocol in order to iterate over tuples
of `(coefficient, variable)` in the `GenericAffExpr`.
"""
struct LinearTermIterator{GAE<:GenericAffExpr}
    aff::GAE
end

"""
    linear_terms(aff::GenericAffExpr{C,V})

Provides an iterator over coefficient-variable tuples `(a_i::C, x_i::V)` in the
linear part of the affine expression.
"""
linear_terms(aff::GenericAffExpr) = LinearTermIterator(aff)

_reverse_pair_to_tuple(p::Pair) = (p.second, p.first)

function Base.iterate(lti::LinearTermIterator)
    ret = iterate(lti.aff.terms)
    if ret === nothing
        return
    else
        return _reverse_pair_to_tuple(ret[1]), ret[2]
    end
end

function Base.iterate(lti::LinearTermIterator, state)
    ret = iterate(lti.aff.terms, state)
    if ret === nothing
        return
    else
        return _reverse_pair_to_tuple(ret[1]), ret[2]
    end
end
Base.length(lti::LinearTermIterator) = length(lti.aff.terms)

function Base.eltype(::LinearTermIterator{GenericAffExpr{C,V}}) where {C,V}
    return Tuple{C,V}
end

"""
    add_to_expression!(expression, terms...)

Updates `expression` in-place to `expression + (*)(terms...)`.

This is typically much more efficient than `expression += (*)(terms...)` because
it avoids the temorary allocation of the right-hand side term.

For example, `add_to_expression!(expression, a, b)` produces the same result as
`expression += a*b`, and `add_to_expression!(expression, a)` produces the same result as
`expression += a`.

## When to implement

Only a few methods are defined, mostly for internal use, and only for the cases
when:

 1. they can be implemented efficiently
 2. `expression` is capable of storing the result. For example,
    `add_to_expression!(::AffExpr, ::GenericVariableRef, ::GenericVariableRef)`
    is not defined because a `GenericAffExpr` cannot store the product of two
    variables.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> expr = 2 + x
x + 2

julia> add_to_expression!(expr, 3, x)
4 x + 2

julia> expr
4 x + 2
```
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
    aff::GenericAffExpr{S,V},
    other::GenericAffExpr{T,V},
) where {S,T,V}
    # Note: merge!() doesn't appear to call sizehint!(). Is this important?
    merge!(+, aff.terms, other.terms)
    aff.constant += other.constant
    return aff
end

# With two factors.

function add_to_expression!(
    expr::GenericAffExpr{C,V},
    α::_Constant,
    β::_Constant,
) where {C,V}
    return add_to_expression!(expr, *(α, β))
end

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
    new_coef::V,
    new_var::V,
) where {C,V<:_Constant}
    _add_or_set!(aff.terms, new_var, convert(C, _constant_to_number(new_coef)))
    return aff
end

function add_to_expression!(
    aff::GenericAffExpr{S,V},
    coef::_Constant,
    other::GenericAffExpr{T,V},
) where {S,T,V}
    sizehint!(aff, length(linear_terms(aff)) + length(linear_terms(other)))
    for (term_coef, var) in linear_terms(other)
        _add_or_set!(aff.terms, var, convert(S, coef * term_coef))
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

function _drop_zeros!(terms::OrderedDict)
    for (var, coef) in terms
        if iszero(coef)
            delete!(terms, var)
        elseif coef isa Complex && iszero(imag(coef))
            terms[var] = real(coef)
        end
    end
    return
end

function SparseArrays.dropzeros(aff::GenericAffExpr)
    result = copy(aff)
    _drop_zeros!(result.terms)
    if iszero(result.constant)
        # This is to work around isequal(0.0, -0.0) == false.
        result.constant = zero(typeof(result.constant))
    elseif result.constant isa Complex && iszero(imag(result.constant))
        result.constant = real(result.constant)
    end
    return result
end

function isequal_canonical(
    aff::GenericAffExpr{C,V},
    other::GenericAffExpr{C,V},
) where {C,V}
    aff_nozeros = SparseArrays.dropzeros(aff)
    other_nozeros = SparseArrays.dropzeros(other)
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

function Base.convert(
    ::Type{GenericAffExpr{T,V}},
    aff::GenericAffExpr{T,V},
) where {T,V}
    return aff
end

function Base.convert(
    ::Type{GenericAffExpr{T,V}},
    aff::GenericAffExpr{S,V},
) where {S,T,V}
    return GenericAffExpr{T,V}(
        convert(T, aff.constant),
        convert(OrderedDict{V,T}, aff.terms),
    )
end

"""
    AffExpr

Alias for `GenericAffExpr{Float64,VariableRef}`, the specific
[`GenericAffExpr`](@ref) used by JuMP.
"""
const AffExpr = GenericAffExpr{Float64,VariableRef}

# Check all coefficients are finite, that is, not NaN, not Inf, not -Inf
function _assert_isfinite(a::GenericAffExpr)
    for (coef, var) in linear_terms(a)
        if !isfinite(coef)
            error("Invalid coefficient $coef on variable $var.")
        end
    end
    if isnan(a.constant)
        error(
            "Expression contains an invalid NaN constant. This could be " *
            "produced by `Inf - Inf`.",
        )
    end
    return
end

"""
    value(v::GenericAffExpr; result::Int = 1)

Return the value of the `GenericAffExpr` `v` associated with result index
`result` of the most-recent solution returned by the solver.

See also: [`result_count`](@ref).
"""
function value(a::GenericAffExpr; result::Int = 1)
    return value(a) do x
        return value(x; result = result)
    end
end

function check_belongs_to_model(a::GenericAffExpr, model::AbstractModel)
    for variable in keys(a.terms)
        check_belongs_to_model(variable, model)
    end
    return
end

# Note: No validation is performed that the variables in the AffExpr belong to
# the same model. The verification is done in `check_belongs_to_model` which
# should be called before calling `MOI.ScalarAffineFunction`.
function MOI.ScalarAffineFunction(
    a::GenericAffExpr{C,<:GenericVariableRef},
) where {C}
    _assert_isfinite(a)
    terms = MOI.ScalarAffineTerm{C}[
        MOI.ScalarAffineTerm(t[1], index(t[2])) for t in linear_terms(a)
    ]
    return MOI.ScalarAffineFunction(terms, a.constant)
end

"""
    moi_function(x::AbstractJuMPScalar)
    moi_function(x::AbstractArray{<:AbstractJuMPScalar})

Given a JuMP object `x`, return the MathOptInterface equivalent.

See also: [`jump_function`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> f = 2.0 * x + 1.0
2 x + 1

julia> moi_function(f)
1.0 + 2.0 MOI.VariableIndex(1)
```
"""
function moi_function end

function moi_function(x::AbstractArray{AbstractJuMPScalar})
    return error(
        "Unable to convert array of type `::$(typeof(x))` to an equivalent " *
        "function in MathOptInterface because the array has the abstract " *
        "element type `AbstractJuMPScalar`. To fix this error, convert every " *
        "element in the array to the same concrete element type.\n\n" *
        """For example, instead of:
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

moi_function(a::GenericAffExpr) = MOI.ScalarAffineFunction(a)

function moi_function_type(::Type{<:GenericAffExpr{T}}) where {T}
    return MOI.ScalarAffineFunction{T}
end

function GenericAffExpr{C,GenericVariableRef{T}}(
    m::GenericModel{T},
    f::MOI.ScalarAffineFunction,
) where {C,T}
    aff = GenericAffExpr{C,GenericVariableRef{T}}(f.constant)
    for t in f.terms
        add_to_expression!(
            aff,
            t.coefficient,
            GenericVariableRef(m, t.variable),
        )
    end
    return aff
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.ScalarAffineFunction{C}},
) where {C,T}
    return GenericAffExpr{C,GenericVariableRef{T}}
end

function jump_function(
    model::GenericModel{T},
    f::MOI.ScalarAffineFunction{C},
) where {C,T}
    return GenericAffExpr{C,GenericVariableRef{T}}(model, f)
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VectorAffineFunction{C}},
) where {C,T}
    return Vector{GenericAffExpr{C,GenericVariableRef{T}}}
end

function jump_function(
    model::GenericModel{T},
    f::MOI.VectorAffineFunction{C},
) where {T,C}
    ret = GenericAffExpr{C,GenericVariableRef{T}}[]
    for scalar_f in MOIU.eachscalar(f)
        g = GenericAffExpr{C,GenericVariableRef{T}}(scalar_f.constant)
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

"""
    _fill_vaf!(
        terms::Vector{<:MOI.VectorAffineTerm},
        offset::Int,
        oi::Int,
        aff::AbstractJuMPScalar,
    )

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

function MOI.VectorAffineFunction(
    affs::Vector{GenericAffExpr{C,GenericVariableRef{T}}},
) where {C,T}
    len = 0
    for aff in affs
        len += length(linear_terms(aff))
    end
    terms = Vector{MOI.VectorAffineTerm{C}}(undef, len)
    constant = Vector{C}(undef, length(affs))
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

"""
    _eval_as_variable(f::F, x::GenericAffExpr, args...) where {F}

In many cases, `@variable` can return a `GenericAffExpr` instead of a
`GenericVariableRef`. This is particularly the case for complex-valued
expressions. To make common operations like `lower_bound(x)` work, we should
forward the method if and only if `x` is convertable to a `GenericVariableRef`.
"""
function _eval_as_variable(f::F, x::GenericAffExpr, args...) where {F}
    if length(x.terms) != 1
        error(
            "Cannot call $f with $x because it is not an affine expression " *
            "of one variable.",
        )
    end
    variable, coefficient = first(x.terms)
    if !isone(coefficient)
        error(
            "Cannot call $f with $x because the variable has a coefficient " *
            "that is different to `+1`.",
        )
    end
    return f(variable, args...)
end

# start_value(::GenericAffExpr)

start_value(x::GenericAffExpr) = _eval_as_variable(start_value, x)

function set_start_value(x::GenericAffExpr, value)
    _eval_as_variable(set_start_value, x, value)
    return
end

# lower_bound(::GenericAffExpr)

has_lower_bound(x::GenericAffExpr) = _eval_as_variable(has_lower_bound, x)

lower_bound(x::GenericAffExpr) = _eval_as_variable(lower_bound, x)

delete_lower_bound(x::GenericAffExpr) = _eval_as_variable(delete_lower_bound, x)

function set_lower_bound(x::GenericAffExpr, value)
    return _eval_as_variable(set_lower_bound, x, value)
end

# upper_bound(::GenericAffExpr)

has_upper_bound(x::GenericAffExpr) = _eval_as_variable(has_upper_bound, x)

upper_bound(x::GenericAffExpr) = _eval_as_variable(upper_bound, x)

delete_upper_bound(x::GenericAffExpr) = _eval_as_variable(delete_upper_bound, x)

function set_upper_bound(x::GenericAffExpr, value)
    return _eval_as_variable(set_upper_bound, x, value)
end
