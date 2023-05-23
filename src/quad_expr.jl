#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# src/quad_expr.jl
# Defines all types relating to expressions with a quadratic and affine part
# - GenericQuadExpr             ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c
#   - QuadExpr                  Alias for (Float64, VariableRef)
# - QuadExprConstraint       ∑qᵢⱼ xᵢⱼ  +  ∑ aᵢ xᵢ  +  c  in set
# Operator overloads in src/operators.jl
#############################################################################

"""
    UnorderedPair(a::T, b::T)

A wrapper type used by [`GenericQuadExpr`](@ref) with fields `.a` and `.b`.
"""
struct UnorderedPair{T}
    a::T
    b::T
end

Base.hash(p::UnorderedPair, h::UInt) = hash(hash(p.a) + hash(p.b), h)
function Base.isequal(p1::UnorderedPair, p2::UnorderedPair)
    return (p1.a == p2.a && p1.b == p2.b) || (p1.a == p2.b && p1.b == p2.a)
end

"""
    mutable struct GenericQuadExpr{CoefType,VarType} <: AbstractJuMPScalar
        aff::GenericAffExpr{CoefType,VarType}
        terms::OrderedDict{UnorderedPair{VarType}, CoefType}
    end

An expression type representing an quadratic expression of the form:
``\\sum q_{i,j} x_i x_j + \\sum a_i x_i + c``.

## Fields

 * `.aff`: an [`GenericAffExpr`](@ref) representing the affine portion of the
   expression.
 * `.terms`: an `OrderedDict`, with keys of `UnorderedPair{VarType}` and
   values of `CoefType`, describing the sparse list of terms `q`.
"""
mutable struct GenericQuadExpr{CoefType,VarType} <: AbstractJuMPScalar
    aff::GenericAffExpr{CoefType,VarType}
    terms::OrderedDict{UnorderedPair{VarType},CoefType}
end

variable_ref_type(::Type{GenericQuadExpr{C,V}}) where {C,V} = V

"""
    GenericQuadExpr(
        aff::GenericAffExpr{V,K},
        kv::AbstractArray{Pair{UnorderedPair{K},V}}
    ) where {K,V}

Create a [`GenericQuadExpr`](@ref) by passing a [`GenericAffExpr`](@ref) and a
vector of ([`UnorderedPair`](@ref), coefficient) pairs.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> GenericQuadExpr(GenericAffExpr(1.0, x => 2.0), [UnorderedPair(x, x) => 3.0])
3 x² + 2 x + 1
```
"""
function GenericQuadExpr(
    aff::GenericAffExpr{V,K},
    kv::AbstractArray{Pair{UnorderedPair{K},V}},
) where {K,V}
    return GenericQuadExpr{V,K}(aff, _new_ordered_dict(UnorderedPair{K}, V, kv))
end

function GenericQuadExpr(
    aff::GenericAffExpr{V,K},
    kv::Pair{UnorderedPair{K},V}...,
) where {K,V}
    return GenericQuadExpr{V,K}(
        aff,
        _new_ordered_dict(UnorderedPair{K}, V, kv...),
    )
end

function GenericQuadExpr{V,K}(aff::GenericAffExpr{V,K}, kv::Pair...) where {K,V}
    return GenericQuadExpr{V,K}(
        aff,
        _new_ordered_dict(UnorderedPair{K}, V, kv...),
    )
end

function Base.iszero(expr::GenericQuadExpr)
    return iszero(expr.aff) && all(iszero, values(expr.terms))
end
function Base.zero(::Type{GenericQuadExpr{C,V}}) where {C,V}
    return GenericQuadExpr(
        zero(GenericAffExpr{C,V}),
        OrderedDict{UnorderedPair{V},C}(),
    )
end
function Base.one(::Type{GenericQuadExpr{C,V}}) where {C,V}
    return GenericQuadExpr(
        one(GenericAffExpr{C,V}),
        OrderedDict{UnorderedPair{V},C}(),
    )
end
Base.zero(q::GenericQuadExpr) = zero(typeof(q))
Base.one(q::GenericQuadExpr) = one(typeof(q))
Base.copy(q::GenericQuadExpr) = GenericQuadExpr(copy(q.aff), copy(q.terms))
Base.broadcastable(q::GenericQuadExpr) = Ref(q)

Base.conj(a::GenericQuadExpr{<:Real}) = a
Base.real(a::GenericQuadExpr{<:Real}) = a
Base.imag(a::GenericQuadExpr{<:Real}) = a
Base.isreal(::GenericQuadExpr{<:Real}) = true

Base.conj(a::GenericQuadExpr{<:Complex}) = map_coefficients(conj, a)
Base.real(a::GenericQuadExpr{<:Complex}) = map_coefficients(real, a)
Base.imag(a::GenericQuadExpr{<:Complex}) = map_coefficients(imag, a)

function Base.isreal(x::GenericQuadExpr{<:Complex})
    return isreal(x.aff) && all(isreal, values(x.terms))
end

# Needed for cases when Julia uses `x == 0` instead of `iszero(x)` (e.g., in the
# stdlib).
Base.:(==)(x::GenericQuadExpr, y::Number) = isempty(x.terms) && x.aff == y

"""
    coefficient(a::GenericAffExpr{C,V}, v1::V, v2::V) where {C,V}

Return the coefficient associated with the term `v1 * v2` in the quadratic expression `a`.

Note that `coefficient(a, v1, v2)` is the same as `coefficient(a, v2, v1)`.
"""
function coefficient(q::GenericQuadExpr{C,V}, v1::V, v2::V) where {C,V}
    return get(q.terms, UnorderedPair(v1, v2), zero(C))
end

"""
    coefficient(a::GenericQuadExpr{C,V}, v::V) where {C,V}

Return the coefficient associated with variable `v` in the affine component of `a`.
"""
coefficient(q::GenericQuadExpr{C,V}, v::V) where {C,V} = coefficient(q.aff, v)

"""
    drop_zeros!(expr::GenericQuadExpr)

Remove terms in the quadratic expression with `0` coefficients.
"""
function drop_zeros!(expr::GenericQuadExpr)
    drop_zeros!(expr.aff)
    _drop_zeros!(expr.terms)
    return
end

"""
    map_coefficients_inplace!(f::Function, a::GenericQuadExpr)

Apply `f` to the coefficients and constant term of an [`GenericQuadExpr`](@ref)
`a` and update them in-place.

See also: [`map_coefficients`](@ref)

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> a = @expression(model, x^2 + x + 1)
x² + x + 1

julia> map_coefficients_inplace!(c -> 2 * c, a)
2 x² + 2 x + 2

julia> a
2 x² + 2 x + 2
```
"""
function map_coefficients_inplace!(f::Function, q::GenericQuadExpr)
    # The iterator remains valid if existing elements are updated.
    for (key, value) in q.terms
        q.terms[key] = f(value)
    end
    map_coefficients_inplace!(f, q.aff)
    return q
end

"""
    map_coefficients(f::Function, a::GenericQuadExpr)

Apply `f` to the coefficients and constant term of an [`GenericQuadExpr`](@ref)
`a` and return a new expression.

See also: [`map_coefficients_inplace!`](@ref)

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> a = @expression(model, x^2 + x + 1)
x² + x + 1

julia> map_coefficients(c -> 2 * c, a)
2 x² + 2 x + 2

julia> a
x² + x + 1
```
"""
function map_coefficients(f::Function, q::GenericQuadExpr)
    # `map_coefficients(f, q.aff)` infers the coefficient type
    # which is then picked up in the method signature of `_map_quad`
    # and then used to build the `OrderedDict`.
    return _map_quad(f, map_coefficients(f, q.aff), q)
end
function _map_quad(
    f::Function,
    aff::GenericAffExpr{C,V},
    q::GenericQuadExpr,
) where {C,V}
    terms = OrderedDict{UnorderedPair{V},C}()
    sizehint!(terms, length(q.terms))
    for (key, value) in q.terms
        terms[key] = f(value)
    end
    return GenericQuadExpr(aff, terms)
end

function _affine_coefficient(f::GenericQuadExpr{C,V}, variable::V) where {C,V}
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

"""
    QuadTermIterator{GQE<:GenericQuadExpr}

A struct that implements the `iterate` protocol in order to iterate over tuples
of `(coefficient, variable, variable)` in the `GenericQuadExpr`.
"""
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
function Base.eltype(qti::QuadTermIterator{GenericQuadExpr{C,V}}) where {C,V}
    return Tuple{C,V,V}
end

# With one factor.

function add_to_expression!(quad::GenericQuadExpr, other::_Constant)
    add_to_expression!(quad.aff, other)
    return quad
end

function add_to_expression!(quad::GenericQuadExpr{C,V}, other::V) where {C,V}
    add_to_expression!(quad.aff, other)
    return quad
end

function add_to_expression!(
    q::GenericQuadExpr{T,S},
    other::GenericAffExpr{T,S},
) where {T,S}
    add_to_expression!(q.aff, other)
    return q
end

function add_to_expression!(
    q::GenericQuadExpr{T,V},
    other::GenericQuadExpr{S,V},
) where {T,S,V}
    merge!(+, q.terms, other.terms)
    add_to_expression!(q.aff, other.aff)
    return q
end

# With two factors.

function add_to_expression!(
    expr::GenericQuadExpr{C,V},
    α::_Constant,
    β::_Constant,
) where {C,V}
    return add_to_expression!(expr, *(α, β))
end

function add_to_expression!(
    quad::GenericQuadExpr{C,V},
    new_coef::_Constant,
    new_var::V,
) where {C,V}
    add_to_expression!(quad.aff, new_coef, new_var)
    return quad
end

function add_to_expression!(
    quad::GenericQuadExpr{C,V},
    new_var::Union{V,GenericAffExpr{C,V}},
    new_coef::_Constant,
) where {C,V}
    return add_to_expression!(quad, new_coef, new_var)
end

function add_to_expression!(
    quad::GenericQuadExpr{C,V},
    new_coef::V,
    new_var::V,
) where {C,V<:Union{Number,LinearAlgebra.UniformScaling}}
    add_to_expression!(quad.aff, new_coef, new_var)
    return quad
end

function add_to_expression!(
    quad::GenericQuadExpr,
    new_coef::_Constant,
    new_aff::GenericAffExpr,
)
    add_to_expression!(quad.aff, new_coef, new_aff)
    return quad
end

function add_to_expression!(
    quad::GenericQuadExpr{S,V},
    coef::_Constant,
    other::GenericQuadExpr{T,V},
) where {S,T,V}
    for (key, term_coef) in other.terms
        _add_or_set!(quad.terms, key, convert(S, coef * term_coef))
    end
    return add_to_expression!(quad, coef, other.aff)
end

function add_to_expression!(
    quad::GenericQuadExpr,
    other::GenericQuadExpr,
    coef::_Constant,
)
    return add_to_expression!(quad, coef, other)
end

function add_to_expression!(
    quad::GenericQuadExpr{C},
    var_1::AbstractVariableRef,
    var_2::AbstractVariableRef,
) where {C}
    return add_to_expression!(quad, one(C), var_1, var_2)
end

function add_to_expression!(
    quad::GenericQuadExpr{C,V},
    var::V,
    aff::GenericAffExpr{C,V},
) where {C,V}
    for (coef, term_var) in linear_terms(aff)
        key = UnorderedPair(var, term_var)
        _add_or_set!(quad.terms, key, coef)
    end
    return add_to_expression!(quad, var, aff.constant)
end

function add_to_expression!(
    quad::GenericQuadExpr{C,V},
    aff::GenericAffExpr{C,V},
    var::V,
) where {C,V}
    return add_to_expression!(quad, var, aff)
end

function add_to_expression!(
    quad::GenericQuadExpr{C,V},
    aff::GenericAffExpr{C,V},
    var::V,
) where {C,V<:Union{Number,LinearAlgebra.UniformScaling}}
    return add_to_expression!(quad, var, aff)
end

function add_to_expression!(
    quad::GenericQuadExpr{C,V},
    lhs::GenericAffExpr{S,V},
    rhs::GenericAffExpr{T,V},
) where {C,S,T,V}
    lhs_length = length(linear_terms(lhs))
    rhs_length = length(linear_terms(rhs))

    # Quadratic terms
    for (lhscoef, lhsvar) in linear_terms(lhs)
        for (rhscoef, rhsvar) in linear_terms(rhs)
            add_to_expression!(quad, lhscoef * rhscoef, lhsvar, rhsvar)
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
            add_to_expression!(quad.aff, c * rhscoef, rhsvar)
        end
    end

    # [RHS constant] * [LHS linear terms]
    if !iszero(rhs.constant)
        c = rhs.constant
        for (lhscoef, lhsvar) in linear_terms(lhs)
            add_to_expression!(quad.aff, c * lhscoef, lhsvar)
        end
    end

    quad.aff.constant += lhs.constant * rhs.constant

    return quad
end

# With three factors.

function add_to_expression!(
    quad::GenericQuadExpr{C,V},
    new_coef::_Constant,
    new_var1::V,
    new_var2::V,
) where {C,V}
    # Node: OrderedDict updates the *key* as well. That is, if there was a
    # previous value for UnorderedPair(new_var2, new_var1), it's key will now be
    # UnorderedPair(new_var1, new_var2) (because these are defined as equal).
    key = UnorderedPair(new_var1, new_var2)
    _add_or_set!(quad.terms, key, convert(C, new_coef))
    return quad
end

function _assert_isfinite(q::GenericQuadExpr)
    _assert_isfinite(q.aff)
    for (coef, var1, var2) in quad_terms(q)
        isfinite(coef) ||
            error("Invalid coefficient $coef on quadratic term $var1*$var2.")
    end
end

function Base.isequal(
    q::GenericQuadExpr{T,S},
    other::GenericQuadExpr{T,S},
) where {T,S}
    return isequal(q.aff, other.aff) && isequal(q.terms, other.terms)
end

Base.hash(quad::GenericQuadExpr, h::UInt) = hash(quad.aff, hash(quad.terms, h))

function SparseArrays.dropzeros(quad::GenericQuadExpr)
    quad_terms = copy(quad.terms)
    _drop_zeros!(quad_terms)
    return GenericQuadExpr(SparseArrays.dropzeros(quad.aff), quad_terms)
end

# Check if two QuadExprs are equal regardless of the order, and after dropping zeros.
# Mostly useful for testing.
function isequal_canonical(
    quad::GenericQuadExpr{CoefType,VarType},
    other::GenericQuadExpr{CoefType,VarType},
) where {CoefType,VarType}
    quad_nozeros = SparseArrays.dropzeros(quad)
    other_nozeros = SparseArrays.dropzeros(other)
    return isequal(quad_nozeros, other_nozeros)
end

# Alias for (Float64, VariableRef)
"""
    QuadExpr

An alias for `GenericQuadExpr{Float64,VariableRef}`, the specific
[`GenericQuadExpr`](@ref) used by JuMP.
"""
const QuadExpr = GenericQuadExpr{Float64,VariableRef}

function Base.convert(
    ::Type{GenericQuadExpr{C,V}},
    v::Union{_Constant,AbstractVariableRef,GenericAffExpr},
) where {C,V}
    return GenericQuadExpr(convert(GenericAffExpr{C,V}, v))
end

function Base.convert(
    ::Type{GenericQuadExpr{C,V}},
    quad::GenericQuadExpr{C,V},
) where {C,V}
    return quad
end

function Base.convert(
    ::Type{GenericQuadExpr{T,V}},
    quad::GenericQuadExpr{S,V},
) where {T,S,V}
    return GenericQuadExpr{T,V}(
        convert(GenericAffExpr{T,V}, quad.aff),
        convert(OrderedDict{UnorderedPair{V},T}, quad.terms),
    )
end

GenericQuadExpr{C,V}() where {C,V} = zero(GenericQuadExpr{C,V})

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
    return MOI.ScalarQuadraticTerm(
        t[2] == t[3] ? 2t[1] : t[1],
        index(t[2]),
        index(t[3]),
    )
end
function MOI.ScalarQuadraticFunction(
    q::GenericQuadExpr{C,GenericVariableRef{T}},
) where {C,T}
    _assert_isfinite(q)
    qterms = MOI.ScalarQuadraticTerm{C}[
        _moi_quadratic_term(t) for t in quad_terms(q)
    ]
    moi_aff = MOI.ScalarAffineFunction(q.aff)
    return MOI.ScalarQuadraticFunction(qterms, moi_aff.terms, moi_aff.constant)
end
function moi_function(aff::GenericQuadExpr)
    return MOI.ScalarQuadraticFunction(aff)
end
function moi_function_type(::Type{<:GenericQuadExpr{T}}) where {T}
    return MOI.ScalarQuadraticFunction{T}
end

function GenericQuadExpr{C,GenericVariableRef{T}}(
    m::GenericModel{T},
    f::MOI.ScalarQuadraticFunction,
) where {C,T}
    quad = GenericQuadExpr{C,GenericVariableRef{T}}(
        GenericAffExpr{C,GenericVariableRef{T}}(
            m,
            MOI.ScalarAffineFunction(f.affine_terms, f.constant),
        ),
    )
    for t in f.quadratic_terms
        v1 = t.variable_1
        v2 = t.variable_2
        coef = t.coefficient
        if v1 == v2
            coef /= 2
        end
        add_to_expression!(
            quad,
            coef,
            GenericVariableRef{T}(m, v1),
            GenericVariableRef{T}(m, v2),
        )
    end
    return quad
end
function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.ScalarQuadraticFunction{C}},
) where {C,T}
    return GenericQuadExpr{C,GenericVariableRef{T}}
end
function jump_function(
    model::GenericModel{T},
    f::MOI.ScalarQuadraticFunction{C},
) where {C,T}
    return GenericQuadExpr{C,GenericVariableRef{T}}(model, f)
end
function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VectorQuadraticFunction{C}},
) where {C,T}
    return Vector{GenericQuadExpr{C,GenericVariableRef{T}}}
end
function jump_function(
    model::GenericModel{T},
    f::MOI.VectorQuadraticFunction{C},
) where {C,T}
    return GenericQuadExpr{C,GenericVariableRef{T}}[
        GenericQuadExpr{C,GenericVariableRef{T}}(model, f) for
        f in MOIU.eachscalar(f)
    ]
end

"""
    _fill_vqf!(terms::Vector{<:MOI.VectorQuadraticTerm}, offset::Int, oi::Int,
               quad::AbstractJuMPScalar)

Fills the vectors terms at indices starting at `offset+1` with the quadratic
terms of `quad`. The output index for all terms is `oi`. Return the index of the
last term added.
"""
function _fill_vqf!(
    terms::Vector{<:MOI.VectorQuadraticTerm},
    offset::Int,
    oi::Int,
    aff::AbstractJuMPScalar,
)
    i = 1
    for term in quad_terms(aff)
        terms[offset+i] =
            MOI.VectorQuadraticTerm(Int64(oi), _moi_quadratic_term(term))
        i += 1
    end
    return offset + length(quad_terms(aff))
end

function MOI.VectorQuadraticFunction(
    quads::Vector{GenericQuadExpr{C,GenericVariableRef{T}}},
) where {C,T}
    num_quadratic_terms = sum(quad -> length(quad_terms(quad)), quads)
    quadratic_terms =
        Vector{MOI.VectorQuadraticTerm{C}}(undef, num_quadratic_terms)
    num_lin_terms = sum(quad -> length(linear_terms(quad)), quads)
    lin_terms = Vector{MOI.VectorAffineTerm{C}}(undef, num_lin_terms)
    constants = Vector{C}(undef, length(quads))
    quad_offset = 0
    lin_offset = 0
    for (i, quad) in enumerate(quads)
        quad_offset = _fill_vqf!(quadratic_terms, quad_offset, i, quad)
        lin_offset = _fill_vaf!(lin_terms, lin_offset, i, quad)
        constants[i] = constant(quad)
    end
    return MOI.VectorQuadraticFunction(quadratic_terms, lin_terms, constants)
end

moi_function(a::Vector{<:GenericQuadExpr}) = MOI.VectorQuadraticFunction(a)

function moi_function_type(::Type{<:Vector{<:GenericQuadExpr{T}}}) where {T}
    return MOI.VectorQuadraticFunction{T}
end

"""
    value(var_value::Function, ex::GenericQuadExpr)

Evaluate `ex` using `var_value(v)` as the value for each variable `v`.
"""
function value(
    var_value::Function,
    ex::GenericQuadExpr{CoefType,VarType},
) where {CoefType,VarType}
    RetType = Base.promote_op(
        (ctype, vtype) -> ctype * var_value(vtype) * var_value(vtype),
        CoefType,
        VarType,
    )
    ret = convert(RetType, value(var_value, ex.aff))
    for (vars, coef) in ex.terms
        ret += coef * var_value(vars.a) * var_value(vars.b)
    end
    return ret
end

"""
    value(v::GenericQuadExpr; result::Int = 1)

Return the value of the `GenericQuadExpr` `v` associated with result index
`result` of the most-recent solution returned by the solver.

Replaces `getvalue` for most use cases.

See also: [`result_count`](@ref).
"""
function value(ex::GenericQuadExpr; result::Int = 1)
    return value(ex) do x
        return value(x; result = result)
    end
end
