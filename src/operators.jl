#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

# Overloads
#
# Different objects that must all interact:
# 1. Number
# 2. Variable
# 3. AffExpr
# 4. QuadExpr

#########
# Generic
#########
Base.promote_rule{T,S}(::Type{T},::Type{GenericAffExpr {T,S}}) = GenericAffExpr {T,S}
Base.promote_rule{T,S}(::Type{T},::Type{GenericQuadExpr{T,S}}) = GenericQuadExpr{T,S}
Base.promote_rule{T,S}(::Type{S},::Type{GenericAffExpr {T,S}}) = GenericAffExpr {T,S}
Base.promote_rule{T,S}(::Type{S},::Type{GenericQuadExpr{T,S}}) = GenericQuadExpr{T,S}
Base.convert{T,S}(::Type{GenericAffExpr{T,S}}, x::T) = GenericAffExpr(S[],  T[], x)
Base.convert{T,S}(::Type{GenericAffExpr{T,S}}, x::S) = GenericAffExpr(S[x], T[one(T)], zero(T))
Base.convert{T,S}(::Type{GenericQuadExpr{T,S}}, x::T) = GenericQuadExpr(S[], S[], T[], convert(GenericAffExpr, x))
Base.convert{T,S}(::Type{GenericQuadExpr{T,S}}, x::S) = GenericQuadExpr(S[], S[], T[], convert(GenericAffExpr, x))
Base.convert{T,S}(::Type{GenericQuadExpr{T,S}}, x::GenericAffExpr{T,S}) = GenericQuadExpr(S[], S[], T[], x)
jtypes = [:GenericAffExpr,:GenericQuadExpr]
for type1 in jtypes
    @eval begin
        (+)(x::$(type1)) = x
        (-)(x::$(type1)) = zero(x) - x
    end
    for op in (:+, :-, :*, :/)
        @eval begin
            $(op){T,V}(lhs::T, rhs::$type1{T,V}) = $(op)(promote(lhs,rhs)...)
            # $(op){T,V}(lhs::$type1{T,V}, rhs::T) = $(op)(promote(lhs,rhs)...)
        end
        for type2 in jtypes
            @eval begin
                $(op){T,V}(lhs::$(type1){T,V}, rhs::$(type2){T,V}) = $(op)(promote(lhs,rhs)...)
                $(op){T,V}(lhs::$(type2){T,V}, rhs::$(type1){T,V}) = $(op)(promote(lhs,rhs)...)
            end
        end
    end
end

# AffExpr--AffExpr
(+){T<:GenericAffExpr}(lhs::T, rhs::T) = T(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs, rhs.coeffs),lhs.constant+rhs.constant)
(-){T<:GenericAffExpr}(lhs::T, rhs::T) = T(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,-rhs.coeffs),lhs.constant-rhs.constant)
function (*){T<:GenericAffExpr}(lhs::T, rhs::T)
    ret = zero(T)
    zero_c = zero(ret.constant)

    # Quadratic terms
    n = length(lhs.coeffs)
    m = length(rhs.coeffs)
    sizehint(ret.qvars1, n*m)
    sizehint(ret.qvars2, n*m)
    sizehint(ret.qcoeffs, n*m)
    for i = 1:n
        for j = 1:m
            push!(ret.qvars1,  lhs.vars[i])
            push!(ret.qvars2,  rhs.vars[j])
            push!(ret.qcoeffs, lhs.coeffs[i]*rhs.coeffs[j])
        end
    end
    
    # Try to preallocate space for aff
    if lhs.constant != zero_c && rhs.constant != zero_c
        sizehint(ret.aff.vars, n+m)
        sizehint(ret.aff.coeffs, n+m)
    elseif lhs.constant != zero_c
        sizehint(ret.aff.vars, n)
        sizehint(ret.aff.coeffs, n)
    elseif rhs.constant != zero_c
        sizehint(ret.aff.vars, m)
        sizehint(ret.aff.coeffs, m)
    end

    # [LHS constant] * RHS
    if lhs.constant != zero_c
        c = lhs.constant
        for j = 1:m
            push!(ret.aff.vars,   rhs.vars[j])
            push!(ret.aff.coeffs, rhs.coeffs[j] * c)
        end
        ret.aff.constant += c * rhs.constant
    end
    
    # Expr 2 constant * Expr 1 terms
    if rhs.constant != zero_c
        c = rhs.constant
        for i = 1:m
            push!(ret.aff.vars,   lhs.vars[i])
            push!(ret.aff.coeffs, lhs.coeffs[i] * c)
        end
    end
    
    return ret
end
(/){T,S}(lhs::GenericAffExpr{T,S}, rhs::GenericAffExpr{T,S}) = 
    error("Cannot divide aff. expression by aff. expression")

# QuadExpr--QuadExpr
(+){T<:GenericQuadExpr}(q1::T, q2::T) = T(vcat(q1.qvars1,  q2.qvars1),  vcat(q1.qvars2, q2.qvars2), 
                                          vcat(q1.qcoeffs, q2.qcoeffs), q1.aff + q2.aff)
(-){T<:GenericQuadExpr}(q1::T, q2::T) = T(vcat(q1.qvars1,   q2.qvars1),   vcat(q1.qvars2, q2.qvars2),
                                          vcat(q1.qcoeffs, -q2.qcoeffs),  q1.aff - q2.aff)
(*){T<:GenericQuadExpr}(q1::T, q2::T) = error("Cannot multiply two quadratic expressions")
(/){T<:GenericQuadExpr}(q1::T, q2::T) = error("Cannot divide a quadratic expression by a quadratic expression")


###################
# Variable-specific
###################
Base.promote_rule{T<:Number}(::Type{T},::Type{Variable}) = AffExpr
# Variable--Variable
(+)(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.0,+1.0], 0.0)
(-)(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.0,-1.0], 0.0)
(*)(lhs::Variable, rhs::Variable) = QuadExpr([lhs],[rhs],[1.0],AffExpr(Variable[],Float64[],0.0))
(/)(lhs::Variable, rhs::Variable) = error("Cannot divide a variable by a variable")

# LinearConstraint
# Number--???
for (sgn, osgn) in ( (:<=,:>=), (:(==),:(==)), (:>=,:<=) )
    for typ in (:Variable, :AffExpr, :QuadExpr)
        @eval $(sgn)(lhs::Number, rhs::$(typ)) = $(osgn)(rhs, lhs)
    end
    # Variable--???
    for typ in (:Number, :Variable, :AffExpr, :QuadExpr)
        @eval $(sgn)(lhs::Variable, rhs::$(typ)) = $(sgn)(lhs-rhs, 0.0)
    end
end
# AffExpr--???
(<=)(lhs::AffExpr, rhs::Number) = LinearConstraint(lhs,            -Inf,rhs-lhs.constant)
(==)(lhs::AffExpr, rhs::Number) = LinearConstraint(lhs,rhs-lhs.constant,rhs-lhs.constant)
(>=)(lhs::AffExpr, rhs::Number) = LinearConstraint(lhs,rhs-lhs.constant,             Inf)
for sgn in (:<=, :(==), :>=)
    for typ in (:Variable, :AffExpr, :QuadExpr)
        @eval $(sgn)(lhs::AffExpr, rhs::$(typ)) = $(sgn)(lhs-rhs, 0.0)
    end
end
# There's no easy way to allow operator overloads for range constraints.
# Use macros instead.

# QuadConstraint
# QuadConstraint--Number
for sgn in (:<=, :(==), :>=)
    @eval $(sgn)(lhs::QuadExpr, rhs::Number) = 
        QuadConstraint( QuadExpr(copy(lhs.qvars1), copy(lhs.qvars2), lhs.qcoeffs,lhs.aff - rhs), $(quot(sgn)))
    for typ in (:Variable, :AffExpr, :QuadExpr)
        @eval $(sgn)(lhs::QuadExpr, rhs::$(typ)) = $(sgn)(lhs-rhs, 0)
    end
end

#############################################################################
# High-level operators
# Currently supported
#  - sum
#  - dot
#############################################################################

Base.sum{T<:Real}(j::JuMPDict{T}) = sum(j.innerArray)
Base.sum(j::JuMPDict{Variable}) = AffExpr(vec(j.innerArray), ones(length(j.innerArray)), 0.0)
Base.sum(j::Array{Variable}) = AffExpr(vec(j), ones(length(j)), 0.0)
function Base.sum{S,T}(affs::Array{GenericAffExpr{S,T}})
    new_aff = GenericAffExpr{S,T}()
    for aff in affs
        append!(new_aff, aff)
    end
    return new_aff
end

function Base.dot{T,S}(lhs::Array{T}, rhs::JuMPDict{S})
    size(lhs) == size(rhs.innerArray) || error("Incompatible dimensions")
    dot(lhs,rhs.innerArray)
end
Base.dot{S,T}(lhs::JuMPDict{S},rhs::Array{T}) = dot(rhs,lhs)

Base.dot{T<:Real}(lhs::Array{T}, rhs::JuMPDict{Float64}) = dot(vec(lhs), vec(rhs.innerArray))
Base.dot{T<:Real}(lhs::JuMPDict{Float64}, rhs::Array{T}) = dot(vec(rhs), vec(lhs.innerArray))
Base.dot{T<:Real}(lhs::Array{T}, rhs::Array{Variable})   = AffExpr(vec(rhs), vec(float(lhs)), 0.0)
Base.dot{T<:Real}(rhs::Array{Variable}, lhs::Array{T})   = AffExpr(vec(rhs), vec(float(lhs)), 0.0)

function Base.dot(lhs::JuMPDict{Variable},rhs::JuMPDict{Variable})
    size(lhs.innerArray) == size(rhs.innerArray) || error("Incompatible dimensions") 
    return QuadExpr(vec(lhs.innerArray), vec(rhs.innerArray), ones(length(lhs.innerArray)), AffExpr())
end

function Base.dot(lhs::JuMPDict{Float64},rhs::JuMPDict{Float64})
    size(lhs.innerArray) == size(rhs.innerArray) || error("Incompatible dimensions") 
    return sum(lhs.innerArray .* rhs.innerArray)
end

#############################################################################
# JuMPDict comparison operators (all errors)

for sgn in (:<=, :(==), :>=)
    for term in (:Real, :Variable, :AffExpr)
        @eval $(sgn)(a::JuMPDict, b::$(term)) = error("Cannot construct constraint with a JuMPDict term")
        @eval $(sgn)(a::$(term), b::JuMPDict) = error("Cannot construct constraint with a JuMPDict term")
    end
end
