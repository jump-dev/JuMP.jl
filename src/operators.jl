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

# Number
# Number--Number obviously already taken care of!
# Number--Variable
(+)(lhs::Number, rhs::Variable) = AffExpr([rhs],[+1.],convert(Float64,lhs))
(-)(lhs::Number, rhs::Variable) = AffExpr([rhs],[-1.],convert(Float64,lhs))
(*)(lhs::Number, rhs::Variable) = AffExpr([rhs],[convert(Float64,lhs)], 0.)
(/)(lhs::Number, rhs::Variable) = error("Cannot divide by variable")
# Number--GenericAffExpr
(+)(lhs::Number, rhs::GenericAffExpr)  = GenericAffExpr(copy(rhs.vars),copy(rhs.coeffs),lhs+rhs.constant)
(-)(lhs::Number, rhs::GenericAffExpr)  = GenericAffExpr(copy(rhs.vars),    -rhs.coeffs ,lhs-rhs.constant)
(*)(lhs::Number, rhs::GenericAffExpr)  = GenericAffExpr(copy(rhs.vars), lhs*rhs.coeffs ,lhs*rhs.constant)
(/)(lhs::Number, rhs::GenericAffExpr)  = error("Cannot divide by an affine expression")
# Number--QuadExpr
(+)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2),copy(rhs.qcoeffs),lhs+rhs.aff)
(-)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2),    -rhs.qcoeffs ,lhs-rhs.aff)
(*)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2), lhs*rhs.qcoeffs ,lhs*rhs.aff)
(/)(lhs::Number, rhs::QuadExpr) = error("Cannot divide by a quadratic expression")


# Variable
# Variable--Number
(+)(lhs::Variable, rhs::Number) = (+)( rhs,lhs)
(-)(lhs::Variable, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::Variable, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::Variable, rhs::Number) = (*)(1./rhs,lhs)
# Variable--Variable
(+)(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.,+1.], 0.)
(-)(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.,-1.], 0.)
(*)(lhs::Variable, rhs::Variable) = QuadExpr([lhs],[rhs],[1.],AffExpr(Variable[],Float64[],0.))
(/)(lhs::Variable, rhs::Variable) = error("Cannot divide a variable by a variable")
# Variable--AffExpr
(+)(lhs::Variable, rhs::AffExpr) = AffExpr(vcat(rhs.vars,lhs),vcat( rhs.coeffs,1.), rhs.constant)
(-)(lhs::Variable, rhs::AffExpr) = AffExpr(vcat(rhs.vars,lhs),vcat(-rhs.coeffs,1.),-rhs.constant)
function (*)(lhs::Variable, rhs::AffExpr)
    n = length(rhs.vars)
    if rhs.constant != 0.      
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),AffExpr([lhs], [rhs.constant], 0.))
    else
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),AffExpr())
    end
end
(/)(lhs::Variable, rhs::AffExpr) = error("Cannot divide a variable by an affine expression")
# Variable--QuadExpr
(+)(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),v+q.aff)
(-)(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),    -q.qcoeffs ,v-q.aff)
(*)(v::Variable, q::QuadExpr) = error("Cannot multiply a variable by a quadratic expression")
(/)(v::Variable, q::QuadExpr) = error("Cannot divide a variable by a quadratic expression")


# AffExpr
# AffExpr--Number
(+)(lhs::GenericAffExpr, rhs::Number) = (+)(+rhs,lhs)
(-)(lhs::GenericAffExpr, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::GenericAffExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::GenericAffExpr, rhs::Number) = (*)(1.0/rhs,lhs)
# AffExpr--Variable
(+)(lhs::AffExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::AffExpr, rhs::Variable) = AffExpr(vcat(lhs.vars,rhs),vcat(+lhs.coeffs,-1.),lhs.constant)
(*)(lhs::AffExpr, rhs::Variable) = (*)(rhs,lhs)
(/)(lhs::AffExpr, rhs::Variable) = error("Cannot divide affine expression by a variable")
# AffExpr--AffExpr
(+)(lhs::GenericAffExpr, rhs::GenericAffExpr) = GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs, rhs.coeffs),lhs.constant+rhs.constant)
(-)(lhs::GenericAffExpr, rhs::GenericAffExpr) = GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,-rhs.coeffs),lhs.constant-rhs.constant)
function (*)(lhs::AffExpr, rhs::AffExpr)
    ret = QuadExpr(Variable[],Variable[],Float64[],AffExpr(Variable[],Float64[],0.))

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
    if lhs.constant != 0 && rhs.constant != 0
        sizehint(ret.aff.vars, n+m)
        sizehint(ret.aff.coeffs, n+m)
    elseif lhs.constant != 0
        sizehint(ret.aff.vars, n)
        sizehint(ret.aff.coeffs, n)
    elseif rhs.constant != 0
        sizehint(ret.aff.vars, m)
        sizehint(ret.aff.coeffs, m)
    end

    # [LHS constant] * RHS
    if lhs.constant != 0
        c = lhs.constant
        for j = 1:m
            push!(ret.aff.vars,   rhs.vars[j])
            push!(ret.aff.coeffs, rhs.coeffs[j] * c)
        end
        ret.aff.constant += c * rhs.constant
    end
    
    # Expr 2 constant * Expr 1 terms
    if rhs.constant != 0
        c = rhs.constant
        for i = 1:m
            push!(ret.aff.vars,   lhs.vars[i])
            push!(ret.aff.coeffs, lhs.coeffs[i] * c)
        end
        # Don't do the following line
        #ret.aff.constant += c * lhs.constant
        # If lhs.constant is 0, its a waste of time
        # If lhs.constant is non-zero, its already done
    end
    
    return ret
end
(/)(lhs::AffExpr, rhs::AffExpr) = error("Cannot divide aff. expression by aff. expression")
# AffExpr--QuadExpr
(+)(a::AffExpr, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),a+q.aff)
(-)(a::AffExpr, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),    -q.qcoeffs ,a-q.aff)
(*)(a::AffExpr, q::QuadExpr) = error("Cannot multiply an aff. expression by a quadratic expression")
(/)(a::AffExpr, q::QuadExpr) = error("Cannot divide an aff. expression by a quadratic expression")


# QuadExpr
# QuadExpr--Number
(+)(lhs::QuadExpr, rhs::Number) = (+)(+rhs,lhs)
(-)(lhs::QuadExpr, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::QuadExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::QuadExpr, rhs::Number) = (*)(1.0/rhs,lhs)
# QuadExpr--Variable
(+)(q::QuadExpr, v::Variable) = (+)(v,q)
(-)(q::QuadExpr, v::Variable) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff-v)
(*)(q::QuadExpr, v::Variable) = error("Cannot multiply a quadratic expression by a variable")
(/)(q::QuadExpr, v::Variable) = error("Cannot divide a quadratic expression by a variable")
# QuadExpr--AffExpr
(+)(q::QuadExpr, a::AffExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff+a)
(-)(q::QuadExpr, a::AffExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff-a)
(*)(q::QuadExpr, a::AffExpr) = error("Cannot multiply a quadratic expression by an aff. expression")
(/)(q::QuadExpr, a::AffExpr) = error("Cannot divide a quadratic expression by an aff. expression")
# QuadExpr--QuadExpr
(+)(q1::QuadExpr, q2::QuadExpr) = QuadExpr(vcat(q1.qvars1,   q2.qvars1),
                                                                                     vcat(q1.qvars2,   q2.qvars2),
                                                                                     vcat(q1.qcoeffs,  q2.qcoeffs),
                                                                                     q1.aff + q2.aff)
(-)(q1::QuadExpr, q2::QuadExpr) = QuadExpr(vcat(q1.qvars1,   q2.qvars1),
                                                                                     vcat(q1.qvars2,   q2.qvars2),
                                                                                     vcat(q1.qcoeffs, -q2.qcoeffs),
                                                                                     q1.aff - q2.aff)
(*)(q1::QuadExpr, q2::QuadExpr) = error("Cannot multiply two quadratic expressions")
(/)(q1::QuadExpr, q2::QuadExpr) = error("Cannot divide a quadratic expression by a quadratic expression")

# LinearConstraint
# Number--???
(<=)(lhs::Number, rhs::Variable) = (>=)(rhs, lhs)
(==)(lhs::Number, rhs::Variable) = (==)(rhs, lhs)
(>=)(lhs::Number, rhs::Variable) = (<=)(rhs, lhs)

(<=)(lhs::Number, rhs::AffExpr) = (>=)(rhs, lhs)
(==)(lhs::Number, rhs::AffExpr) = (==)(rhs, lhs)
(>=)(lhs::Number, rhs::AffExpr) = (<=)(rhs, lhs)
(<=)(lhs::Number, rhs::QuadExpr) = (>=)(rhs, lhs)
(==)(lhs::Number, rhs::QuadExpr) = (==)(rhs, lhs)
(>=)(lhs::Number, rhs::QuadExpr) = (<=)(rhs, lhs)
# Variable--???
(<=)(lhs::Variable, rhs::Number) = (<=)(lhs - rhs, 0.0)
(==)(lhs::Variable, rhs::Number) = (==)(lhs - rhs, 0.0)
(>=)(lhs::Variable, rhs::Number) = (>=)(lhs - rhs, 0.0)

(<=)(lhs::Variable, rhs::Variable) = (<=)(lhs - rhs, 0.0)
(==)(lhs::Variable, rhs::Variable) = (==)(lhs - rhs, 0.0)
(>=)(lhs::Variable, rhs::Variable) = (>=)(lhs - rhs, 0.0)

(<=)(lhs::Variable, rhs::AffExpr) = (<=)(lhs - rhs, 0.0)
(==)(lhs::Variable, rhs::AffExpr) = (==)(lhs - rhs, 0.0)
(>=)(lhs::Variable, rhs::AffExpr) = (>=)(lhs - rhs, 0.0)
(<=)(lhs::Variable, rhs::QuadExpr) = (<=)(lhs - rhs, 0.0)
(==)(lhs::Variable, rhs::QuadExpr) = (==)(lhs - rhs, 0.0)
(>=)(lhs::Variable, rhs::QuadExpr) = (>=)(lhs - rhs, 0.0)
# AffExpr--???
(<=)(lhs::AffExpr, rhs::Number) = LinearConstraint(lhs,            -Inf,rhs-lhs.constant)
(==)(lhs::AffExpr, rhs::Number) = LinearConstraint(lhs,rhs-lhs.constant,rhs-lhs.constant)
(>=)(lhs::AffExpr, rhs::Number) = LinearConstraint(lhs,rhs-lhs.constant,             Inf)
(<=)(lhs::AffExpr, rhs::Variable) = (<=)(lhs-rhs, 0.0)
(==)(lhs::AffExpr, rhs::Variable) = (==)(lhs-rhs, 0.0)
(>=)(lhs::AffExpr, rhs::Variable) = (>=)(lhs-rhs, 0.0)
(<=)(lhs::AffExpr, rhs::AffExpr) = (<=)(lhs-rhs, 0.0)
(==)(lhs::AffExpr, rhs::AffExpr) = (==)(lhs-rhs, 0.0)
(>=)(lhs::AffExpr, rhs::AffExpr) = (>=)(lhs-rhs, 0.0)
(<=) (lhs::AffExpr, rhs::QuadExpr) = (<=)(lhs-rhs, 0)
(==) (lhs::AffExpr, rhs::QuadExpr) = (==)(lhs-rhs, 0)
(>=) (lhs::AffExpr, rhs::QuadExpr) = (>=)(lhs-rhs, 0)

# There's no easy way to allow operator overloads for range constraints.
# Use macros instead.

# QuadConstraint
# QuadConstraint--Number
(<=) (lhs::QuadExpr, rhs::Number)   = QuadConstraint( QuadExpr(copy(lhs.qvars1), copy(lhs.qvars2), lhs.qcoeffs,lhs.aff - rhs), :<=   )
(==) (lhs::QuadExpr, rhs::Number)   = QuadConstraint( QuadExpr(copy(lhs.qvars1), copy(lhs.qvars2), lhs.qcoeffs,lhs.aff - rhs), :(==) )
(>=) (lhs::QuadExpr, rhs::Number)   = QuadConstraint( QuadExpr(copy(lhs.qvars1), copy(lhs.qvars2), lhs.qcoeffs,lhs.aff - rhs), :>=   )
(<=) (lhs::QuadExpr, rhs::Variable) = (<=)(lhs-rhs, 0)
(==) (lhs::QuadExpr, rhs::Variable) = (==)(lhs-rhs, 0)
(>=) (lhs::QuadExpr, rhs::Variable) = (>=)(lhs-rhs, 0)
(<=) (lhs::QuadExpr, rhs::AffExpr)  = (<=)(lhs-rhs, 0)
(==) (lhs::QuadExpr, rhs::AffExpr)  = (==)(lhs-rhs, 0)
(>=) (lhs::QuadExpr, rhs::AffExpr)  = (>=)(lhs-rhs, 0)
(<=) (lhs::QuadExpr, rhs::QuadExpr) = (<=)(lhs-rhs, 0)
(==) (lhs::QuadExpr, rhs::QuadExpr) = (==)(lhs-rhs, 0)
(>=) (lhs::QuadExpr, rhs::QuadExpr) = (>=)(lhs-rhs, 0)

#Vectorization Stuff
function dot{T<:Real}(lhs::JuMPDict{Variable}, rhs::Vector{T})
    @assert length(lhs.indexsets) == 1
    @assert length(lhs.indexsets[1]) == length(rhs)
    return AffExpr(lhs.innerArray, rhs, 0.0)
end
dot{T<:Real}(lhs::Vector{T}, rhs::JuMPDict{Variable}) = dot(rhs,lhs)

function bigdot{T<:Real}(lhs::Array{T,2},rhs::JuMPDict{Variable})
    matsize = size(lhs)
    @assert length(matsize) == length(rhs.indexsets)

    coeffs = Float64[]; sizehint(coeffs, matsize[1]*matsize[2])
    vars  = Variable[]; sizehint(vars,   matsize[1]*matsize[2])

    for i = 1:matsize[1]
        for j = 1:matsize[2]
            push!(coeffs, lhs[i,j])
            push!(vars,   rhs.innerArray[i,j])
        end
    end

    return AffExpr(vars, coeffs, 0.0)
end
bigdot{T<:Real}(lhs::JuMPDict{Variable},rhs::Array{T,2}) = bigdot(rhs,lhs)

function bigdot{T<:Real}(lhs::Array{T,3},rhs::JuMPDict{Variable})
    matsize = size(lhs)
    @assert length(matsize) == length(rhs.indexsets)

    coeffs = Float64[]; sizehint(coeffs, matsize[1]*matsize[2]*matsize[3])
    vars  = Variable[]; sizehint(vars,   matsize[1]*matsize[2]*matsize[3])

    for i = 1:matsize[1]
        for j = 1:matsize[2]
            for k = 1:matsize[3]
                push!(coeffs, lhs[i,j,k])
                push!(vars,   rhs.innerArray[i,j,k])
            end
        end
    end

    return AffExpr(vars, coeffs, 0.0)
end
bigdot{T<:Real}(lhs::JuMPDict{Variable},rhs::Array{T,3}) = bigdot(rhs,lhs)