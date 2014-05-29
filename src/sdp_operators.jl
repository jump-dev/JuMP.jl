# DualExpr
(-)(x::DualExpr) = DualExpr(copy(x.vars), map(-,x.coeffs), -x.constant)
# Number--DualExpr
(*)(lhs::Number, rhs::DualExpr) = DualExpr(copy(rhs.vars), lhs*rhs.coeffs, lhs*rhs.constant)
# Number--MatrixVar
(*)(lhs::Number,rhs::MatrixVar) = MatrixExpr(concat(rhs), concat(lhs*𝕀), concat(𝕀), spzeros(size(rhs)...))
# Number--SDPVar
(*)(lhs::Number, rhs::SDPVar) = MatrixExpr(concat(rhs), concat(lhs*𝕀), concat(𝕀), spzeros(size(rhs)...))
# Number--MatrixExpr
# function (+)(lhs::Number, rhs::MatrixExpr)
#     size(rhs) == (1,1) || error("Cannot add scalar to matrix expression")
#     return lhs + rhs[1]
# end
# function (-)(lhs::Number, rhs::MatrixExpr)
#     size(rhs) == (1,1) || error("Cannot add scalar to matrix expression")
#     return lhs - rhs[1]
# end
(*)(lhs::Number, rhs::MatrixExpr)  = MatrixExpr(copy(rhs.elem), lhs*rhs.pre, copy(rhs.post), lhs*rhs.constant)
# Number--MatrixFuncVar
(+)(lhs::Number, rhs::MatrixFuncVar) = ScalarExpr(concat( rhs),AffExpr(lhs),NormExpr())
(-)(lhs::Number, rhs::MatrixFuncVar) = ScalarExpr(concat(-rhs),AffExpr(lhs),NormExpr())
(*)(lhs::Number, rhs::MatrixFuncVar) = MatrixFuncVar(lhs*rhs.expr,rhs.func)
# Number--NormExpr
(+)(lhs::Number,rhs::NormExpr)  = ScalarExpr(MatrixFuncVar[], AffExpr(lhs),  rhs)
(-)(lhs::Number,rhs::NormExpr)  = ScalarExpr(MatrixFuncVar[], AffExpr(lhs), -rhs)
(*)(lhs::Number,rhs::NormExpr)  =
    ( lhs == 0 ? NormExpr() : NormExpr(copy(rhs.vars), rhs.coeffs*lhs, rhs.form) )
# Number--ScalarExpr
(+)(lhs::Number, rhs::ScalarExpr) = ScalarExpr(        copy(rhs.matvars), lhs+rhs.aff, copy(rhs.normexpr))
(-)(lhs::Number, rhs::ScalarExpr) = ScalarExpr(       map(-,rhs.matvars), lhs-rhs.aff,     -rhs.normexpr )
(*)(lhs::Number, rhs::ScalarExpr) = ScalarExpr(map(x->lhs*x,rhs.matvars), lhs*rhs.aff,  lhs*rhs.normexpr )

# Variable--AbstractArray{T,2}
function (*){T<:Number}(lhs::Variable, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix coefficient must be symmetric")
    DualExpr(concat(lhs), concat(rhs), spzeros(size(rhs)...))
end
# Variable--MatrixExpr
function (+)(lhs::Variable, rhs::MatrixExpr)
    size(rhs) == (1,1) || error("Cannot add scalar to matrix expression")
    return lhs + rhs[1]
end
function (-)(lhs::Variable, rhs::MatrixExpr)
    size(rhs) == (1,1) || error("Cannot add scalar to matrix expression")
    return lhs - rhs[1]
end
# Variable--MatrixFuncVar
(+)(lhs::Variable, rhs::MatrixFuncVar) = ScalarExpr(concat( rhs), convert(AffExpr,lhs), NormExpr())
(-)(lhs::Variable, rhs::MatrixFuncVar) = ScalarExpr(concat(-rhs), convert(AffExpr,lhs), NormExpr())
# Variable--NormExpr
(+)(lhs::Variable,rhs::NormExpr)  = ScalarExpr(MatrixFuncVar[], convert(AffExpr,lhs),  rhs)
(-)(lhs::Variable,rhs::NormExpr)  = ScalarExpr(MatrixFuncVar[], convert(AffExpr,lhs), -rhs)
# Variable--ScalarExpr
(+)(lhs::Variable, rhs::ScalarExpr) = ScalarExpr(copy(rhs.matvars), lhs+rhs.aff, copy(rhs.normexpr))
(-)(lhs::Variable, rhs::ScalarExpr) = ScalarExpr(    -rhs.matvars , lhs-rhs.aff,     -rhs.normexpr )

# AffExpr--AbstractArray{T,2}
function (*){T<:Number}(lhs::AffExpr, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix coefficient must be symmetric")
    DualExpr(copy(lhs.vars), map(x->x*rhs,lhs.coeffs), lhs.constant*rhs)
end
# AffExpr--MatrixExpr
function (+)(lhs::AffExpr, rhs::MatrixExpr)
    size(rhs) == (1,1) || error("Cannot add scalar to matrix expression")
    return lhs + rhs[1]
end
function (-)(lhs::AffExpr, rhs::MatrixExpr)
    size(rhs) == (1,1) || error("Cannot add scalar to matrix expression")
    return lhs - rhs[1]
end
# AffExpr--MatrixFuncVar
(+)(lhs::AffExpr, rhs::MatrixFuncVar) = ScalarExpr(concat( rhs), copy(lhs), NormExpr())
(-)(lhs::AffExpr, rhs::MatrixFuncVar) = ScalarExpr(concat(-rhs), copy(lhs), NormExpr())
# AffExpr--NormExpr
(+)(lhs::AffExpr,rhs::NormExpr) = ScalarExpr(MatrixFuncVar[], lhs,  rhs)
(-)(lhs::AffExpr,rhs::NormExpr) = ScalarExpr(MatrixFuncVar[], lhs, -rhs)
# AffExpr--ScalarExpr
(+)(lhs::AffExpr, rhs::ScalarExpr) = ScalarExpr(copy(rhs.matvars), lhs+rhs.aff, copy(rhs.normexpr))
(-)(lhs::AffExpr, rhs::ScalarExpr) = ScalarExpr(    -rhs.matvars , lhs-rhs.aff,     -rhs.normexpr )

# QuadExpr--SDPVar
# QuadExpr--MatrixExpr

# DualExpr
# DualExpr--Number
(*)(lhs::DualExpr, rhs::Number) = DualExpr(copy(lhs.vars), lhs.coeffs*rhs, lhs.constant*rhs)
(/)(lhs::DualExpr, rhs::Number) = DualExpr(copy(lhs.vars), lhs.coeffs/rhs, lhs.constant/rhs)
# DualExpr--AbstractArray{T,2}
function (+){T<:Number}(lhs::DualExpr, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix must be symmetric")
    DualExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant+rhs)
end
function (-){T<:Number}(lhs::DualExpr, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix must be symmetric")
    DualExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant-rhs)
end
function (*){T<:Number}(lhs::DualExpr, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix coefficient must be symmetric")
    coeffs = copy(lhs.coeffs)
    for it in 1:length(coeffs)
        coeffs[it] *= rhs
    end
    DualExpr(copy(lhs.vars), coeffs, lhs.constant*rhs)
end
# DualExpr--DualExpr
(+)(lhs::DualExpr, rhs::DualExpr) = DualExpr(concat(lhs.vars,rhs.vars), concat(lhs.coeffs, rhs.coeffs), lhs.constant+rhs.constant)
(-)(lhs::DualExpr, rhs::DualExpr) = DualExpr(concat(lhs.vars,rhs.vars), concat(lhs.coeffs,-rhs.coeffs), lhs.constant-rhs.constant)
# AbstractArray{T,2}
# AbstractArray{T,2}--Variable
function (*){T<:Number}(lhs::AbstractArray{T,2}, rhs::Variable)
    issym(lhs) || error("Matrix coefficient must be symmetric")
    DualExpr(concat(rhs), concat(lhs), spzeros(size(lhs)...))
end
# AbstractArray{T,2}--AffExpr
function (*){T<:Number}(lhs::AbstractArray{T,2}, rhs::AffExpr)
    issym(lhs) || error("Matrix coefficient must be symmetric")
    DualExpr(copy(rhs.vars), map(x->lhs*x,rhs.coeffs), lhs*rhs.constant)
end
# AbstractArray{T,2}--DualExpr
function (+){T<:Number}(lhs::AbstractArray{T,2}, rhs::DualExpr)
    issym(rhs) || error("Matrix must be symmetric")
    DualExpr(copy(rhs.vars), copy(rhs.coeffs), lhs+rhs.constant)
end
function (-){T<:Number}(lhs::AbstractArray{T,2}, rhs::DualExpr)
    issym(rhs) || error("Matrix must be symmetric")
    DualExpr(copy(rhs.vars), -rhs.coeffs, lhs-rhs.constant)
end
function (*){T<:Number}(lhs::AbstractArray{T,2}, rhs::DualExpr)
    issym(rhs) || error("Matrix coefficient must be symmetric")
    coeffs = copy(rhs.coeffs)
    for it in 1:length(coeffs)
        coeffs[it] = lhs * coeffs[it]
    end
    DualExpr(copy(rhs.vars), coeffs, lhs*rhs.constant)
end
# AbstractArray{T}--MatrixVar
function (+){T<:Number}(lhs::AbstractArray{T},rhs::MatrixVar)
    ( (size(lhs,1) == size(rhs,1)) && (size(lhs,2) == size(rhs,2)) ) || error("Cannot add matrices of unequal size")
    MatrixExpr(concat(rhs), concat(𝕀), concat(𝕀), lhs)
end
function (-){T<:Number}(lhs::AbstractArray{T},rhs::MatrixVar)
    ( (size(lhs,1) == size(rhs,1)) && (size(lhs,2) == size(rhs,2)) ) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(rhs), concat(-𝕀), concat(I), lhs)
end
function (*){T<:Number}(lhs::AbstractArray{T},rhs::MatrixVar)
    (size(lhs,2) == size(rhs,1)) || error("Matrix multiplication with incompatible sizes")
    MatrixExpr(concat(rhs), concat(lhs), concat(𝕀), spzeros(size(lhs,1),size(rhs,2)))
end
# AbstractArray{T,2}--SDPVar
function (+){T<:Number}(lhs::AbstractArray{T,2}, rhs::SDPVar)
    (size(lhs) == size(rhs)) || error("Cannot add matrices of unequal size")
    MatrixExpr(concat(rhs), concat(𝕀), concat(𝕀), lhs)
end
function (-){T<:Number}(lhs::AbstractArray{T,2}, rhs::SDPVar)
    (size(lhs) == size(rhs)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(rhs), concat(-𝕀), concat(𝕀), lhs)
end
function (*){T<:Number}(lhs::AbstractArray{T}, rhs::SDPVar)
    size(lhs,2) == size(rhs,1) || error("Incompatible dimensions")
    MatrixExpr(concat(rhs), concat(lhs), concat(𝕀), spzeros(size(lhs,1),size(rhs,2)))
end
# AbstractArray{T}--MatrixExpr
function (+){T<:Number}(lhs::AbstractArray{T}, rhs::MatrixExpr)
    (size(lhs,1) == size(rhs,1) && size(lhs,2) == size(rhs,2)) || error("Cannot add matrices of unequal size")
    MatrixExpr(copy(rhs.elem), copy(rhs.pre), copy(rhs.post), lhs+rhs.constant)
end
function (-){T<:Number}(lhs::AbstractArray{T}, rhs::MatrixExpr)
    (size(lhs,1) == size(rhs,1) && size(lhs,2) == size(rhs,2)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(copy(rhs.elem), -rhs.pre, copy(rhs.post), lhs-rhs.constant)
end
function (*){T<:Number}(lhs::AbstractArray{T}, rhs::MatrixExpr)
    size(lhs,2) == size(rhs,1) || error("Incompatible dimensions")
    MatrixExpr(copy(rhs.elem), map(x->lhs*x, rhs.pre), copy(rhs.post), lhs*rhs.constant)
end
# AbstractArray{T,2}--MatrixFuncVar
function (*){T<:Number}(lhs::AbstractArray{T,2},rhs::MatrixFuncVar)
    issym(lhs) || error("Matrix must be symmetric")
    DualExpr(concat(rhs), concat(lhs), spzeros(size(lhs)...))
end
# AbstractArray{T,2}--ScalarExpr
function (*){T<:Number}(lhs::AbstractArray{T,2},rhs::ScalarExpr)
    isequal(rhs.normexpr,NormExpr()) || error("Cannot multiply norm expression with a matrix")
    issym(lhs) || error("Matrix must by symmetric")
    DualExpr(concat(rhs.matvars,rhs.aff.vars), concat([lhs for c in rhs.matvars],[lhs*c for c in rhs.aff.coeffs]), lhs*rhs.aff.constant)
end
# # UniformScaling
# # UniformScaling--SDPVar
# (+)(lhs::UniformScaling,rhs::SDPVar) = MatrixExpr(concat(rhs),concat( 𝕀) ,concat(𝕀),lhs)
# (-)(lhs::UniformScaling,rhs::SDPVar) = MatrixExpr(concat(rhs),concat(-𝕀) ,concat(𝕀),lhs)
# (*)(lhs::UniformScaling,rhs::SDPVar) = MatrixExpr(concat(rhs),concat(lhs),concat(𝕀),spzeros(size(rhs)...))
# # UniformScaling--MatrixExpr
# (+)(lhs::UniformScaling,rhs::MatrixExpr) = MatrixExpr(copy(rhs.elem),copy(rhs.pre),copy(rhs.post),rhs.constant+lhs)
# (-)(lhs::UniformScaling,rhs::MatrixExpr) = MatrixExpr(copy(rhs.elem),-rhs.pre,copy(rhs.post),lhs-rhs.constant)
# (*)(lhs::UniformScaling,rhs::MatrixExpr) = MatrixExpr(copy(rhs.elem),lhs.λ*rhs.pre,copy(rhs.post),lhs.λ*rhs.constant)

# MatrixVar
(-)(var::MatrixVar) = MatrixExpr(concat(var),concat(-𝕀),concat(𝕀),spzeros(size(var)...))
# MatrixVar--Number
(*)(lhs::MatrixVar, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::MatrixVar, rhs::Number) = (*)(1/rhs,lhs)
# MatrixVar--AbstractArray{T}
function (+){T<:Number}(lhs::MatrixVar, rhs::AbstractArray{T})
    ( (size(lhs,1) == size(rhs,1)) && (size(lhs,2) == size(rhs,2)) ) || error("Cannot add matrices of unequal size")
    MatrixExpr(concat(lhs), concat(𝕀), concat(𝕀), rhs)
end
function (-){T<:Number}(lhs::MatrixVar, rhs::AbstractArray{T})
    ( (size(lhs,1) == size(rhs,1)) && (size(lhs,2) == size(rhs,2)) ) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(lhs), concat(𝕀), concat(𝕀), -rhs)
end
function (*){T<:Number}(lhs::MatrixVar, rhs::AbstractArray{T})
    (size(lhs,2) == size(rhs,1)) || error("Matrix multiplication with incompatible sizes")
    MatrixExpr(concat(lhs), concat(𝕀), concat(rhs), spzeros(size(lhs,1),size(rhs,2)))
end
# MatrixVar--MatrixVar
function (+)(lhs::MatrixVar,rhs::MatrixVar)
    size(lhs) == size(rhs) || error("Cannot add matrices of unequal size")
    MatrixExpr(concat(lhs,rhs), concat(𝕀,𝕀), concat(𝕀,𝕀), spzeros(size(lhs)...))
end
function (-)(lhs::MatrixVar,rhs::MatrixVar)
    size(lhs) == size(rhs) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(lhs,rhs), concat(𝕀,-𝕀), concat(𝕀,𝕀), spzeros(size(lhs)...))
end
function (*)(lhs::MatrixVar,rhs::MatrixVar)
    (m,p),(q,n) = size(lhs), size(rhs)
    p == q || error("Imcompatible dimensions for matrix multiplication")
    arr = Array(QuadExpr, m, n)
    for i in 1:m, j in 1:n
        arr[i,j] = QuadExpr(vec(lhs[i,:]),rhs[:,j],ones(p), AffExpr())
    end
    bl = BlockMatrix(arr)
    MatrixExpr(concat(bl), concat(𝕀), concat(𝕀), spzeros(size(lhs,1),size(rhs,2)))
end
# MatrixVar--SDPVar
function (+)(lhs::MatrixVar,rhs::SDPVar)
    size(lhs) == size(rhs) || error("Cannot add matrices of unequal size")
    MatrixExpr(concat(lhs,rhs), concat(𝕀,𝕀), concat(𝕀,𝕀), spzeros(size(lhs)...))
end
function (-)(lhs::MatrixVar,rhs::SDPVar)
    size(lhs) == size(rhs) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(lhs,rhs), concat(𝕀,-𝕀), concat(𝕀,𝕀), spzeros(size(lhs)...))
end
# MatrixVar--MatrixExpr
function (+)(lhs::MatrixVar,rhs::MatrixExpr)
    size(lhs) == size(rhs) || error("Cannot add matrices of unequal size")
    MatrixExpr(concat(lhs,rhs.elem), concat(𝕀,rhs.pre), concat(𝕀,rhs.post), copy(rhs.constant))
end
function (-)(lhs::MatrixVar,rhs::MatrixExpr)
    size(lhs) == size(rhs) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(lhs,rhs.elem), concat(𝕀,-rhs.pre), concat(𝕀,rhs.post), -rhs.constant)
end
function (*)(lhs::MatrixVar,rhs::MatrixExpr)
    (m,p),(q,n) = size(lhs), size(rhs)
    p == q || error("Imcompatible dimensions for matrix multiplication")
    m == n == 1 || error("Only quadratic forms are allowed")
    ret = QuadExpr()
    for k in 1:p
        ret += lhs[k]*rhs[k]
    end
    return ret
end

# SDPVar
(-)(var::SDPVar) = MatrixExpr(concat(var),concat(-𝕀),concat(𝕀),spzeros(size(var)...))
# SDPVar--Number
(*)(lhs::SDPVar, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::SDPVar, rhs::Number) = (*)(1/rhs,lhs)
# SDPVar--AbstractArray{T}
function (+){T<:Number}(lhs::SDPVar, rhs::AbstractArray{T})
    (size(lhs,1) == size(rhs,1) && size(lhs,2) == size(lhs,2)) || error("Cannot add matrices of unequal size")
    MatrixExpr(concat(lhs), concat(𝕀), concat(𝕀), rhs)
end
function (-){T<:Number}(lhs::SDPVar, rhs::AbstractArray{T})
    (size(lhs,1) == size(rhs,1) && size(lhs,2) == size(lhs,2)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(lhs), concat(𝕀), concat(I), -rhs)
end
function (*){T<:Number}(lhs::SDPVar, rhs::AbstractArray{T})
    (size(lhs,2) == size(rhs,1)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(lhs), concat(𝕀), concat(rhs), spzeros(size(lhs,1),size(rhs,2)))
end
# # SDPVar--UniformScaling
# (+)(lhs::SDPVar,rhs::UniformScaling) = MatrixExpr(concat(lhs),concat(𝕀),concat(𝕀), rhs)
# (-)(lhs::SDPVar,rhs::UniformScaling) = MatrixExpr(concat(lhs),concat(𝕀),concat(𝕀),-rhs)
# (*)(lhs::SDPVar,rhs::UniformScaling) = MatrixExpr(concat(lhs),concat(𝕀),concat(rhs),spzeros(size(lhs)...))
# SDPVar--Variable
# SDPVar--AffExpr
# SDPVar--QuadExpr
# SDPVar--MatrixVar
function (+)(lhs::SDPVar,rhs::MatrixVar)
    size(lhs) == size(rhs) || error("Cannot add matrices of unequal size")
    MatrixExpr(concat(lhs,rhs), concat(𝕀,𝕀), concat(𝕀,𝕀), spzeros(size(lhs)...))
end
function (-)(lhs::SDPVar,rhs::MatrixVar)
    size(lhs) == size(rhs) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(lhs,rhs), concat(𝕀,-𝕀), concat(𝕀,𝕀), spzeros(size(lhs)...))
end
# SDPVar--SDPVar
function (+)(lhs::SDPVar, rhs::SDPVar)
    (size(lhs) == size(rhs)) || error("Cannot add matrix variables of unequal size")
    MatrixExpr(concat(lhs,rhs), concat(𝕀,𝕀), concat(𝕀,𝕀), spzeros(size(lhs)...))
end
function (-)(lhs::SDPVar, rhs::SDPVar)
    (size(lhs) == size(rhs)) || error("Cannot subtract matrix variables of unequal size")
    MatrixExpr(concat(lhs,rhs),concat(𝕀,-𝕀), concat(𝕀,𝕀), spzeros(size(lhs)...))
end
# SDPVar--MatrixExpr
(+)(lhs::SDPVar, rhs::MatrixExpr) = MatrixExpr(concat(lhs,rhs.elem),concat(𝕀, rhs.pre),concat(𝕀,rhs.post), rhs.constant)
(-)(lhs::SDPVar, rhs::MatrixExpr) = MatrixExpr(concat(lhs,rhs.elem),concat(𝕀,-rhs.pre),concat(𝕀,rhs.post),-rhs.constant)
function (*)(lhs::SDPVar, rhs::MatrixExpr) 
    (length(rhs.elem) == 0) || error("Cannot multiply matrix variables")
    (size(lhs) == size(rhs.constant)) || error("Cannot multiply matrixes of incompatible sizes")
    MatrixExpr(copy(lhs), copy(rhs.constant), concat(𝕀), spzeros(size(lhs)...))
end

# MatrixExpr
(-)(v::MatrixExpr) = MatrixExpr(copy(v.elem), -v.pre, copy(v.post), copy(v.constant))
# MatrixExpr--Number
# function (+)(lhs::MatrixExpr, rhs::Number)
#     size(lhs) == (1,1) || error("Cannot add scalar to matrix expression")
#     return lhs[1] + rhs
# end
# function (-)(lhs::MatrixExpr, rhs::Number)
#     size(lhs) == (1,1) || error("Cannot add scalar to matrix expression")
#     return lhs[1] - rhs
# end
(*)(lhs::MatrixExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::MatrixExpr, rhs::Number) = (*)(1/rhs,lhs)
# MatrixExpr--AbstractArray{T,2}
function (+){T<:Number}(lhs::MatrixExpr, rhs::AbstractArray{T})
    (size(lhs,1) == size(rhs,1) && size(lhs,2) == size(rhs,2)) || error("Cannot add matrices of unequal size")
    MatrixExpr(copy(lhs.elem), copy(lhs.pre), copy(lhs.post), lhs.constant+rhs)
end
function (-){T<:Number}(lhs::MatrixExpr, rhs::AbstractArray{T})
    (size(lhs,1) == size(rhs,1) && size(lhs,2) == size(rhs,2)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(copy(lhs.elem), copy(lhs.pre), copy(lhs.post), lhs.constant-rhs)
end
function (*){T<:Number}(lhs::MatrixExpr, rhs::AbstractArray{T})
    size(lhs,2) == size(rhs,1) || error("Incompatible dimensions") 
    MatrixExpr(copy(lhs.elem), copy(lhs.pre), map(x->x*rhs, lhs.post), lhs.constant*rhs)
end
# # MatrixExpr--UniformScaling
# (+)(lhs::MatrixExpr,rhs::UniformScaling) = MatrixExpr(copy(lhs.elem),copy(lhs.pre),copy(lhs.post),lhs.constant+rhs)
# (-)(lhs::MatrixExpr,rhs::UniformScaling) = MatrixExpr(copy(lhs.elem),copy(lhs.pre),copy(lhs.post),lhs.constant-rhs)
# (*)(lhs::MatrixExpr,rhs::UniformScaling) = MatrixExpr(copy(lhs.elem),copy(lhs.pre),lhs.post*rhs.λ,lhs.constant*rhs.λ)
# MatrixExpr--Variable
function (+)(lhs::MatrixExpr, rhs::Variable)
    size(lhs) == (1,1) || error("Cannot add scalar to matrix expression")
    return lhs[1] + rhs
end
function (-)(lhs::MatrixExpr, rhs::Variable)
    size(lhs) == (1,1) || error("Cannot add scalar to matrix expression")
    return lhs[1] - rhs
end
# MatrixExpr--AffExpr
function (+)(lhs::MatrixExpr, rhs::AffExpr)
    size(lhs) == (1,1) || error("Cannot add scalar to matrix expression")
    return lhs[1] + rhs
end
function (-)(lhs::MatrixExpr, rhs::AffExpr)
    size(lhs) == (1,1) || error("Cannot add scalar to matrix expression")
    return lhs[1] - rhs
end
# MatrixExpr--QuadExpr
# MatrixExpr--MatrixVar
function (+)(lhs::MatrixExpr,rhs::MatrixVar)
    size(lhs) == size(rhs) || error("Cannot add matrices of unequal size")
    MatrixExpr(concat(lhs.elem,rhs), concat(lhs.pre,𝕀), concat(lhs.post,𝕀), copy(lhs.constant))
end
function (-)(lhs::MatrixExpr,rhs::MatrixVar)
    size(lhs) == size(rhs) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(concat(lhs.elem,rhs), concat(lhs.pre,-𝕀), concat(lhs.post,𝕀), lhs.constant)
end
function (*)(lhs::MatrixExpr,rhs::MatrixVar)
    (m,p),(q,n) = size(lhs), size(rhs)
    p == q || error("Imcompatible dimensions for matrix multiplication")
    m == n == 1 || error("Only quadratic forms are allowed")
    ret = QuadExpr()
    for k in 1:p
        ret += lhs[k]*rhs[k]
    end
    return ret
end
# MatrixExpr--SDPVar
function (+)(lhs::MatrixExpr, rhs::SDPVar)
    (size(lhs.constant) == size(rhs)) || error("Cannot add matrix variables of unequal size")
    MatrixExpr(concat(lhs.elem,rhs), concat(lhs.pre,𝕀), concat(lhs.post,𝕀), lhs.constant)
end
function (-)(lhs::MatrixExpr, rhs::SDPVar)
    (size(lhs.constant) == size(rhs)) || error("Cannot subtract matrix variables of unequal size")
    MatrixExpr(concat(lhs.elem,rhs), concat(lhs.pre,-𝕀), concat(lhs.post,𝕀), lhs.constant)
end
# MatrixExpr--MatrixExpr
(+)(lhs::MatrixExpr, rhs::MatrixExpr) = MatrixExpr(concat(lhs.elem,rhs.elem),concat(lhs.pre, rhs.pre),concat(lhs.post,rhs.post),lhs.constant+rhs.constant)
(-)(lhs::MatrixExpr, rhs::MatrixExpr) = MatrixExpr(concat(lhs.elem,rhs.elem),concat(lhs.pre,-rhs.pre),concat(lhs.post,rhs.post),lhs.constant-rhs.constant)
function (*)(lhs::MatrixExpr,rhs::MatrixExpr)
    (m,p),(q,n) = size(lhs), size(rhs)
    p == q || error("Imcompatible dimensions for matrix multiplication")
    m == n == 1 || error("Only quadratic forms are allowed")
    ret = QuadExpr()
    for k in 1:p
        ret += lhs[k]*rhs[k]
    end
    return ret
end
# MatrixFuncVar
(-)(lhs::MatrixFuncVar) = MatrixFuncVar(-lhs.expr, lhs.func)
# MatrixFuncVar--Number
(+)(lhs::MatrixFuncVar,rhs::Number) = ScalarExpr(concat(lhs), AffExpr( rhs), NormExpr())
(-)(lhs::MatrixFuncVar,rhs::Number) = ScalarExpr(concat(lhs), AffExpr(-rhs), NormExpr())
(*)(lhs::MatrixFuncVar,rhs::Number) = MatrixFuncVar(lhs.expr*rhs,lhs.func)
(/)(lhs::MatrixFuncVar,rhs::Number) = MatrixFuncVar(lhs.expr/rhs,lhs.func)
# MatrixFuncVar--Variable
(+)(lhs::MatrixFuncVar,rhs::Variable) = ScalarExpr(concat(lhs),convert(AffExpr,rhs), NormExpr())
(-)(lhs::MatrixFuncVar,rhs::Variable) = ScalarExpr(concat(lhs),               -rhs , NormExpr())
# MatrixFuncVar--AffExpr
(+)(lhs::MatrixFuncVar,rhs::AffExpr) = ScalarExpr(concat(lhs), copy(rhs), NormExpr())
(-)(lhs::MatrixFuncVar,rhs::AffExpr) = ScalarExpr(concat(lhs),     -rhs , NormExpr())
# MatrixFuncVar--AbstractArray{T,2}
function (*){T<:Number}(lhs::MatrixFuncVar,rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix must be symmetric")
    DualExpr(concat(lhs), concat(rhs), spzeros(size(rhs)...))
end
# MatrixFuncVar--MatrixFuncVar
(+)(lhs::MatrixFuncVar,rhs::MatrixFuncVar) = ScalarExpr(concat(lhs, rhs), AffExpr(), NormExpr())
(-)(lhs::MatrixFuncVar,rhs::MatrixFuncVar) = ScalarExpr(concat(lhs,-rhs), AffExpr(), NormExpr())
# MatrixFuncVar--NormExpr
(+)(lhs::MatrixFuncVar,rhs::NormExpr) = ScalarExpr(concat(lhs), AffExpr(), copy(rhs))
(-)(lhs::MatrixFuncVar,rhs::NormExpr) = ScalarExpr(concat(lhs), AffExpr(),     -rhs )
# MatrixFuncVar--ScalarExpr
(+)(lhs::MatrixFuncVar,rhs::ScalarExpr) = ScalarExpr(concat(lhs,      rhs.matvars ), copy(rhs.aff), copy(rhs.normexpr))
(-)(lhs::MatrixFuncVar,rhs::ScalarExpr) = ScalarExpr(concat(lhs,map(-,rhs.matvars)),     -rhs.aff ,     -rhs.normexpr )

# NormExpr
(-)(v::NormExpr) = NormExpr(copy(v.vars), -v.coeffs, copy(v.form))
# NormExpr--Number
(+)(lhs::NormExpr,rhs::Number)  = ScalarExpr(MatrixFuncVar[], AffExpr( rhs), lhs)
(-)(lhs::NormExpr,rhs::Number)  = ScalarExpr(MatrixFuncVar[], AffExpr(-rhs), lhs)
(*)(lhs::NormExpr,rhs::Number)  = 
    ( rhs == 0 ? NormExpr(lhs.form) : NormExpr(copy(lhs.vars), lhs.coeffs*rhs, lhs.form) )
(/)(lhs::NormExpr,rhs::Number)  = NormExpr(lhs.vars, lhs.coeffs/rhs, lhs.form)
# NormExpr--Variable
(+)(lhs::NormExpr,rhs::Variable) = ScalarExpr(MatrixFuncVar[], convert(AffExpr, rhs), lhs)
(-)(lhs::NormExpr,rhs::Variable) = ScalarExpr(MatrixFuncVar[],                 -rhs,  lhs)
# NormExpr--AffExpr
(+)(lhs::NormExpr,rhs::AffExpr) = ScalarExpr(MatrixFuncVar[],  rhs, lhs)
(-)(lhs::NormExpr,rhs::AffExpr) = ScalarExpr(MatrixFuncVar[], -rhs, lhs)
# NormExpr--NormExpr
(+)(lhs::NormExpr,rhs::NormExpr) = NormExpr(concat(lhs.vars,rhs.vars), vcat(lhs.coeffs, rhs.coeffs), vcat(lhs.form,rhs.form))
(-)(lhs::NormExpr,rhs::NormExpr) = NormExpr(concat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), vcat(lhs.form,rhs.form))
# NormExpr--MatrixFuncVar
(+)(lhs::NormExpr,rhs::MatrixFuncVar) = ScalarExpr(concat(+rhs), AffExpr(), copy(lhs))
(-)(lhs::NormExpr,rhs::MatrixFuncVar) = ScalarExpr(concat(-rhs), AffExpr(), copy(lhs))
# NormExpr-ScalarExpr
(+)(lhs::NormExpr,rhs::ScalarExpr) = ScalarExpr( copy(rhs.matvars), copy(rhs.aff), lhs+rhs.normexpr)
(-)(lhs::NormExpr,rhs::ScalarExpr) = ScalarExpr(map(-,rhs.matvars),     -rhs.aff , lhs-rhs.normexpr)

# ScalarExpr
(-)(lhs::ScalarExpr) = ScalarExpr(map(-,lhs.matvars), -lhs.aff, -lhs.normexpr)
# ScalarExpr--Number
(+)(lhs::ScalarExpr,rhs::Number) = ScalarExpr(        copy(lhs.matvars), lhs.aff+rhs, copy(lhs.normexpr))
(-)(lhs::ScalarExpr,rhs::Number) = ScalarExpr(        copy(lhs.matvars), lhs.aff-rhs, copy(lhs.normexpr))
(*)(lhs::ScalarExpr,rhs::Number) = ScalarExpr(map(x->x*rhs,lhs.matvars), lhs.aff*rhs, lhs.normexpr*rhs)
(/)(lhs::ScalarExpr,rhs::Number) = ScalarExpr(map(x->x/rhs,lhs.matvars), lhs.aff/rhs, lhs.normexpr/rhs)
# ScalarExpr--Variable
(+)(lhs::ScalarExpr,rhs::Variable) = ScalarExpr(copy(lhs.matvars), lhs.aff+rhs, copy(lhs.normexpr))
(-)(lhs::ScalarExpr,rhs::Variable) = ScalarExpr(copy(lhs.matvars), lhs.aff-rhs, copy(lhs.normexpr))
# ScalarExpr--AffExpr
(+)(lhs::ScalarExpr,rhs::AffExpr) = ScalarExpr(copy(lhs.matvars), lhs.aff+rhs, copy(lhs.normexpr))
(-)(lhs::ScalarExpr,rhs::AffExpr) = ScalarExpr(copy(lhs.matvars), lhs.aff-rhs, copy(lhs.normexpr))
# ScalarExpr--AbstractArray{T,2}
function (*){T<:Number}(lhs::ScalarExpr,rhs::AbstractArray{T,2})
    isequal(lhs.normexpr,NormExpr()) || error("Cannot multiply norm expression by a matrix")
    issym(rhs) || error("Matrix must by symmetric")
    DualExpr(concat(lhs.matvars,lhs.aff.vars), concat([rhs for c in lhs.matvars],[c*rhs for c in lhs.aff.coeffs]), lhs.aff.constant*rhs)
end
# ScalarExpr--MatrixFuncVar
(+)(lhs::ScalarExpr,rhs::MatrixFuncVar) = ScalarExpr(concat(lhs.matvars, rhs), copy(lhs.aff), copy(lhs.normexpr))
(-)(lhs::ScalarExpr,rhs::MatrixFuncVar) = ScalarExpr(concat(lhs.matvars,-rhs), copy(lhs.aff), copy(lhs.normexpr))
# ScalarExpr-NormExpr
(+)(lhs::ScalarExpr,rhs::NormExpr) = ScalarExpr(copy(lhs.matvars), copy(lhs.aff), lhs.normexpr+rhs)
(-)(lhs::ScalarExpr,rhs::NormExpr) = ScalarExpr(copy(lhs.matvars), copy(lhs.aff), lhs.normexpr-rhs)

# ScalarExpr--ScalarExpr
(+)(lhs::ScalarExpr,rhs::ScalarExpr) = ScalarExpr(concat(lhs.matvars, rhs.matvars), lhs.aff+rhs.aff, lhs.normexpr+rhs.normexpr)
(-)(lhs::ScalarExpr,rhs::ScalarExpr) = ScalarExpr(concat(lhs.matvars,-rhs.matvars), lhs.aff-rhs.aff, lhs.normexpr-rhs.normexpr)

for sgn in (:<=, :(==), :>=, :(.<=), :(.>=))
    @eval begin 
        # Number
        for rtype in [:MatrixVar, :SDPVar, :MatrixExpr]
            @eval begin
                function $(sgn)(lhs::Number, rhs::MatrixVar)
                    lhs == 0.0 || error("Cannot use scalar bound unless it is zero")
                    MatrixConstraint(-rhs, $(quot(sgn)))
                end
            end
        end

        # Number--DualExpr
        function $(sgn)(lhs::Number,rhs::DualExpr)
            lhs == 0.0 || error("Cannot use scalar bound unless it is zero")
            DualConstraint(-rhs, $(quot(sgn)))
        end
        # AbstractArray{T,2}
        # AbstractArray{T,2}--DualExpr
        $(sgn){T<:Number}(lhs::AbstractArray{T,2}, rhs::DualExpr) = DualConstraint(lhs-rhs, $(quot(sgn)))
        # AbstractArray{T,2}--MatrixVar
        $(sgn){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixVar) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # AbstractArray{T,2}--SDPVar
        $(sgn){T<:Number}(lhs::AbstractArray{T,2}, rhs::SDPVar) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # AbstractArray{T,2}--MatrixExpr
        $(sgn){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixExpr) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # MatrixVar
        # MatrixVar--Number
        function $(sgn)(lhs::MatrixVar, rhs::Number)
            rhs == 0.0 || error("Incompatible dimensions")
            MatrixConstraint(lhs, $(quot(sgn)))
        end
        # MatrixVar--AbstractArray{T}
        $(sgn){T<:Number}(lhs::MatrixVar, rhs::AbstractArray{T}) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # MatrixVar--MatrixVar
        $(sgn)(lhs::MatrixVar, rhs::MatrixVar) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # MatrixVar--SDPVar
        $(sgn)(lhs::MatrixVar, rhs::SDPVar) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # MatrixVar--MatrixExpr
        # $(sgn)(lhs::MatrixVar, rhs::MatrixExpr) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # SDPVar
        # SDPVar--Number
        function $(sgn)(lhs::SDPVar, rhs::Number)
            rhs == 0.0 || error("Incompatible dimensions")
            MatrixConstraint(lhs, $(quot(sgn)))
        end
        function $(sgn)(lhs::MatrixVar, rhs::Number)
            rhs == 0.0 || error("Cannot use scalar bound unless it is zero")
            MatrixConstraint(lhs, $(quot(sgn)))
        end
        # SDPVar--AbstractArray{T,2}
        $(sgn){T<:Number}(lhs::SDPVar, rhs::AbstractArray{T,2}) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # SDPVar--UniformScaling
         $(sgn)(lhs::SDPVar, rhs::UniformScaling) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # SDPVar--SDPVar
        $(sgn)(lhs::SDPVar, rhs::SDPVar) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # SDPVar--MatrixExpr
        $(sgn)(lhs::SDPVar, rhs::MatrixExpr) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # MatrixExpr
        # MatrixExpr--Number
        # function $(sgn)(lhs::MatrixExpr,rhs::Number)
            # rhs == 0.0 || error("Cannot use scalar bound unless it is zero")
            # MatrixConstraint(lhs, $(quot(sgn)))
        # end
        # MatrixExpr--AbstractArray{T,2}
        # $(sgn){T<:Number}(lhs::MatrixExpr, rhs::AbstractArray{T,2}) = 
            # MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # MatrixExpr--UniformScaling
        $(sgn)(lhs::MatrixExpr, rhs::UniformScaling) =
            MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # MatrixExpr--SDPVar
        $(sgn)(lhs::MatrixExpr, rhs::SDPVar) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # MatrixExpr--MatrixExpr
        # $(sgn)(lhs::MatrixExpr, rhs::MatrixExpr) = MatrixConstraint(lhs-rhs, $(quot(sgn)))
        # DualExpr--Number
        function $(sgn)(lhs::DualExpr,rhs::Number)
            rhs == 0.0 || error("Cannot use scalar bound unless it is zero")
            DualConstraint(lhs, $(quot(sgn)))
        end
        # DualExpr--AbstractArray{T,2}
        $(sgn){T<:Number}(lhs::DualExpr, rhs::AbstractArray{T,2}) = DualConstraint(lhs-rhs, $(quot(sgn)))
        # DualExpr--DualExpr
        $(sgn)(lhs::DualExpr, rhs::DualExpr) = DualConstraint(lhs-rhs, $(quot(sgn)))
    end
end

for ltype in (:Number, :Variable, :AffExpr, :MatrixFuncVar, :NormExpr, :ScalarExpr)
    for rtype in (:Number, :Variable, :AffExpr, :MatrixFuncVar, :NormExpr, :ScalarExpr)
        if !(ltype in [:Number, :Variable, :AffExpr]) ||
           !(rtype in [:Number, :Variable, :AffExpr]) # don't overwrite existing ones
            @eval begin
                (<=)(lhs::$(ltype), rhs::$(rtype)) = PrimalConstraint(lhs-rhs, -Inf, 0.0)
                (==)(lhs::$(ltype), rhs::$(rtype)) = PrimalConstraint(lhs-rhs,  0.0, 0.0)
                (>=)(lhs::$(ltype), rhs::$(rtype)) = PrimalConstraint(lhs-rhs,  0.0, Inf)

            end
        end
    end
end

function (<=)(lhs::Number, rhs::MatrixExpr)
    !(lhs == 0) && !(issym(rhs)) && error("Incompatible dimensions")
    return MatrixConstraint(rhs, :>=)
end
function (<=)(lhs::MatrixExpr, rhs::Number)
    !(rhs == 0) && !(issym(lhs)) && error("Incompatible dimensions")
    return MatrixConstraint(-lhs, :>=)
end
function (==)(lhs::Number, rhs::MatrixExpr)
    !(lhs == 0) && !(issym(rhs)) && error("Incompatible dimensions")
    return MatrixConstraint(rhs, :(==))
end
function (==)(lhs::MatrixExpr, rhs::Number)
    rhs == 0 || error("Incompatible dimensions")
    return MatrixConstraint(lhs, :(==))
end
function (>=)(lhs::Number, rhs::MatrixExpr)
    !(lhs == 0) && !(issym(rhs)) && error("Incompatible dimensions")
    return MatrixConstraint(-rhs, :>=)
end
function (>=)(lhs::MatrixExpr, rhs::Number)
    !(rhs == 0) && !(issym(lhs)) && error("Incompatible dimensions")
    return MatrixConstraint(lhs, :>=)
end
function (.<=)(lhs::Number, rhs::MatrixExpr)
    lhs == 0 || error("Incompatible dimensions")
    return MatrixConstraint(rhs, :.>=)
end
function (.<=)(lhs::MatrixExpr, rhs::Number)
    rhs == 0 || error("Incompatible dimensions")
    return MatrixConstraint(lhs, :.<=)
end
function (.>=)(lhs::Number, rhs::MatrixExpr)
    lhs == 0 || error("Incompatible dimensions")
    return MatrixConstraint(rhs, :.<=)
end
function (.>=)(lhs::MatrixExpr, rhs::Number)
    rhs == 0 || error("Incompatible dimensions")
    return MatrixConstraint(lhs, :.>=)
end

# >= and <= reserved for SDP constraints, .>= and .<= for componentwise
for ltype in (:AbstractArray, :Variable, :AffExpr, :MatrixVar, :MatrixExpr, :MatrixFuncVar, :NormExpr, :ScalarExpr)
    @eval begin
        function (<=)(lhs::$(ltype), rhs::MatrixExpr)
            terms = rhs - lhs
            issym(terms) || error("Semidefinite constraint must be symmetric")
            return MatrixConstraint(terms, :>=)
        end
        function (<=)(lhs::MatrixExpr, rhs::$(ltype))
            terms = rhs - lhs
            issym(terms) || error("Semidefinite constraint must be symmetric")
            return MatrixConstraint(terms, :>=)
        end
        function (>=)(lhs::$(ltype), rhs::MatrixExpr)
            terms = lhs - rhs
            issym(terms) || error("Semidefinite constraint must be symmetric")
            return MatrixConstraint(terms, :>=)
        end
        function (>=)(lhs::MatrixExpr, rhs::$(ltype))
            terms = lhs - rhs
            issym(terms) || error("Semidefinite constraint must be symmetric")
            return MatrixConstraint(terms, :>=)
        end
        (==)(lhs::$(ltype), rhs::MatrixExpr) = MatrixConstraint(rhs-lhs, :(==))
        (==)(lhs::MatrixExpr, rhs::$(ltype)) = MatrixConstraint(rhs-lhs, :(==))
        (.<=)(lhs::$(ltype), rhs::MatrixExpr) = MatrixConstraint(rhs-lhs, :.>=)
        (.<=)(lhs::MatrixExpr, rhs::$(ltype)) = MatrixConstraint(lhs-rhs, :.<=)
        (.>=)(lhs::$(ltype), rhs::MatrixExpr) = MatrixConstraint(rhs-lhs, :.<=)
        (.>=)(lhs::MatrixExpr, rhs::$(ltype)) = MatrixConstraint(lhs-rhs, :.>=)
    end
end

function Base.dot(x::AbstractArray, y::MatrixExpr)
    (size(x,1) == size(y,1) && size(x,2) == size(y,2)) || error("Incompatible sizes")
    ret = ScalarExpr()
    for i in 1:size(x,1), j in 1:size(x,2)
        ret += x[i,j] * y[i,j]
    end
    return ret
end
Base.dot(x::MatrixExpr, y::AbstractArray) = dot(y,x)
Base.dot(x::MatrixVar,y::AbstractArray) = dot(convert(MatrixExpr,x),y)
Base.dot(x::AbstractArray, y::MatrixVar) = dot(y,x)
