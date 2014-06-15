export SDPData, 
       SDPModel, 
       SDPVar, 
       MatrixVar,
       BlockMatrix,
       MatrixExpr, 
       MatrixFuncVar,
       NormExpr,
       ScalarExpr,
       DualExpr,
       PrimalConstraint,
       DualConstraint,
       MatrixConstraint,
       addConstraints,
       @defSDPVar,
       @defMatrixVar

type SolverInfo
    id # internal solver IDs for SDPVar
       # column indices for MatrixVar (Range)
    psd::Bool # psd==true, nsd==false
    offset # constant C in X≽C
end

SolverInfo() = SolverInfo(0,true,0.0)

type SDPData
    sdpvar#::Vector{SDPVar}
    varname::Vector{String}
    lb
    ub
    solverinfo::Vector{SolverInfo}
    matrixconstr#::Vector{MatrixConstraint}
    primalconstr#::Vector{PrimalConstraint}
    dualconstr#::Vector{DualConstraint}
    sdpobj
    sdpval
end

SDPData() = SDPData({}, String[], {}, {}, SolverInfo[], MatrixConstraint[], PrimalConstraint[], DualConstraint[], ScalarExpr(), {})

initSDP(m::Model) = (m.sdpdata == nothing) && (m.sdpdata = SDPData())

const snsmap = [:(>=) => "≽", :(==) => "=", :(<=) => "≼", :(.>=) => "≥", :(.<=) => "≤"]

if Pkg.installed("Mosek") != nothing
    eval(Expr(:using,:Mosek))
end

const 𝕀 = UniformScaling(1)

# Useful type hierarchy for vcat, hcat, etc.
abstract SDPMatrix 

###############################################################################
# Semidefinite Variable class
# Pointer to model, with solver index and dimension
type SDPVar <: SDPMatrix
    m::Model
    index::Int64
    dim::Int64
end

Base.transpose(a::SDPVar)  = a
Base.ctranspose(a::SDPVar) = a
Base.conj(a::SDPVar) = a

Base.size(d::SDPVar) = (d.dim,d.dim)
Base.size(d::SDPVar, slice::Int64) = (0 <= slice <= 2) ? d.dim : 1
Base.ndims(d::SDPVar) = 2
Base.eye(d::SDPVar)  = eye(d.dim)
Base.issym(d::SDPVar) = true
Base.isequal(a::SDPVar, b::SDPVar) = isequal(a.m,b.m) && isequal(a.index,b.index)

function Base.diagm(d::SDPVar)
    m,n = size(d)
    return DualExpr(MatrixFuncVar[d[it,it] for it in 1:m],
                   SparseMatrixCSC[sparse([it],[it],[1.0],m,n) for it in 1:m],
                   spzeros(m,n))
#    return mapreduce(it->d[it,it]*sparse([it],[it],[1.0],m,n), +, 1:m)
end
Base.spdiagm(d::SDPVar) = diagm(d)

getValue(d::SDPVar) = d.m.sdpdata.sdpval[d.index] 

getLower(d::SDPVar) = d.m.sdpdata.lb[d.index]
getUpper(d::SDPVar) = d.m.sdpdata.ub[d.index]

getName(d::SDPVar) = d.m.sdpdata.varname[d.index]
setName(d::SDPVar, name::String) = (d.m.sdpdata.varname[d.index] = name)

Base.trace(c::SDPVar)  = trace(convert(MatrixExpr, c))
Base.dot(c::SDPVar,d::AbstractArray) = trace(c*d)
Base.dot(c::AbstractArray,d::SDPVar) = trace(c*d)
Base.norm(c::SDPVar)   = norm(convert(MatrixExpr, c))
Base.sum(c::SDPVar)    = sum(convert(MatrixExpr, c))

Base.show(io::IO,d::SDPVar)  = print(io, "$(d.m.sdpdata.varname[d.index])")
Base.print(io::IO,d::SDPVar) = println(io, "$(d.m.sdpdata.varname[d.index]) ∈ 𝒮₊($(d.dim))")

Base.getindex(d::SDPVar, x::Int64, y::Int64) = 
    MatrixFuncVar(MatrixExpr(concat(d),concat(sparse([x,y],[y,x],[0.5,0.5],d.dim,d.dim)),concat(𝕀),spzeros(d.dim,d.dim)),:ref)

macro defSDPVar(m, x, extra...)
    m = esc(m)
    if isexpr(x,:comparison)
        # we have some bounds
        if x.args[2] == :>=
            if length(x.args) == 5
                error("Use the form lb <= var <= ub instead of ub >= var >= lb")
            end
            @assert length(x.args) == 3
            # lower bounds, no upper
            lb = esc(x.args[3])
            ub = Inf
            var = x.args[1]
        elseif x.args[2] == :<=
            if length(x.args) == 5
                # lb <= x <= u
                lb = esc(x.args[1])
                if (x.args[4] != :<=)
                    error("Expected <= operator")
                end
                ub = esc(x.args[5])
                var = x.args[3]
            else
                # x <= u
                ub = esc(x.args[3])
                lb = -Inf
                var = x.args[1]
            end
        end
    else
        var = x
        lb = 0.0
        ub = Inf
    end
    if length(extra) > 0
        # TODO: allow user to specify matrix properties here (diagonal, etc.)
    end

    if length(var.args) == 3
        var.args[2] == var.args[3] || error("SDP variable must be square")
    elseif length(var.args) > 3
        error("SDPVar must be specificed by one or two (equal) dimensions")
    end

    if !isexpr(var,:ref) # TODO: infer size via syntax like @defSDPVar(m, X >= ones(3,3))
        error("Syntax error: Need to specify matrix size (e.g. $var[5])")
    else
        varname = esc(var.args[1])
        sz = esc(var.args[2])
        code = quote
            issym($lb) || error("Lower bound is not symmetric")
            issym($ub) || error("Upper bound is not symmetric")
            (isa($lb,AbstractArray) && !($sz == Base.size($lb,1))) && error("Lower bound is not of same size as variable")
            (isa($ub,AbstractArray) && !($sz == Base.size($ub,1))) && error("Upper bound is not of same size as variable")
            isa($lb,Number) && !($lb == 0.0 || $lb == -Inf) && error("Bounds must be of same size as variable")
            isa($ub,Number) && !($ub == 0.0 || $ub ==  Inf) && error("Bounds must be of same size as variable")
            $lb == -Inf && $ub == Inf && error("Replace unrestricted SDP variable with regular, scalar variables")
            $lb == $ub && error("Replace with constant matrix")
            initSDP($m)
            sdp = $(m).sdpdata
            $(varname) = SDPVar($m, length(sdp.sdpvar)+1, $(esc(var.args[2])))
            push!(sdp.sdpvar, $(varname))
            push!(sdp.lb, $lb)
            push!(sdp.ub, $ub)
            push!(sdp.varname, $(string(var.args[1])))
            push!(sdp.solverinfo, SolverInfo())
            nothing
        end
        return code
    end
end

###############################################################################
# Matrix Variable class
# Pointer to model, with solver index and dimension
type MatrixVar <: SDPMatrix
    m::Model
    index::Int64
    dim::(Int64,Int64)
end

function MatrixVar(m::Model, sz::(Int,Int), lb, ub, varname::String)
    sdp = m.sdpdata
    push!(sdp.lb, lb)
    push!(sdp.ub, ub)
    push!(sdp.varname, varname)
    rng = (m.numCols+1):(m.numCols+prod(sz))
    patt = (sz[2] == 1) ? :col : 
           (sz[1] == 1) ? :row : :mat
    for j in 1:sz[2], i in 1:sz[1]
        lower = (lb == -Inf || lb == 0.0) ? lb : lb[i,j]
        upper = (ub ==  Inf || ub == 0.0) ? ub : ub[i,j]
        vname = (patt == :col) ? "$varname[$i]" : 
                (patt == :row) ? "$varname[1,$j]" : "$varname[$i,$j]"
        Variable(m, lower, upper, JuMP.CONTINUOUS, vname) # TODO: make this more efficient
    end
    push!(sdp.solverinfo, SolverInfo(rng,true,0.0))
    return MatrixVar(m, length(sdp.sdpvar)+1, sz)
end

Base.transpose(a::MatrixVar)  = MatrixVar(a.m,a.index,reverse(a.dim))
Base.ctranspose(a::MatrixVar) = MatrixVar(a.m,a.index,reverse(a.dim))
Base.conj(a::MatrixVar) = a

Base.size(a::MatrixVar) = a.dim
Base.size(d::MatrixVar, slice::Int64) = (0 <= slice <= 2) ? d.dim[slice] : 1
Base.ndims(d::MatrixVar) = 2
Base.eye(d::MatrixVar)  = eye(d.dim[1],d.dim[2])
Base.issym(d::MatrixVar) = (d.dim == (1,1))
Base.isequal(a::MatrixVar, b::MatrixVar) = (a.m == b.m) && isequal(a.index,b.index) && isequal(a.dim,b.dim)

function Base.dot{T<:Number}(a::Array{T},b::MatrixVar)
    (size(a,1) == size(b,1) && size(a,2) == size(b,2)) || error("Incompatible dimensions")
    return AffExpr(Variable[b[it] for it in 1:length(a)],
                   Float64[a[it] for it in 1:length(a)],
                   0.0)
    # return mapreduce(it->a[it]*b[it], +, 1:length(a))
end
Base.dot{T<:Number}(a::MatrixVar,b::Array{T}) = dot(b,a)

function Base.diagm(d::MatrixVar)
    m,n = size(d)
    m == n || error("Cannot construct diagonal of nonsquare matrix")
    return DualExpr(Variable[d[it,it] for it in 1:m],
               SparseMatrixCSC[sparse([it],[it],[1.0],m,n) for it in 1:m],
               spzeros(m,n))
    # return mapreduce(it->d[it,it]*sparse([it],[it],[1.0],m,n), +, 1:m)
end
Base.spdiagm(d::MatrixVar) = diagm(d)

getValue(d::MatrixVar) = reshape(d.m.sdpdata.sdpval[d.index], size(d))

getLower(d::MatrixVar) = d.m.sdpdata.lb[d.index]
getUpper(d::MatrixVar) = d.m.sdpdata.ub[d.index]

getName(d::MatrixVar) = d.m.sdpdata.varname[d.index]
setName(d::MatrixVar, name::String) = (d.m.sdpdata.varname[d.index] = name)

Base.show(io::IO,d::MatrixVar)  = print(io, "$(d.m.sdpdata.varname[d.index])")
Base.print(io::IO,d::MatrixVar) = println(io, "$(d.m.sdpdata.varname[d.index]) ∈ ℝ($(d.dim[1])×$(d.dim[2]))")

function Base.getindex(d::MatrixVar, x::Int)
    rng = d.m.sdpdata.solverinfo[d.index].id
    Variable(d.m, rng[x])
end

function Base.getindex(d::MatrixVar, x::Int, y::Int)
    sx,sy = d.dim
    rng   = d.m.sdpdata.solverinfo[d.index].id
    off = (y-1)*sx + x
    return Variable(d.m, rng[off])
end

function Base.getindex(d::MatrixVar, x::Int, y::Range{Int})
    sx,sy = d.dim
    rng   = d.m.sdpdata.solverinfo[d.index].id
    arr = Array(Variable, 1, length(y))
    for (it,val) in enumerate(y)
        off = (val-1)*sx + x
        arr[it] = Variable(d.m, rng[off])
    end
    return arr
end

function Base.getindex(d::MatrixVar, x::Range{Int}, y::Int)
    sx,sy = d.dim
    rng   = d.m.sdpdata.solverinfo[d.index].id
    arr = Array(Variable, length(x))
    for (it,val) in enumerate(x)
        off = (y-1)*sx + val
        arr[it] = Variable(d.m, rng[off])
    end
    return arr
end

macro defMatrixVar(m, x, extra...)
    m = esc(m)
    if isexpr(x,:comparison)
        # we have some bounds
        if x.args[2] == :>=
            if length(x.args) == 5
                error("Use the form lb <= var <= ub instead of ub >= var >= lb")
            end
            @assert length(x.args) == 3
            # lower bounds, no upper
            lb = esc(x.args[3])
            ub = Inf
            var = x.args[1]
        elseif x.args[2] == :<=
            if length(x.args) == 5
                # lb <= x <= u
                lb = esc(x.args[1])
                if (x.args[4] != :<=)
                    error("Expected <= operator")
                end
                ub = esc(x.args[5])
                var = x.args[3]
            else
                # x <= u
                ub = esc(x.args[3])
                lb = -Inf
                var = x.args[1]
            end
        end
    else
        var = x
        lb = -Inf
        ub =  Inf
    end
    if length(extra) > 0
        # TODO: allow user to specify matrix properties here (diagonal, etc.)
    end

    if isa(var,Symbol)
        error("Consider using defVar for scalar variables")
    elseif length(var.args) == 2
        Base.warn_once("Single dimension defaults to column vector")
        sx = esc(var.args[2])
        sy = 1
    elseif length(var.args) == 3
        sx = esc(var.args[2])
        sy = esc(var.args[3])
    else
        error("MatrixVar must be specificed by less than two dimensions")
    end

    if !isexpr(var,:ref)
        error("Syntax error: Need to specify matrix size (e.g. $var[5])")
    else
        varname = esc(var.args[1])
        code = quote
            (isa($lb,AbstractArray) && !($sx == Base.size($lb,1) && $sy == Base.size($lb,2))) && error("Lower bound is not of same size as variable")
            (isa($ub,AbstractArray) && !($sx == Base.size($ub,1) && $sy == Base.size($ub,2))) && error("Upper bound is not of same size as variable")
            isa($lb,Number) && !($lb == 0.0 || $lb == -Inf) && error("Bounds must be of same size as variable")
            isa($ub,Number) && !($ub == 0.0 || $ub ==  Inf) && error("Bounds must be of same size as variable")
            $lb == $ub && error("Replace with constant matrix")
            initSDP($m)
            sdp = $(m).sdpdata
            $(varname) = MatrixVar($m, ($sx,$sy), $lb, $ub, $(string(var.args[1])))
            # $(varname) = MatrixVar($m, length(sdp.sdpvar)+1, ($sx,$sy))
            push!(sdp.sdpvar, $(varname))
            # push!(sdp.lb, $lb)
            # push!(sdp.ub, $ub)
            # push!(sdp.varname, $(string(var.args[1])))
            # push!(sdp.solverinfo, SolverInfo())
            nothing
        end
        return code
    end
end

###############################################################################
# Block Matrix Class
# Dense container for matrix objects. Elements can be a matrix expression, 
# matrix variable, SDP variable, Variable, array, or number. 
type BlockMatrix <: SDPMatrix
    elem::Array
    sz::(Int,Int)
end

function BlockMatrix(array::Array)
    bm, bn = blocksize(array) # do this to assert dimensions are consistent
    return BlockMatrix(array,(bm,bn))
end

function blocksize(array::Array)
    m, n = size(array)
    sx = Array(Int64, m, n)
    sy = Array(Int64, m, n)
    for i in 1:m
        for j in 1:n
            elem = array[i,j]
            if isa(elem, UniformScaling)
                error("Not yet support for UniformScaling")
            elseif isa(elem, Number) || isa(elem, Variable) || isa(elem, AffExpr) || isa(elem, QuadExpr)
                sx[i,j] = 1
                sy[i,j] = 1
            else
                sx[i,j] = size(elem,2)
                sy[i,j] = size(elem,1)
            end
        end
    end
    for i in 1:m
        for j in 2:n
            if sy[i,j] != sy[i,j-1]
                error("Incompatible number of rows")
            end
        end
    end
    for j in 1:n
        for i in 2:m
            if sx[i,j] != sx[i-1,j]
                error("Incompatible number of columns")
            end
        end
    end
    return sum(sy[:,1]), sum(sx[1,:])
end

Base.size(b::BlockMatrix) = b.sz
blocksize(b::BlockMatrix) = size(b.elem)

Base.transpose(b::BlockMatrix)  = BlockMatrix(transpose(map(transpose,b.elem)))
Base.ctranspose(b::BlockMatrix) = BlockMatrix(transpose(map(transpose,b.elem)))

Base.isequal(a::BlockMatrix, b::BlockMatrix) = isequal(a.elem,b.elem)

function Base.issym(d::BlockMatrix)
    m,n = size(d.elem)
    m == n || return false
    for i in 1:n # check that diagonal is symmetric
        if !issym(d.elem[i,i])
            return false
        end
    end
    for i = 1:(n-1), j = (i+1):n # check off-diagonal blocks are transposes of each other
        if !isequal(d.elem[i,j], transpose(d.elem[j,i]))
            return false
        end
    end
    return true
end

Base.size(::Variable) = (1,)
Base.size(::Variable,::Int) = 1
Base.size(::AffExpr) = (1,1)
Base.size(::AffExpr,::Int) = 1
Base.size(::QuadExpr) = (1,)
Base.size(::QuadExpr,::Int) = 1

function Base.getindex(d::BlockMatrix, x::Int64, y::Int64)
    m,n = size(d)
    curr = 0
    it = 1
    while x > (curr+size(d.elem[it,1],1))
        curr += size(d.elem[it,1],1)
        it += 1
    end
    idx = it
    offx = x - curr
    curr = 0
    it = 1
    while y > (curr+size(d.elem[1,it],2))
        curr += size(d.elem[1,it],2)
        it += 1
    end
    idy = it
    offy = y - curr
    if isa(d.elem[idx,idy],Number) || isa(d.elem[idx,idy],Variable) || isa(d.elem[idx,idy],AffExpr) || isa(d.elem[idx,idy],QuadExpr)
        return d.elem[idx,idy]
    else
        return getindex(d.elem[idx,idy],offx,offy)
    end
end

getValue(b::BlockMatrix) = hvcat(size(b.elem), map(getValue, b.elem)'...)

function getnames(c::BlockMatrix,d::Set)
    for el in c.elem
        if isa(el,SDPVar) || isa(el,MatrixVar)
            push!(d,el.m.sdpdata.varname[el.index])
        elseif isa(el, BlockMatrix) || isa(el, MatrixExpr)
            getnames(el,d)
        end
    end
    return d
end

###############################################################################
# Matrix Expression class
# Expressions of the form ΣᵢAᵢ*Xᵢ*Bᵢ+C for matrices A, B, C, and where
# X may be a block matrix or a matrix variable. 
type MatrixExpr{T<:Number} <: SDPMatrix
    elem::Vector
    pre::Array
    post::Array
    constant::AbstractArray{T,2}
end

MatrixExpr(n::Int) = MatrixExpr({},{},{},spzeros(n,n))

MatrixExpr(arr::BlockMatrix) = MatrixExpr(concat(arr), concat(𝕀), concat(𝕀), spzeros(size(arr)...))
function MatrixExpr(arr::Array)
    bl = BlockMatrix(arr)
    MatrixExpr(concat(bl), concat(𝕀), concat(𝕀), spzeros(size(bl)...))
end

# really ugly hack to coerce vectors to matrix...
MatrixExpr{T<:Number}(elem::Vector,pre::Array,post::Array,constant::AbstractArray{T}) = 
    MatrixExpr(elem,pre,post,transpose(transpose(constant)))

Base.size(d::MatrixExpr) = (size(d.constant,1),size(d.constant,2))
Base.size(d::MatrixExpr, slice::Int64) = (0 <= slice <= 2) ? size(d)[slice] : 1
Base.ndims(d::MatrixExpr) = 2
Base.eye(d::MatrixExpr)  = eye(size(d)...)

Base.convert(::Type{MatrixExpr}, v::SDPVar)    = MatrixExpr(concat(v), concat(𝕀), concat(𝕀), spzeros(v.dim,v.dim))
Base.convert(::Type{MatrixExpr}, v::MatrixVar) = MatrixExpr(concat(v), concat(𝕀), concat(𝕀), spzeros(v.dim...))

Base.transpose(d::MatrixExpr)  = MatrixExpr(map(transpose, d.elem), map(transpose, d.post), map(transpose, d.pre), transpose(d.constant))
Base.ctranspose(d::MatrixExpr) = MatrixExpr(map(transpose, d.elem), map(transpose, d.post), map(transpose, d.pre), transpose(d.constant))

Base.isequal(a::MatrixExpr, b::MatrixExpr) = isequal(a.elem,b.elem) && isequal(a.pre,b.pre) && isequal(a.post,b.post) && isequal(a.constant,b.constant)

Base.trace(c::MatrixExpr) = MatrixFuncVar(c, :trace)
Base.norm(c::MatrixExpr)  = MatrixFuncVar(c, :norm)
Base.sum(c::MatrixExpr)   = MatrixFuncVar(c, :sum)

Base.issym(d::MatrixExpr) = issym(d.constant) && all(issym, d.elem)

Base.copy(d::MatrixExpr) = MatrixExpr(copy(d.elem), copy(d.pre), copy(d.post), copy(d.constant))

function Base.getindex(d::MatrixExpr, x::Int64)
    m = size(d,1)
    idx,idy = rem(x-1,m)+1, div(x-1,m)+1
    return getindex(d, idx, idy)
end

function Base.getindex(d::MatrixExpr, x::Int64, y::Int64)
    m,n = size(d)
    refer = ScalarExpr(d.constant[x,y])
    for it in 1:length(d.elem)
        !isa(d.pre[it],UniformScaling) && !isa(d.post[it],UniformScaling) && error("Expression not linear in matrix variables")
        if !isa(d.post[it],UniformScaling)
            for k in 1:size(d.elem[it],2)
                refer += d.pre[it].λ*d.elem[it][x,k]*d.post[it][k,y]
            end
            # refer += mapreduce(k->d.pre[it].λ*d.elem[it][x,k]*d.post[it][k,y], +, 1:size(d.elem[it],2))
        elseif !isa(d.pre[it],UniformScaling)
            for k in 1:size(d.elem[it],1)
                refer += d.pre[it][x,k]*d.elem[it][k,y]*d.post[it].λ
            end
            # refer += mapreduce(k->d.pre[it][x,k]*d.elem[it][k,y]*d.post[it].λ, +, 1:size(d.elem[it],1))
        else
            refer += d.pre[it].λ * d.elem[it][x,y] * d.post[it].λ
        end
    end
    return refer
end

function Base.getindex(d::MatrixExpr, x::Int, y::Range{Int})
    sx,sy = size(d)
    arr = Array(AffExpr, 1, length(y))
    for (it,val) in enumerate(y)
        arr[it] = d[x,val]
    end
    return arr
end

function Base.getindex(d::MatrixExpr, x::Range{Int}, y::Int)
    sx,sy = size(d)
    arr = Array(AffExpr, length(x))
    for (it,val) in enumerate(x)
        arr[it] = d[val,y]
    end
    return arr
end

# helper for getValue(c::MatrixExpr)
getValue{T<:Number}(c::AbstractArray{T,2}) = c

getValue(c::MatrixExpr) =
    mapreduce(it->c.pre[it]*getValue(c.elem[it])*c.post[it], +, 1:length(c.elem)) + c.constant

function getnames(c::MatrixExpr,d::Set)
    for el in c.elem
        if isa(el,SDPVar) || isa(el,MatrixVar)
            push!(d,el.m.sdpdata.varname[el.index])
        elseif isa(el, BlockMatrix)
            getnames(el,d)
        end
    end
    return d
end

Base.show(io::IO,d::MatrixExpr)  = print(io, "Matrix expression")
function Base.print(io::IO,d::MatrixExpr)
    n = getnames(d,Set())
    str = join([chomp(string(v)) for v in n], ", ")
    println(io, string("Matrix expression in ", str))
end

###############################################################################
# (Linear function) of a matrix expression
# Represents a linear function acting on a matrix expression of the form
# 𝒮₊ⁿ→ℝ. Current types include a trace operator or element reference.
type MatrixFuncVar
    expr::MatrixExpr
    func::Symbol
end

Base.isequal(x::MatrixFuncVar,y::MatrixFuncVar) = isequal(x.expr,y.expr) && isequal(x.func,y.func)

Base.size(::MatrixFuncVar) = 1
Base.size(::MatrixFuncVar,::Int) = 1

function getValue(v::MatrixFuncVar)
    if v.func == :trace
        return trace(getValue(v.expr))
    elseif v.func == :ref
        (i,j,val) = findnz(v.expr.pre[1])
        if length(i) == 2
            return getValue(MatrixExpr(v.expr.elem,concat(𝕀),concat(𝕀),spzeros(size(v.expr)...)))[i[1],j[1]] * 2val[1]
        elseif length(i) == 1
            return getValue(MatrixExpr(v.expr.elem,concat(𝕀),concat(𝕀),spzeros(size(v.expr)...)))[i[1],j[1]] * val[1]
        else
            error("Shouldn't reach here")
        end
        # return getValue(MatrixExpr(v.expr.elem,concat(𝕀),concat(𝕀),spzeros(size(v.expr)...)))[findn(v.expr.pre[1])...][1] # this is real hacky and ugly
    end
end

setObjective(m::Model, sense::Symbol, c::MatrixFuncVar) = setObjective(m, sense, convert(ScalarExpr,c))

getnames(c::MatrixFuncVar,d::Set)  = getnames(c.expr,d)

###############################################################################
# Norm Expression class
# Stores a normed expression of the form ‖X‖ for some norm ‖·‖, matrix
# expression X. 
type NormExpr
    vars::Vector{Union(AffExpr,MatrixExpr)}
    coeffs::Vector{Float64}
    form::Vector{Symbol}
end

Base.isequal(x::NormExpr,y::NormExpr) = isequal(x.vars,y.vars) && isequal(x.coeffs,y.coeffs) && isequal(x.form,y.form)

NormExpr() = NormExpr(Union(AffExpr,MatrixExpr)[], Float64[], Symbol[])
NormExpr(form::Symbol) = NormExpr(Union(AffExpr,MatrixExpr)[], [form])

Base.abs(d::Variable) = NormExpr(concat(convert(AffExpr,d)), [1.0], [:abs])
Base.abs(d::AffExpr)  = NormExpr(concat(d), [1.0], [:abs], AffExpr())

Base.norm(d::MatrixVar)  = NormExpr(concat(convert(MatrixExpr,d)), [1.0], [:norm2])
Base.norm(d::MatrixExpr) = NormExpr(concat(d), [1.0], [:norm2])

Base.vecnorm(d::MatrixVar)  = NormExpr(concat(convert(MatrixExpr,d)), [1.0], [:normfrob])
Base.vecnorm(d::MatrixExpr) = NormExpr(concat(d), [1.0], [:normfrob])

Base.copy(d::NormExpr) = NormExpr(copy(d.vars), copy(d.coeffs), copy(d.form))

function getValue(d::NormExpr)
    ret = 0.0
    for it in 1:length(d.vars)
        if d.form[it] == :abs
            ret += coeffs[it] * abs(getValue(d.vars[it]))
        elseif d.form[it] == :norm2
            ret += coeffs[it] * norm(getValue(d.vars[it]))
        elseif d.form[it] == :normfrob
            ret += coeffs[it] * vecnorm(getValue(d.vars[it]))
        end
    end
    return ret
end

setObjective(m::Model, sense::Symbol, c::NormExpr) = setObjective(m, sense, convert(ScalarExpr,c))

getnames(c::NormExpr,d::Set) = map(x->getnames(x,d), c.vars)

function exprToStr(c::NormExpr)
    d = Set()
    for var in c.vars
        getnames(var,d)
    end
    str = join([chomp(string(v)) for v in d], ", ")
    return string("Norm expression in ", str)
end

Base.show(io::IO, d::NormExpr)  = print(io, "Norm expression")
Base.print(io::IO, d::NormExpr) = print(io, exprToStr(d))

###############################################################################
# Scalar Expression class
# Combination of AffExpr, MatrixFuncExpr, and NormExpr.
type ScalarExpr
    matvars::Vector{MatrixFuncVar}
    aff::AffExpr
    normexpr::NormExpr
end

Base.isequal(x::ScalarExpr,y::ScalarExpr) = isequal(x.matvars,y.matvars) && isequal(x.aff,y.aff) && isequal(x.normexpr,y.normexpr)

ScalarExpr() = ScalarExpr(MatrixFuncVar[], AffExpr(), NormExpr())
ScalarExpr(v::Number) = ScalarExpr(MatrixFuncVar[], AffExpr(v), NormExpr())

Base.convert(::Type{ScalarExpr}, v::MatrixFuncVar) = ScalarExpr(concat(v), AffExpr(), NormExpr())
Base.convert(::Type{ScalarExpr}, v::AffExpr) = ScalarExpr(MatrixFuncVar[], v, NormExpr())
Base.convert(::Type{ScalarExpr}, v::NormExpr) = ScalarExpr(MatrixFuncVar[], AffExpr(), v)

function Base.convert(::Type{ScalarExpr}, v::MatrixExpr)
    size(v) == (1,1) || error("Cannot coerce matrix expression to scalar expression")
    return v[1]
end

Base.size(::ScalarExpr) = 1
Base.size(::ScalarExpr,::Int) = 1

getValue(v::ScalarExpr) = getValue(v.aff) + mapreduce(getValue, +, v.matvars) + getValue(v.normexpr)

setObjective(m::Model, sense::Symbol, c::MatrixExpr) = setObjective(m, sense, convert(ScalarExpr,c))
function setObjective(m::Model, sense::Symbol, c::ScalarExpr)
    initSDP(m)
    setObjectiveSense(m, sense)
    m.sdpdata.sdpobj = c
end

function getnames(c::ScalarExpr,d::Set)
    map(x->getnames(x,d), c.matvars)
    getnames(c.aff,d)
    getnames(c.normexpr,d)
    return d
end
getnames(c::ScalarExpr) = getnames(c,Set())
getnames(c::Variable,d::Set) = (push!(d,c.m.colNames[c.col]))
function getnames(c::Variable)
    d = Set()
    push!(d, c.m.colNames[c.col])
end
getnames(a::AffExpr,d::Set) = map(x->getnames(x,d), a.vars)
function getnames(a::AffExpr)
    d = Set()
    getnames(a,d)
    return d
end

function exprToStr(a::ScalarExpr)
    d = Set()
    for it in a.matvars
        getnames(it,d)
    end
    getnames(a.aff,d)
    str = join([chomp(string(v)) for v in d], ", ")
    return string("Scalar expression in ", str)
end

Base.show(io::IO, c::ScalarExpr)  = print(io, "Scalar expression")
Base.print(io::IO, c::ScalarExpr) = print(io, exprToStr(c))

###############################################################################
# Dual Expression class
# Expressions of the form ΣᵢyᵢAᵢ + C, where Aᵢ and C are (symmetric) matrices
# and yᵢ are scalar variables. Used in dual SDP constraints.
typealias DualExpr JuMP.GenericAffExpr{AbstractArray{Float64,2},Union(Variable,MatrixFuncVar)}

Base.isequal(x::DualExpr,y::DualExpr) = isequal(x.vars,y.vars) && isequal(x.coeffs,y.coeffs) && (x.constant == y.constant)

DualExpr(n::Integer) = DualExpr({},AbstractArray[],spzeros(n,n))

Base.size(d::DualExpr) = size(d.constant)
Base.size(d::DualExpr, slice::Int64) = size(d.constant,slice)

Base.issym(d::DualExpr) =  issym(d.constant) && all(issym, d.coeffs)

function Base.getindex(d::DualExpr, x::Int, y::Int)
    m,n = size(d)
    ret = ScalarExpr(d.constant[x,y])
    for it in 1:m
        ret += d.vars[it]*d.coeffs[it][x,y]
    end
    return ret
    #return mapreduce(it->d.vars[it]*d.coeffs[it][x,y], +, 1:m) + d.constant[x,y]
end

getValue(d::DualExpr) = mapreduce(it->d.coeffs[it]*getValue(d.vars[it]), +, 1:length(d.vars)) + d.constant

function getnames(c::DualExpr)
    d = Set()
    for v in c.vars
        if isa(v,Variable)
            push!(d, v.m.colNames[v.col])
        elseif isa(v,MatrixFuncVar)
            getnames(v.expr,d)
        end
    end
    return d
end

Base.show(io::IO, c::DualExpr)  = print(io, "Dual expression in ", join(getnames(c),", "))
Base.print(io::IO, c::DualExpr) = println(io, "Dual expression in ", join(getnames(c),", "))

###############################################################################
# Primal Constraint class
# Stores a constraint of the type ub ≤ X ≤ lb, where X is a matrix function
# expression.

typealias PrimalConstraint JuMP.GenericRangeConstraint{ScalarExpr}

function addConstraint(m::Model, c::PrimalConstraint)
    initSDP(m)
    push!(m.sdpdata.primalconstr,c)
    return ConstraintRef{PrimalConstraint}(m,length(m.sdpdata.primalconstr))
end

function conToStr(c::PrimalConstraint)
    d = Set()
    getnames(c.terms,d)
    str = join([chomp(string(v)) for v in d], ", ")
    return string("Primal constraint in ", str) 
end

Base.show(io::IO, c::ConstraintRef{PrimalConstraint})  = print(io, conToStr(c.m.sdpdata.primalconstr[c.idx]))
Base.print(io::IO, c::ConstraintRef{PrimalConstraint}) = print(io, conToStr(c.m.sdpdata.primalconstr[c.idx]))

###############################################################################
# Dual Constraint class
# Stores a constraint of the type X ? 0, where X is a dual expression and ? is
# inequality or equality.
type DualConstraint <: JuMPConstraint
    terms::DualExpr
    sense::Symbol
end

function addConstraint(m::Model, c::DualConstraint)
    initSDP(m)
    push!(m.sdpdata.dualconstr,c)
    return ConstraintRef{DualConstraint}(m,length(m.sdpdata.dualconstr))
end

conToStr(c::DualConstraint) = string("Dual constraint in ", join(getnames(c.terms),", ")) 

Base.show(io::IO, c::ConstraintRef{DualConstraint})  = print(io, conToStr(c.m.sdpdata.dualconstr[c.idx]))
Base.print(io::IO, c::ConstraintRef{DualConstraint}) = print(io, conToStr(c.m.sdpdata.dualconstr[c.idx]))

###############################################################################
# Matrix Constraint class
# Stores a constraint of the type X ? 0, where X is a matrix expression and ? 
# is a a semidefinite inequality (>=, <=, or ==) or an entrywise inequality
# (.>= or .<=).
type MatrixConstraint <: JuMPConstraint
    terms::MatrixExpr
    sense::Symbol
end

function addConstraint(m::Model, c::MatrixConstraint)
    initSDP(m)
    push!(m.sdpdata.matrixconstr,c)
    return ConstraintRef{MatrixConstraint}(m,length(m.sdpdata.matrixconstr))
end

function conToStr(c::MatrixConstraint)
    d = Set()
    getnames(c.terms,d)
    str = join([chomp(string(v)) for v in d], ", ")
    return string("SDP matrix constraint in ", str) 
end

Base.show(io::IO, c::ConstraintRef{MatrixConstraint})  = print(io, conToStr(c.m.sdpdata.matrixconstr[c.idx]))
Base.print(io::IO, c::ConstraintRef{MatrixConstraint}) = print(io, conToStr(c.m.sdpdata.matrixconstr[c.idx]))
