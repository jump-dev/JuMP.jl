typealias Matrices Union(SDPMatrix, Variable, AffExpr, MatrixFuncVar, ScalarExpr, DualExpr, AbstractArray, Number)
# typealias Matrices Union(SDPMatrix, AbstractArray, Number)

# add these to deal with cases like [1 ones(1,3)]
# They are copied from generic version in Base
vcat(args::Union(AbstractArray,Number)...) = cat(1, args...)
hcat(args::Union(AbstractArray,Number)...) = cat(2, args...)
function hvcat(rows::(Int...), as::Union(AbstractArray,Number)...)
    nbr = length(rows)  # number of block rows
    rs = cell(nbr)
    a = 1
    for i = 1:nbr
        rs[i] = hcat(as[a:a-1+rows[i]]...)
        a += rows[i]
    end
    vcat(rs...)
end

function hcat(args::Matrices...)
    n = length(args)
    tmp = Array(Matrices, 1, n)
    for i in 1:n
        tmp[1,i] = args[i]
    end
    return MatrixExpr(tmp)
end

function vcat(args::Matrices...)
    m = length(args)
    tmp = Array(Matrices, m, 1)
    for i in 1:m
        tmp[i,1] = args[i]
    end
    return MatrixExpr(tmp)
end

function hvcat(dims::(Int64...,), args::Matrices...)
    m, n = dims
    tmp = Array(Matrices, m, n)
    cnt = 1
    for i in 1:m
        for j in 1:n
            tmp[i,j] = args[cnt]
            cnt += 1
        end
    end
    return MatrixExpr(tmp)
end

# internal concatenation function so that vcat/hcat/hvcat can remain user-facing
function concat(X...)
    nargs = length(X)
    isList = map(a->isa(a,Vector) & !issubtype(eltype(a),Number), X)
    dimsX = map((a->isList[a] ? size(X[a],1) : 1), 1:nargs)
    dimC = sum(dimsX)

    types = Set{DataType}()
    for (k,x) in enumerate(X)
        if isList[k]
            elt = eltype(x)
            if isa(elt, UnionType)
                for it in elt.types
                    push!(types, it)
                end
            else
                push!(types, elt)
            end
        else
            push!(types, typeof(x))
        end
    end

    # I'm not sure what the best approach is here, but I think that since this
    # is all internal code, I think this will work fine with ScalarExpr etc.
    # typeC = promote_type(types...)
    typeC = Union(types...)
    C = Array(typeC, dimC)

    range = 1
    for k=1:nargs
        nextrange = range + dimsX[k]
        if dimsX[k] == 1 && !isList[k]
            C[range] = X[k]
        else
            C[range:(nextrange-1)] = X[k]
        end
        range = nextrange
    end
    return C
end
