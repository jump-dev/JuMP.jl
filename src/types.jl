

type ExprNode
    ex
    parents
    value # value of this expression in the tree
    deriv # derivative wrt this expression (when computed)
    linear_so_far::Bool # true if the result of the expression is known to be linear with respect to the value of this node
end

ExprNode(ex, parents, linear_so_far) = ExprNode(ex, parents, nothing, nothing, linear_so_far)

abstract Placeholder

immutable BasicPlaceholder <: Placeholder
    idx::Int
end

getplaceindex(x::BasicPlaceholder) = x.idx

function placeholders(n::Int)
    v = Array(BasicPlaceholder, n)
    for i in 1:n
        v[i] = BasicPlaceholder(i)
    end
    return v
end

export Placeholder, BasicPlaceholder, placeholders

function Base.dump(io::IO, x::ExprNode)
    dump(io, x.ex)
    dump(io, x.value)
    dump(io, x.deriv)
end

function Base.show(io::IO, x::ExprNode)
    show(io, x.ex)
end


getvalue(x::ExprNode) = x.value
getvalue(x) = x # for constants

function cleargraph(x::ExprNode)
    x.value = nothing
    x.deriv = nothing
    if isa(x.ex,Expr)
        for ex in x.ex.args
            cleargraph(ex)
        end
    end
end

cleargraph(x) = nothing

# set over finite indices with O(1) lookup, addition, and deletion.
# similar to IndexedVector in JuMP
type IndexedSet
    idx::Vector{Int}
    set::BitArray{1} # true if element is in set
    nel::Int # number of elements in set
end

IndexedSet(n::Integer) = IndexedSet(zeros(Int,n),falses(n),0)

function Base.push!(s::IndexedSet,i::Integer)
    s.set[i] && return # already there
    s.set[i] = true
    s.nel += 1
    s.idx[s.nel] = i
    return
end
Base.length(s::IndexedSet) = s.nel
function Base.empty!(s::IndexedSet)
    for i in 1:s.nel
        s.set[s.idx[i]] = false
    end
    s.nel = 0
    return
end

function Base.unique(s::IndexedSet)
    arr = Array(Int,length(s))
    for i in 1:length(s)
        arr[i] = s.idx[i]
    end
    return arr
end
