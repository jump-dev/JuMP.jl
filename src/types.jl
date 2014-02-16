

import Base: dump, show

type ExprNode
    ex
    parents
    value # value of this expression in the tree
    deriv # derivative wrt this expression (when computed)
end

abstract Placeholder

immutable BasicPlaceholder <: Placeholder
    idx::Int
end

getindex(x::BasicPlaceholder) = x.idx

function placeholders(n::Int)
    v = Array(BasicPlaceholder, n)
    for i in 1:n
        v[i] = BasicPlaceholder(i)
    end
    return v
end

export Placeholder, BasicPlaceholder, placeholders

function dump(io::IO, x::ExprNode)
    dump(io, x.ex)
    dump(io, x.value)
    dump(io, x.deriv)
end

function show(io::IO, x::ExprNode)
    show(io, x.ex)
end


getvalue(x::ExprNode) = x.value
getvalue(x) = x # for constants

function cleargraph(x::ExprNode)
    x.value = nothing
    x.deriv = nothing
    if isexpr(x.ex,:call)
        for i in 2:length(x.ex.args)
            cleargraph(x.ex.args[i])
        end
    elseif isexpr(x.ex,:curly)
        @assert x.ex.args[1] == :sum
        cleargraph(x.ex.args[2])
    end
end

cleargraph(x) = nothing
