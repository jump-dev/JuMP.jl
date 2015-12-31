
# convert from Julia expression into NodeData form


function expr_to_nodedata(ex::Expr)
    nd = NodeData[]
    values = Float64[]
    expr_to_nodedata(ex,nd,values,-1,-1)
    return nd,values
end

function expr_to_nodedata(ex::Expr,nd::Vector{NodeData},values::Vector{Float64},parentid,whichchild)

    myid = length(nd) + 1
    if isexpr(ex,:call)
        op = ex.args[1]
        if length(ex.args) == 2
            push!(nd,NodeData(CALLUNIVAR, univariate_operator_to_id[op], parentid, whichchild))
        else
            push!(nd,NodeData(CALL, operator_to_id[op], parentid,whichchild))
        end
        for k in 2:length(ex.args)
            expr_to_nodedata(ex.args[k],nd,values,myid,k-1)
        end
    elseif isexpr(ex, :ref)
        @assert ex.args[1] == :x
        push!(nd,NodeData(VARIABLE, ex.args[2], parentid, whichchild))
    elseif isexpr(ex, :comparison)
        op = ex.args[2]
        opid = comparison_operator_to_id[op]
        for k in 2:2:length(ex.args)-1
            @assert ex.args[k] == op
        end
        push!(nd, NodeData(COMPARISON, opid, parentid, whichchild))
        for k in 1:2:length(ex.args)
            expr_to_nodedata(ex.args[k],nd,values,myid,div(k+1,2))
        end
    elseif isexpr(ex,:&&) || isexpr(ex,:||)
        @assert length(ex.args) == 2
        op = ex.head
        opid = logic_operator_to_id[op]
        push!(nd, NodeData(LOGIC, opid, parentid, whichchild))
        expr_to_nodedata(ex.args[1],nd,values,myid,1)
        expr_to_nodedata(ex.args[2],nd,values,myid,2)
    else
        error("Unrecognized expression $ex: $(ex.head)")
    end
    nothing
end

function expr_to_nodedata(ex::Number,nd::Vector{NodeData},values::Vector{Float64},parentid,whichchild)
    valueidx = length(values)+1
    push!(values,ex)
    push!(nd, NodeData(VALUE, valueidx, parentid,whichchild))
    nothing
end

export expr_to_nodedata

# (i,j) nonzero means there's an edge *from* j to *i*
# since we get a column-oriented matrix, this gives us a fast way to look up the
# edges leaving any node (i.e., the children)
function adjmat(nd::Vector{NodeData})
    len = length(nd)
    I = Array(Int,len)
    J = Array(Int,len)
    realnz = 0
    for nz in 1:len
        par = nd[nz].parent
        par > 0 || continue
        realnz += 1
        I[realnz] = nz
        J[realnz] = par
    end
    resize!(I,realnz)
    resize!(J,realnz)
    return sparse(I,J,ones(Bool,realnz),len,len)
end

export adjmat
