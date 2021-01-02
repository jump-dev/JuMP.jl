
# convert from Julia expression into NodeData form

function expr_to_nodedata(
    ex::Expr,
    r::UserOperatorRegistry = UserOperatorRegistry(),
)
    nd = NodeData[]
    values = Float64[]
    expr_to_nodedata(ex, nd, values, -1, r)
    return nd, values
end

function expr_to_nodedata(
    ex::Expr,
    nd::Vector{NodeData},
    values::Vector{Float64},
    parentid,
    r::UserOperatorRegistry,
)
    myid = length(nd) + 1
    if isexpr(ex, :call)
        op = ex.args[1]
        if length(ex.args) == 2
            id =
                haskey(univariate_operator_to_id, op) ?
                univariate_operator_to_id[op] :
                r.univariate_operator_to_id[op] +
                USER_UNIVAR_OPERATOR_ID_START - 1
            push!(nd, NodeData(CALLUNIVAR, id, parentid))
        elseif op in comparison_operators
            push!(
                nd,
                NodeData(COMPARISON, comparison_operator_to_id[op], parentid),
            )
        else
            id =
                haskey(operator_to_id, op) ? operator_to_id[op] :
                r.multivariate_operator_to_id[op] + USER_OPERATOR_ID_START - 1
            push!(nd, NodeData(CALL, id, parentid))
        end
        for k in 2:length(ex.args)
            expr_to_nodedata(ex.args[k], nd, values, myid, r)
        end
    elseif isexpr(ex, :ref)
        @assert ex.args[1] == :x
        push!(nd, NodeData(VARIABLE, ex.args[2], parentid))
    elseif isexpr(ex, :comparison)
        op = ex.args[2]
        opid = comparison_operator_to_id[op]
        for k in 2:2:length(ex.args)-1
            @assert ex.args[k] == op
        end
        push!(nd, NodeData(COMPARISON, opid, parentid))
        for k in 1:2:length(ex.args)
            expr_to_nodedata(ex.args[k], nd, values, myid, r)
        end
    elseif isexpr(ex, :&&) || isexpr(ex, :||)
        @assert length(ex.args) == 2
        op = ex.head
        opid = logic_operator_to_id[op]
        push!(nd, NodeData(LOGIC, opid, parentid))
        expr_to_nodedata(ex.args[1], nd, values, myid, r)
        expr_to_nodedata(ex.args[2], nd, values, myid, r)
    else
        error("Unrecognized expression $ex: $(ex.head)")
    end
    return nothing
end

function expr_to_nodedata(
    ex::Number,
    nd::Vector{NodeData},
    values::Vector{Float64},
    parentid,
    r::UserOperatorRegistry,
)
    valueidx = length(values) + 1
    push!(values, ex)
    push!(nd, NodeData(VALUE, valueidx, parentid))
    return nothing
end

export expr_to_nodedata

# (i,j) nonzero means there's an edge *from* j to *i*
# since we get a column-oriented matrix, this gives us a fast way to look up the
# edges leaving any node (i.e., the children)
function adjmat(nd::Vector{NodeData})
    len = length(nd)
    I = Vector{Int}(undef, len)
    J = Vector{Int}(undef, len)
    realnz = 0
    for nz in 1:len
        par = nd[nz].parent
        par > 0 || continue
        realnz += 1
        I[realnz] = nz
        J[realnz] = par
    end
    resize!(I, realnz)
    resize!(J, realnz)
    return sparse(I, J, ones(Bool, realnz), len, len)
end

export adjmat
