# compute hessian sparsity pattern given expression graph

# Idea: detect partially separable structure
# e.g., f(x) = sum f_i(x)

function compute_hessian_sparsity(s::SymbolicOutput)
    code = quote end
    compute_hessian_sparsity(s.tree, true, code)
    # compile a function:
    finit = Expr(:block)
    for i in 1:length(s.inputnames)
        push!( finit.args, :( $(s.inputnames[i]) = $(s.inputvals[i]) ))
    end
    fexpr = quote
        function tmp(colorlist)
            $finit
            $code
        end
    end
    
    clist = Dict()
    eval(fexpr)(clist)
    edgelist = {}
    for color in values(clist)
        # for every pair:
        for x in color
            for y in color
                # don't worry about duplicate indices
                push!(edgelist,(x,y))
            end
        end
    end
    return edgelist

end

function compute_hessian_sparsity_IJ(s::SymbolicOutput)
    edgelist = compute_hessian_sparsity(s)
    I = Array(Int,0)
    J = Array(Int,0)
    for k in 1:length(edgelist)
        x,y = edgelist[k]
        i = getindex(x)
        j = getindex(y)
        if j > i
            continue # ignore upper triangle
        else
            # convert to canonical
            push!(I,s.maptocanonical[i])
            push!(J,s.maptocanonical[j])
        end
    end
    return I,J
end

function compute_hessian_sparsity(x::ExprNode, linear_so_far, expr_out)
    if isexpr(x.ex, :call)
        # is this a linear operator?
        if linear_so_far == 1 
            if x.ex.args[1] == :(+) || x.ex.args[1] == :(-) || (x.ex.args[1] == :(*) && sum([isa(t,ExprNode) for t in x.ex.args[2:end]]) <= 1) 
                # we're still at the top
                # give the same color to the children
                for i in 2:length(x.ex.args)
                    compute_hessian_sparsity(x.ex.args[i], true, expr_out)
                end
            else
                # this is a new f_i, make a new color
                code = quote
                    mycolor = gensym()
                    colorlist[mycolor] = Set()
                end
                for i in 2:length(x.ex.args)
                    compute_hessian_sparsity(x.ex.args[i], false, code)
                end
                push!(expr_out.args, quote let 
                            $code 
                        end end)
            end
        else
            # give the same color to the children
            for i in 2:length(x.ex.args)
                compute_hessian_sparsity(x.ex.args[i], false, expr_out)
            end
        end
    elseif isexpr(x.ex, :curly)
        @assert x.ex.args[1] == :sum
        # we are a linear operator
        code = quote end
        compute_hessian_sparsity(x.ex.args[2], linear_so_far, code)
        for level in length(x.ex.args):-1:3
            code = Expr(:for, x.ex.args[level],code)
        end
        push!(expr_out.args, code)
    else
        # we must be at an input expression
        # add this placeholder to the set of nodes with this color
        # (placeholders can have multiple colors if they appear multiple times in the graph)
        if !linear_so_far
            push!(expr_out.args, :( 
                if isa($(x.ex),Placeholder)
                    push!(colorlist[mycolor], $(x.ex))
                end))
        end
    end
    
end

compute_hessian_sparsity(x, linear_so_var, expr_out) = nothing

export compute_hessian_sparsity_IJ

# returns a function that computes a dense hessian matrix
function gen_hessian_dense(s::SymbolicOutput)
    gradexpr = genfgrad(s)
    fgrad = eval(gradexpr)
    function hessian{T}(x::Vector{T}, output_matrix::Matrix{T})
        @assert length(x) == size(output_matrix,1) == size(output_matrix,2)
        fill!(output_matrix, zero(T))
        dualvec = Array(Dual{T}, length(x))
        dualout = Array(Dual{T}, length(x))
        for i in 1:length(x)
            dualvec[i] = Dual(x[i], zero(T))
        end
        for i in 1:length(x)
            dualvec[i] = Dual(x[i], one(T))
            fill!(dualout, dual(zero(T)))
            fgrad(dualvec, IdentityArray(), dualout, IdentityArray())
            dualvec[i] = Dual(x[i], zero(T))
            for j in 1:length(x)
                output_matrix[j,i] = epsilon(dualout[j])
            end
        end
        return output_matrix
    end
    function hwrapper{T}(x::Vector{T})
        H = Array(T,length(x),length(x))
        hessian(x,H)
    end
end

# returns a function that computes a hessian-matrix product
function gen_hessian_matmat(s::SymbolicOutput)
    gradexpr = genfgrad(s)
    fgrad = eval(gradexpr)
    function hessian_matmat!{T}(S, x::Vector{T})
        @assert length(x) == size(S,1)
        dualvec = Array(Dual{T}, length(x))
        dualout = Array(Dual{T}, length(x))
        for k in 1:size(S,2)
            for i in 1:length(x)
                dualvec[i] = Dual(x[i], S[i,k])
            end
            fill!(dualout, dual(zero(T)))
            fgrad(dualvec, IdentityArray(), dualout, IdentityArray())
            for i in 1:length(x)
                S[i,k] = epsilon(dualout[i])
            end
        end
        return S
    end
end

function gen_hessian_matmat_parametric(s::SymbolicOutput)
    fgrad = genfgrad_parametric(s)
    # S uses canonical indices
    function hessian_matmat_p!{T}(S, x::Vector{T}, inputvals, fromcanonical, tocanonical)
        N = size(S,1)
        @assert length(x) >= N
        dualvec = Array(Dual{T}, N)
        dualout = Array(Dual{T}, N)
        for k in 1:size(S,2)
            for i in 1:N
                dualvec[i] = Dual(x[fromcanonical[i]], S[i,k])
            end
            fill!(dualout, dual(zero(T)))
            fgrad(dualvec, tocanonical, dualout, tocanonical,inputvals)
            for i in 1:N
                S[i,k] = epsilon(dualout[i])
            end
        end
        return S
    end
end


# returns sparse matrix containing hessian structure, and
# function to evaluate the hessian into a given sparse matrix
function gen_hessian_sparse_mat(s::SymbolicOutput)
    fgrad = eval(genfgrad(s))
    I,J = compute_hessian_sparsity_IJ(s)
    structuremat = sparse(I,J,ones(length(I)))
    fill!(structuremat.nzval, 0.0)

    function evalhessian{T}(x::Vector{T}, output_matrix::SparseMatrixCSC{T,Int})
        dualvec = Array(Dual{T}, length(x))
        dualout = Array(Dual{T}, length(x))
        for i in 1:length(x)
            dualvec[i] = Dual(x[i], zero(T))
        end
        rowval = output_matrix.rowval
        nzval = output_matrix.nzval
        for i in 1:length(x)
            dualvec[i] = Dual(x[i], one(T))
            fill!(dualout, dual(zero(T)))
            fgrad(dualvec, IdentityArray(), dualout, IdentityArray())
            dualvec[i] = Dual(x[i], zero(T))
            for idx in output_matrix.colptr[i]:(output_matrix.colptr[i+1]-1)
                nzval[idx] = epsilon(dualout[rowval[idx]])
            end
        end
        return output_matrix
    end
    return structuremat, evalhessian
end

function gen_hessian_sparse_raw(s::SymbolicOutput)
    fgrad = eval(genfgrad(s))

    function evalhessian{T}(placevalues::Vector{T}, output_vec::Vector{T}, nzstruct::SparseMatrixCSC{Int,Int})
        dualvec = Array(Dual{T}, length(placevalues))
        dualout = Array(Dual{T}, length(placevalues))
        for i in 1:length(x)
            dualvec[i] = Dual(placevalues[i], zero(T))
        end
        rowval = nzstruct.rowval
        nzidx = nzstruct.nzval
        for i in 1:length(placevalues)
            dualvec[i] = Dual(placevalues[i], one(T))
            fill!(dualout, dual(zero(T)))
            fgrad(dualvec, IdentityArray(), dualout, IdentityArray())
            dualvec[i] = Dual(placevalues[i], zero(T))
            for idx in nzstruct.colptr[i]:(nzstruct.colptr[i+1]-1)
                output_vec[nzidx[idx]] += epsilon(dualout[rowval[idx]])
            end
        end
        return nothing
    end
    return evalhessian
end

export gen_hessian_dense, gen_hessian_sparse_mat, gen_hessian_sparse_raw, gen_hessian_matmat

            



