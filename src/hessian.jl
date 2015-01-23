# compute hessian sparsity pattern given expression graph

# convert to lower triangular indices
normalize(i,j) = (j > i) ? (j,i) : (i,j)
normalize(e) = normalize(e...)

# Idea: detect partially separable structure
# e.g., f(x) = sum f_i(x)

function compute_hessian_sparsity_IJ_parametric(s::SymbolicOutput)

    code = quote end
    compute_hessian_sparsity(s.tree, true, code)
    # compile a function:
    fexpr = quote
        local _SPARSITY_GEN_
        function _SPARSITY_GEN_(edgelist__)
            mycolor__ = Set{Int}()
            $code
            return
        end
        _SPARSITY_GEN_
    end
    # add arguments for inputnames -- local data
    for i in 1:length(s.inputnames)
        push!(fexpr.args[4].args[1].args,s.inputnames[i])
    end

    f = eval(fexpr)

    return (x::SymbolicOutput)-> (edgelist__ = Array((Int,Int),0); f(edgelist__,x.inputvals...); edgelist_to_IJ(edgelist__,x))

end

compute_hessian_sparsity_IJ(s::SymbolicOutput) = (compute_hessian_sparsity_IJ_parametric(s))(s)

function edgelist_to_IJ(edgelist,s::SymbolicOutput)
    I = Array(Int,0)
    J = Array(Int,0)
    sizehint(I,div(length(edgelist),2))
    sizehint(J,div(length(edgelist),2))
    for k in 1:length(edgelist)
        i,j = edgelist[k]
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
    nonlinear_cleanup = quote
        # for every pair:
        for x1___ in mycolor__
            for x2___ in mycolor__
                # don't worry about duplicate indices
                push!(edgelist__,(x1___,x2___))
            end
        end
    end
    if isexpr(x.ex, :call)
        # is this a linear operator?
        if linear_so_far == 1
            # pre-generate the code for the two cases

            code_linear = quote end
            # we're still at the top
            # give the same color to the children
            for i in 2:length(x.ex.args)
                compute_hessian_sparsity(x.ex.args[i], true, code_linear)
            end

            # this is a new f_i, make a new color
            code_nonlinear = quote
                empty!(mycolor__)
            end
            for i in 2:length(x.ex.args)
                compute_hessian_sparsity(x.ex.args[i], false, code_nonlinear)
            end


            if x.ex.args[1] == :(+) || x.ex.args[1] == :(-) || x.ex.args[1] == :ifelse
                push!(expr_out.args, code_linear)
            elseif x.ex.args[1] == :(*)
                # if at most one of the multiplicands is a "complex" expression,
                # try to detect if all others are constants.
                iscomplexexpr(t) = isa(t,ExprNode) && (isexpr(t.ex,:call) || isexpr(t.ex,:curly))
                num_complex = mapreduce(iscomplexexpr, +, x.ex.args[2:end])
                if num_complex <= 1
                    conditions = quote
                        num_nonconstant = $num_complex
                    end
                    for t in x.ex.args[2:end]
                        iscomplexexpr(t) && continue
                        if isa(t,ExprNode)
                            push!(conditions.args,
                            :(num_nonconstant += isa($(t.ex),Placeholder)))
                        end # else, definitely a constant
                    end
                    push!(expr_out.args, conditions)
                    push!(expr_out.args, Expr(:if, :(num_nonconstant <= 1), code_linear,
                        :($code_nonlinear; $nonlinear_cleanup)))
                else
                    # too hard to detect constants, just say nonlinear
                    push!(expr_out.args, quote
                            $code_nonlinear
                            $nonlinear_cleanup
                        end)
                end

            else

                push!(expr_out.args, quote
                            $code_nonlinear
                            $nonlinear_cleanup
                        end)
            end
        else
            # give the same color to the children
            for i in 2:length(x.ex.args)
                compute_hessian_sparsity(x.ex.args[i], false, expr_out)
            end
        end
    elseif isexpr(x.ex, :curly)
        if issum(x.ex.args[1]) || !linear_so_far
            code = quote end
            compute_hessian_sparsity(curlyexpr(x.ex), linear_so_far, code)
            push!(expr_out.args, gencurlyloop(x.ex,code))
        else
            @assert isprod(x.ex.args[1])
            # new color
            code = quote end
            compute_hessian_sparsity(curlyexpr(x.ex), false, code)
            push!(expr_out.args, 
                quote let
                        empty!(mycolor__)
                        $(gencurlyloop(x.ex, code))
                        $nonlinear_cleanup
                end end)
        end
    elseif isexpr(x.ex, :comparison) || isexpr(x.ex, :&&) || isexpr(x.ex, :||)
        # nothing to do here
    else
        # we must be at an input expression
        # add this placeholder to the set of nodes with this color
        # (placeholders can have multiple colors if they appear multiple times in the graph)
        if !linear_so_far
            push!(expr_out.args, :( 
                if isa($(x.ex),Placeholder)
                    push!(mycolor__, getplaceindex($(x.ex)))
                end))
        end
    end
    
end

compute_hessian_sparsity(x, linear_so_far, expr_out) = nothing

export compute_hessian_sparsity_IJ

# returns a function that computes a dense hessian matrix
function gen_hessian_dense(s::SymbolicOutput)
    fgrad = genfgrad_parametric(s)
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
            fgrad(dualvec, IdentityArray(), dualout, IdentityArray(), s.inputvals...)
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




function gen_hessian_matmat_parametric(s::SymbolicOutput, fgrad = genfgrad_parametric(s))
    hexpr = quote
        local _HESS_MATMAT_
        function _HESS_MATMAT_{T,Q}(S, x::Vector{T}, dualvec::Vector{Dual{T}}, dualout::Vector{Dual{T}}, inputvals::Q, fromcanonical)
            # S uses canonical indices
            N = size(S,1)
            @assert length(x) >= N
            for k in 1:size(S,2)
                for i in 1:N
                    r = fromcanonical[i]
                    dualvec[r] = Dual(x[r], S[i,k])
                    dualout[r] = zero(Dual{T})
                end
                $(fgrad)(dualvec, IdentityArray(), dualout, IdentityArray(),inputvals...)
                for i in 1:N
                    S[i,k] = epsilon(dualout[fromcanonical[i]])
                end
            end
            #return S
        end
        _HESS_MATMAT_
    end

    return eval(hexpr)
end


# returns sparse matrix containing hessian structure, and
# function to evaluate the hessian into a given sparse matrix
function gen_hessian_sparse_mat(s::SymbolicOutput)
    fgrad = genfgrad_parametric(s)
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
            fgrad(dualvec, IdentityArray(), dualout, IdentityArray(), s.inputvals...)
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
    fgrad = genfgrad_parametric(s)

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
            fgrad(dualvec, IdentityArray(), dualout, IdentityArray(), s.inputvals...)
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

            



