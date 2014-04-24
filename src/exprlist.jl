# Helps manage a set of expressions
# Tries to share AD across sets of expressions with identical structure
# Useful for computing the Hessian of the lagrangian and related terms


export ExprList, prep_sparse_hessians, eval_g!,
    eval_jac_g!, jac_nz, eval_hess!

type ExprList
    exprs::Vector{SymbolicOutput}
    referenceExpr::Dict
   # matchesCanonical::Vector{Bool}
    valfuncs::Vector{Function}
    gradfuncs::Vector{Function}
    hessfuncs::Vector{Function}
    hessIJ::Vector{Tuple}
end


ExprList() = ExprList(SymbolicOutput[],Dict(), Function[], Function[], Function[], Tuple[])

import Base: push!, getindex
push!(l::ExprList, s::SymbolicOutput) = push!(l.exprs, s)
getindex(l::ExprList, i) = l.exprs[i]

function appendToIJ!(I,J,hI,hJ,x::SymbolicOutput)
    @assert length(hI) == length(hJ)
    mp = x.mapfromcanonical
    for k in 1:length(hI)
        i,j = normalize(mp[hI[k]],mp[hJ[k]])
        push!(I,i)
        push!(J,j)
    end
end

# returns sparsity pattern of combined hessian, with duplicates
function prep_sparse_hessians(l::ExprList)
    @assert length(l.referenceExpr) == 0
    @assert length(l.valfuncs) == 0
    @assert length(l.gradfuncs) == 0
    @assert length(l.hessfuncs) == 0

    N = length(l.exprs)
    sizehint(l.valfuncs, N)
    sizehint(l.gradfuncs, N)
    sizehint(l.hessfuncs, N)

    I = Int[]
    J = Int[]

    for (i,x) in enumerate(l.exprs)
        if !haskey(l.referenceExpr, x.hashval)
            l.referenceExpr[x.hashval] = i
            f = genfval_parametric(x)
            gf = genfgrad_parametric(x)
            hI, hJ, hf = gen_hessian_sparse_color_parametric(x)
            push!(l.valfuncs, f)
            push!(l.gradfuncs, gf)
            push!(l.hessfuncs, hf)
            push!(l.hessIJ, (hI, hJ))
            appendToIJ!(I,J,hI,hJ,x)
        else
            # check if there's a 1-1 mapping from reference indices to
            # indices in this expression
            refidx = l.referenceExpr[x.hashval] 
            ref = l.exprs[refidx]
            matches = (length(x.indexlist) == length(ref.indexlist))
            if matches
                for k in 1:length(x.indexlist)
                    if x.maptocanonical[x.indexlist[k]] != ref.maptocanonical[ref.indexlist[k]]
                        matches = false
                        break
                    end
                end
            end
            if matches
                # re-use AD from previous expression
                push!(l.valfuncs, l.valfuncs[refidx])
                push!(l.gradfuncs, l.gradfuncs[refidx])
                push!(l.hessfuncs, l.hessfuncs[refidx])
                hI,hJ = l.hessIJ[refidx]
                push!(l.hessIJ, (hI,hJ))
                appendToIJ!(I,J,hI,hJ,x)
            else
                # compute from scratch
                # TODO: actually we can share the AD but not the coloring here
                f = genfval_parametric(x)
                gf = genfgrad_parametric(x)
                hI, hJ, hf = gen_hessian_sparse_color_parametric(x)
                push!(l.valfuncs, f)
                push!(l.gradfuncs, gf)
                push!(l.hessfuncs, hf)
                push!(l.hessIJ, (hI, hJ))
                appendToIJ!(I,J,hI,hJ,x)
            end
        end
    end

    return I,J
end

# compute the value of each expression at xval, writing the result into out
function eval_g!(out::AbstractVector, l::ExprList, xval)
    @assert length(l.valfuncs) == length(l.exprs)

    for i in 1:length(l.exprs)
        out[i] = l.valfuncs[i](xval, IdentityArray(), l.exprs[i].inputvals)
    end
end

# evaluate the jacobian, putting the *sparse* result into V
function eval_jac_g!(V::AbstractVector, l::ExprList, xval)
    @assert length(l.gradfuncs) == length(l.exprs)

    idx = 1
    for (i,x) in enumerate(l.exprs)
        k = length(x.maptocanonical)
        l.gradfuncs[i](xval, IdentityArray(), sub(V, idx:(idx+k-1)), x.maptocanonical, x.inputvals)
        idx += k
    end

end

function jac_nz(l::ExprList)

    nnz_jac = 0
    for x in l.exprs
        nnz_jac += length(x.mapfromcanonical)
    end
    I = Int[]; J = Int[]
    sizehint(I,nnz_jac)
    sizehint(J,nnz_jac)
    for (i,x) in enumerate(l.exprs)
        for k in x.mapfromcanonical
            push!(I,i)
            push!(J,k)
        end
    end

    return I,J
end

# evaluate the sum of the hessian of each expression at xval
# each term weighted by lambda
function eval_hess!(V::AbstractVector, l::ExprList, xval, lambda)
    
    idx = 1
    for (i,x) in enumerate(l.exprs)
        nnz = length(l.hessIJ[i][1])
        subarr = sub(V,idx:(idx+nnz-1))
        l.hessfuncs[i](xval, subarr, x)
        scale!(subarr, lambda[i])
        idx += nnz
    end

end

