# Helps manage a set of expressions
# Tries to share AD across sets of expressions with identical structure
# Useful for computing the Hessian of the lagrangian and related terms


export ExprList, prep_sparse_hessians, eval_g!,
    eval_jac_g!, jac_nz, eval_hess!, eval_hessvec!,
    prep_expression_output

type ExprList
    exprs::Vector{SymbolicOutput}
    idxfuncs::Vector{Function}
    valfuncs::Vector{Function}
    gradfuncs::Vector{Function}
    hess_matmat_funcs::Vector{Function}
    hess_IJ_funcs::Vector{Function}
    hessfuncs::Vector{Function} # these can't be shared if coloring is different
    hessIJ::Vector{Tuple{Vector{Int},Vector{Int}}}
    exprfuncs::Vector{Function}
    dualvec::Vector{Dual4{Float64}}
    dualout::Vector{Dual4{Float64}}
end


ExprList() = ExprList(SymbolicOutput[], Function[], Function[], Function[], Function[], Function[], Function[], Array(Tuple{Vector{Int},Vector{Int}},0), Function[], Array(Dual4{Float64},0), Array(Dual4{Float64},0))

Base.push!(l::ExprList, s::SymbolicOutput) = push!(l.exprs, s)
Base.getindex(l::ExprList, i) = l.exprs[i]

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
# if hessvec_only is true, we just prepare the gradients and jacobians
function prep_sparse_hessians(l::ExprList, num_total_vars; hessvec_only::Bool=false)
    # clear current state
    empty!(l.idxfuncs)
    empty!(l.valfuncs)
    empty!(l.gradfuncs)
    empty!(l.hess_matmat_funcs)
    empty!(l.hess_IJ_funcs)
    empty!(l.hessfuncs)
    empty!(l.hessIJ)

    N = length(l.exprs)
    sizehint!(l.valfuncs, N)
    sizehint!(l.gradfuncs, N)
    sizehint!(l.hessfuncs, N)
    referenceExpr = Dict()

    genindexlist_t = 0.0
    indexlist_t = 0.0
    genfgrad_t = 0.0
    coloring_t = 0.0

    # shared temporary storage for hessian-vector products
    l.dualvec = Array(Dual4{Float64}, ceil(Int,num_total_vars/2))
    l.dualout = Array(Dual4{Float64}, ceil(Int,num_total_vars/2))
    #             ""           for sparsity pattern detection
    idxset = IndexedSet(num_total_vars)


    I = Int[]
    J = Int[]

    for (i,x) in enumerate(l.exprs)
        if !haskey(referenceExpr, x.hashval)
            referenceExpr[x.hashval] = i
            tic()
            idxf = genindexlist_parametric(x)
            genindexlist_t += toq()
            tic()
            prepare_indexmap(x, idxf(idxset,x.inputvals...))
            x.indexlist = idxf(Int[],x.inputvals...)
            indexlist_t += toq()
            f = genfval_parametric(x)
            tic()
            gf = genfgrad_parametric(x)
            genfgrad_t += toq()

            push!(l.idxfuncs, idxf)
            push!(l.valfuncs, f)
            push!(l.gradfuncs, gf)
            if !hessvec_only
                hess_matmat = gen_hessian_matmat_parametric(x, gf)
                hess_IJf = compute_hessian_sparsity_IJ_parametric(x)
                tic()
                hI, hJ, hf = gen_hessian_sparse_color_parametric(x, num_total_vars, hess_matmat, hess_IJf, l.dualvec, l.dualout, idxset)
                coloring_t += toq()
                push!(l.hess_matmat_funcs, hess_matmat)
                push!(l.hess_IJ_funcs, hess_IJf)
                push!(l.hessfuncs, hf)
                push!(l.hessIJ, (hI, hJ))
                appendToIJ!(I,J,hI,hJ,x)
            end
        else
            # check if there's a 1-1 mapping from reference indices to
            # indices in this expression
            refidx = referenceExpr[x.hashval]
            ref = l.exprs[refidx]
            tic()
            x.indexlist = l.idxfuncs[refidx](Int[],x.inputvals...)
            process_indexmap(x,x.indexlist)
            indexlist_t += toq()
            # these are always shared
            push!(l.idxfuncs, l.idxfuncs[refidx])
            push!(l.valfuncs, l.valfuncs[refidx])
            push!(l.gradfuncs, l.gradfuncs[refidx])
            hessvec_only && continue

            push!(l.hess_matmat_funcs, l.hess_matmat_funcs[refidx])
            push!(l.hess_IJ_funcs, l.hess_IJ_funcs[refidx])

            matches = (length(x.indexlist) == length(ref.indexlist))
            if matches
                # Could we just hash the indexlist sequence instead of storing it?
                idxlist_x = x.indexlist::Vector{Int}
                idxlist_ref = ref.indexlist::Vector{Int}
                for k in 1:length(idxlist_x)
                    if x.maptocanonical[idxlist_x[k]] != ref.maptocanonical[idxlist_ref[k]]
                        matches = false
                        break
                    end
                end
            end
            if matches
                # re-use AD and coloring from previous expression
                push!(l.hessfuncs, l.hessfuncs[refidx])
                hI,hJ = l.hessIJ[refidx]
                push!(l.hessIJ, (hI,hJ))
                appendToIJ!(I,J,hI,hJ,x)
            else
                # we can share the AD but not the coloring here
                tic()
                hI, hJ, hf = gen_hessian_sparse_color_parametric(x, num_total_vars, l.hess_matmat_funcs[refidx], l.hess_IJ_funcs[refidx], l.dualvec, l.dualout, idxset)
                coloring_t += toq()
                push!(l.hessfuncs, hf)
                push!(l.hessIJ, (hI, hJ))
                appendToIJ!(I,J,hI,hJ,x)
            end
        end
    end

    #=
    @show genindexlist_t
    @show indexlist_t
    @show genfgrad_t
    @show coloring_t
    =#

    return I,J
end

function prep_expression_output(l::ExprList)
    # clear current state
    empty!(l.exprfuncs)

    N = length(l.exprs)
    sizehint!(l.exprfuncs, N)
    referenceExpr = Dict()

    for (i,x) in enumerate(l.exprs)
        if !haskey(referenceExpr, x.hashval)
            referenceExpr[x.hashval] = i
            push!(l.exprfuncs, genfexpr_parametric(x))
        else
            refidx = referenceExpr[x.hashval]
            push!(l.exprfuncs, l.exprfuncs[refidx])
        end
    end

    return
end

# compute the value of each expression at xval, writing the result into out
function eval_g!(out::AbstractVector, l::ExprList, xval)
    @assert length(l.valfuncs) == length(l.exprs)

    for i in 1:length(l.exprs)
        out[i] = l.valfuncs[i](xval, IdentityArray(), l.exprs[i].inputvals...)::eltype(out)
    end
end

# evaluate the jacobian, putting the *sparse* result into V
function eval_jac_g!(V::DenseVector, l::ExprList, xval)
    @assert length(l.gradfuncs) == length(l.exprs)

    idx = 1
    p = pointer(V)
    for i in 1:length(l.exprs)
        x = l.exprs[i]
        k = length(x.maptocanonical)::Int
        #Vsub = subarr(V, idx:(idx+k-1))
        Vsub = VectorView(idx-1, k, p)
        l.gradfuncs[i](xval, IdentityArray(), Vsub, x.maptocanonical, x.inputvals...)
        idx += k
    end

end

# weighted hessian-vector product h += ∑_i μ_i*∇g_i(xval)*v
function eval_hessvec!{T}(h,v::DenseVector{T}, l::ExprList, xval, μ)
    idx = 1
    dualvec = reinterpret(Dual{T}, l.dualvec)
    dualout = reinterpret(Dual{T}, l.dualout)
    @assert length(dualvec) >= length(xval)
    @assert length(dualout) >= length(xval)
    for k in 1:length(l.exprs)
        x = l.exprs[k]
        for i in 1:length(xval)
            dualvec[i] = Dual(xval[i],v[i])
            dualout[i] = zero(Dual{T})
        end
        l.gradfuncs[k](dualvec, IdentityArray(), dualout, IdentityArray(), x.inputvals...)
        for i in 1:length(xval)
            h[i] += μ[k]*epsilon(dualout[i])
        end
    end
end

function jac_nz(l::ExprList)

    nnz_jac = 0
    for x in l.exprs
        nnz_jac += length(x.mapfromcanonical)
    end
    I = Int[]; J = Int[]
    sizehint!(I,nnz_jac)
    sizehint!(J,nnz_jac)
    for (i,x) in enumerate(l.exprs)
        for k in x.mapfromcanonical
            push!(I,i)
            push!(J,k)
        end
    end

    return I,J
end

# lightweight unsafe view for vectors
# it seems that the only way to avoid triggering
# allocations is to have only bitstype fields, so
# we store a pointer.
immutable VectorView{T}
    offset::Int
    len::Int
    ptr::Ptr{T}
end

Base.getindex(v::VectorView,idx) = unsafe_load(v.ptr, idx+v.offset)
Base.setindex!(v::VectorView,value,idx) = unsafe_store!(v.ptr,value,idx+v.offset)
Base.length(v::VectorView) = v.len

# evaluate the sum of the hessian of each expression at xval
# each term weighted by lambda
function eval_hess!(V::DenseVector, l::ExprList, xval, lambda)
    
    idx = 1
    p = pointer(V)
    for i in 1:length(l.exprs)
        nnz = length(l.hessIJ[i][1])
        #Vsub = subarr(V,idx:(idx+nnz-1))
        Vsub = VectorView(idx-1,nnz,p)
        l.hessfuncs[i](xval, Vsub, l.exprs[i])
        for k in idx:(idx+nnz-1)
            V[k] *= lambda[i]
        end
        #scale!(Vsub, lambda[i])
        idx += nnz
    end

end

# convert idx to a flat expression
function to_flat_expr(l::ExprList, idx)
    if length(l.exprfuncs) == 0
        error("Cannot generate expression, feature was not requested")
    end
    return l.exprfuncs[idx](l.exprs[idx].inputvals...)
end

