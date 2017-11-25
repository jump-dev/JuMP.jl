# Used in @constraint m X in SDP
struct PSDCone end

##########################################################################
# SDConstraint is a (dual) semidefinite constraint of the form
# ∑ cᵢ Xᵢ ≥ D, where D is a n×n symmetric data matrix, cᵢ are
# scalars, and Xᵢ are n×n symmetric variable matrices. The inequality
# is taken w.r.t. the psd partial order.
mutable struct SDConstraint <: AbstractConstraint
    terms
end

# Special-case X ≥ 0, which is often convenient
function SDConstraint(lhs::AbstractMatrix, rhs::Number)
    rhs == 0 || error("Cannot construct a semidefinite constraint with nonzero scalar bound $rhs")
    SDConstraint(lhs)
end

"""
    addconstraint(m::Model, c::SDConstraint)

Add a SD constraint to `Model m`.
"""
function addconstraint(m::Model, c::SDConstraint)
    push!(m.sdpconstr,c)
    m.internalModelLoaded = false
    ConstraintRef{Model,SDConstraint}(m,length(m.sdpconstr))
end

# helper method for mapping going on below
Base.copy(x::Number, new_model::Model) = copy(x)

Base.copy(c::SDConstraint, new_model::Model) =
    SDConstraint(map(t -> copy(t, new_model), c.terms))

# Returns the number of rows used by SDP constraints in the MPB conic representation
# (excluding symmetry constraints)
#   Julia seems to not be able to infer the return type (probably because c.terms is Any)
#   so getNumSDPRows tries to call zero(Any)... Using ::Int solves this issue
function getNumRows(c::SDConstraint)::Int
    n = size(c.terms, 1)
    (n * (n+1)) ÷ 2
end
getNumSDPRows(m::Model) = sum(getNumRows.(m.sdpconstr))

# Returns the number of symmetry-enforcing constraints for SDP constraints
function getNumSymRows(m::Model)
    sum(map(length, m.sdpconstrSym))
end

# Let S₊ be the cone of symmetric semidefinite matrices in
# the n*(n+1)/2 dimensional space of symmetric R^{nxn} matrices.
# It is well known that S₊ is a self-dual proper cone.
# Let P₊ be the cone of symmetric semidefinite matrices in
# the n^2 dimensional space of R^{nxn} matrices and
# let D₊ be the cone of matrices A such that A+Aᵀ ∈ P₊.
# P₊ is not proper since it is not solid (as it is not n^2 dimensional) so it is not ensured that (P₊)** = P₊
# However this is the case since, as we will see, (P₊)* = D₊ and (D₊)* = P₊.
# * Let us first see why (P₊)* = D₊.
#   If B is symmetric, then ⟨A,B⟩ = ⟨Aᵀ,Bᵀ⟩ = ⟨Aᵀ,B⟩ so 2⟨A,B⟩ = ⟨A,B⟩ + ⟨Aᵀ,B⟩ = ⟨A+Aᵀ,B⟩
#   Therefore, ⟨A,B⟩ ⩾ 0 for all B ∈ P₊ if and only if ⟨A+Aᵀ,B⟩ ⩾ 0 for all B ∈ P₊
#   Since A+Aᵀ is symmetric and we know that S₊ is self-dual, we have shown that (P₊)*
#   is the set of matrices A such that A+Aᵀ is PSD
# * Let us now see why (D₊)* = P₊.
#   Since A ∈ D₊ implies that Aᵀ ∈ D₊, B ∈ (D₊)* means that ⟨A+Aᵀ,B⟩ ⩾ 0 for any A ∈ D₊ hence B is positive semi-definite.
#   To see why it should be symmetric, simply notice that if B[i,j] < B[j,i] then ⟨A,B⟩ can be made arbitrarily small by setting
#   A[i,j] += s
#   A[j,i] -= s
#   with s arbitrarilly large, and A stays in D₊ as A+Aᵀ does not change.
#
# Typically, SDP primal/dual are presented as
# min ⟨C, X⟩                                                                max ∑ b_ky_k
# ⟨A_k, X⟩ = b_k ∀k                                                         C - ∑ A_ky_k ∈ S₊
#        X ∈ S₊                                                                      y_k free ∀k
# Here, as we allow A_i to be non-symmetric, we should rather use
# min ⟨C, X⟩                                                                max ∑ b_ky_k
# ⟨A_k, X⟩ = b_k ∀k                                                         C - ∑ A_ky_k ∈ P₊
#        X ∈ D₊                                                                      y_k free ∀k
# which is implemented as
# min ⟨C, Z⟩ + (C[i,j]-C[j-i])s[i,j]                                        max ∑ b_ky_k
# ⟨A_k, Z⟩ + (A_k[i,j]-A_k[j,i])s[i,j] = b_k ∀k                   C+Cᵀ - ∑ (A_k+A_kᵀ)y_k ∈ S₊
#       s[i,j] free  1 ⩽ i,j ⩽ n with i > j     C[i,j]-C[j-i] - ∑ (A_k[i,j]-A_k[j,i])y_k = 0  1 ⩽ i,j ⩽ n with i > j
#        Z ∈ S₊                                                                      y_k free ∀k
# where "∈ S₊" only look at the diagonal and upper diagonal part.
# In the last primal program, we have the variables Z = X + Xᵀ and a upper triangular matrix S such that X = Z + S - Sᵀ

"""
    getdual(c::ConstraintRef{Model,SDConstraint})

"""
function getdual(c::ConstraintRef{Model,SDConstraint})
    dual, symdual = getconicdualaux(c.m, c.idx, true)
    n = size(c.m.sdpconstr[c.idx].terms, 1)
    X = Matrix{eltype(dual)}(n, n)
    @assert length(dual) == convert(Int, n*(n+1)/2)
    idx = 0
    for i in 1:n
        for j in i:n
            idx += 1
            if i == j
                X[i,j] = dual[idx]
            else
                X[j,i] = X[i,j] = dual[idx] / sqrt(2)
            end
        end
    end
    if !isempty(symdual)
        @assert length(symdual) == length(c.m.sdpconstrSym[c.idx])
        idx = 0
        for (i,j) in c.m.sdpconstrSym[c.idx]
            idx += 1
            s = symdual[idx]
            X[i,j] -= s
            X[j,i] += s
        end
    end
    X
end

# Returns a boolean vector indicating if variable in the model
# is an off-diagonal element of an SDP matrix.
# This is needed because we have to rescale coefficients that
# touch these variables.
function offdiagsdpvars(m::Model)
    offdiagvars = falses(m.numCols)
    for (name,idx) in m.varCones
        if name == :SDP
            conelen = length(idx)
            n = round(Int,sqrt(1/4+2*conelen)-1/2)
            @assert n*(n+1)/2 == conelen
            r = 1
            for i in 1:n
                for j in i:n
                    if i != j
                        offdiagvars[idx[r]] = true
                    end
                    r += 1
                end
            end
        end
    end
    return offdiagvars
end

function getSDrowsinfo(m::Model)
    # find starting column indices for sdp matrices
    nnz = 0
    numSDPRows = 0
    numSymRows = 0
    for c in m.sdpconstr
        n = size(c.terms,1)
        @assert n == size(c.terms,2)
        @assert ndims(c.terms) == 2
        numSDPRows += convert(Int, n*(n+1)/2)
        for i in 1:n, j in i:n
            nnz += length(c.terms[i,j].coeffs)
        end
        if !issymmetric(c.terms)
            # symmetry constraints
            numSymRows += convert(Int, n*(n-1)/2)
        end
    end
    numSDPRows, numSymRows, nnz
end

function rescaleSDcols!(f, J, V, m)
    # Objective coefficients and columns of A matrix are
    # rescaled for SDP variables
    offdiagvars = offdiagsdpvars(m)
    f[offdiagvars] /= sqrt(2)
    for k in 1:length(J)
        if offdiagvars[J[k]]
            V[k] /= sqrt(2)
        end
    end
end

function fillconstr!(I, J, V, b, con_cones, tmprow::IndexedVector, constr_to_row, c, d, constrs::Vector{SDConstraint}, m::Model, ignore_not_owned::Bool=false)
    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    sdpconstr_sym = Vector{Vector{Tuple{Int,Int}}}(length(constrs))
    sdpidx = 0
    for con in constrs
        sdpidx += 1
        sdp_start = c + 1
        n = size(con.terms,1)
        for i in 1:n, j in i:n
            c += 1
            terms::AffExpr = con.terms[i,j] + con.terms[j,i]
            collect_expr!(m, tmprow, terms, ignore_not_owned)
            nnz = tmprow.nnz
            indices = tmpnzidx[1:nnz]
            append!(I, fill(c, nnz))
            append!(J, indices)
            # scale to svec form
            scale = (i == j) ? 0.5 : 1/sqrt(2)
            append!(V, -scale*tmpelts[indices])
            b[c] = scale*terms.constant
        end
        push!(con_cones, (:SDP, sdp_start:c))
        constr_to_row[d + sdpidx] = collect(sdp_start:c)
        syms = Tuple{Int,Int}[]
        if !issymmetric(con.terms)
            sym_start = c + 1
            # add linear symmetry constraints
            for i in 1:n, j in 1:(i-1)
                collect_expr!(m, tmprow, con.terms[i,j] - con.terms[j,i], ignore_not_owned)
                nnz = tmprow.nnz
                # if the symmetry-enforcing row is empty or has only tiny coefficients due to unintended numerical asymmetry, drop it
                largestabs = 0.0
                for k in 1:nnz
                    largestabs = max(largestabs,abs(tmpelts[tmpnzidx[k]]))
                end
                if largestabs < 1e-10
                    continue
                end
                push!(syms, (i,j))
                c += 1
                indices = tmpnzidx[1:nnz]
                append!(I, fill(c, nnz))
                append!(J, indices)
                append!(V, tmpelts[indices])
                b[c] = 0
            end
            if c >= sym_start
                push!(con_cones, (:Zero, sym_start:c))
            end
            constr_to_row[d + length(constrs) + sdpidx] = collect(sym_start:c)
            @assert length(syms) == length(sym_start:c)
        else
            constr_to_row[d + length(constrs) + sdpidx] = Int[]
        end
        sdpconstr_sym[sdpidx] = syms
    end

    m.sdpconstrSym = sdpconstr_sym

    c, d + 2 * length(constrs)
end
