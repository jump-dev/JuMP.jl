export PSDCone
# Used in @constraint m X in PSDCone
struct PSDCone end

struct SDVariableConstraint <: AbstractConstraint
    Q::Matrix{JuMP.Variable}
end

# Used by the @variable macro. Currently cannot also be used through the @constraint macro because of the underscore
# It needs a larger discussion on whether we want to allow adding VectorOfVariable in cone using the @constraint macro.
function _constructconstraint!(Q::Matrix{JuMP.Variable}, ::PSDCone)
    #@assert issymmetric(Q) # TODO it could be nonsymmetric if used through the @constraint macro
    SDVariableConstraint(Q)
end

"""
    addconstraint(m::Model, c::SDVariableConstraint)

Add a SD variable constraint to `Model m`.
"""
function addconstraint(m::Model, c::SDVariableConstraint)
    @assert issymmetric(c.Q)
    @assert !m.solverinstanceattached # TODO
    n = Base.LinAlg.checksquare(c.Q)
    cref = MOI.addconstraint!(m.instance, MOI.VectorOfVariables([instancereference(c.Q[i, j]) for j in 1:n for i in 1:j]), MOI.PositiveSemidefiniteConeTriangle(n))
    return ConstraintRef(m, cref)
end

##########################################################################
# SDConstraint is a (dual) semidefinite constraint of the form
# ∑ cᵢ Xᵢ ⪰ D, where D is a n×n data matrix, cᵢ are scalars,
# and Xᵢ are n×n symmetric variable matrices. The inequality
# is taken w.r.t. the psd partial order.
struct SDConstraint{T<:AbstractJuMPScalar, MT<:AbstractMatrix{T}} <: AbstractConstraint
    terms::MT
end
SDConstraint(terms::MT) where {T<:AbstractJuMPScalar, MT<:AbstractMatrix{T}} = SDConstraint{T, MT}(terms)

"""
    trimap(i, j)

Returns the vectorized index `k` corresponding to the matrix indices `(i, j)`.
"""
function trimap(i::Integer, j::Integer)
    if i < j
        trimap(j, i)
    else
        div((i-1)*i, 2) + j
    end
end

"""
    invtrimap(i, j)

Returns the matrix indices `(i, j)` (`i <= j`) corresponding to the vectorized index `k`.
"""
function invtrimap(k::Integer)
    i = isqrt(2k)
    j = k - div((i-1)*i, 2)
    i, j
end

struct SDSymInstanceRef{T}
    sd::MOICON{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
    # symmetry enforcing constraints, it is null if the matrix is symmetric
    symref::Nullable{MOICON{MOI.VectorAffineFunction{T}, MOI.Zeros}}
    # if k in symidx then invtrimap(k) has a symmetry enforcing constraint in symref
    symidx::IntSet
end

constructconstraint!(x::AbstractMatrix, ::PSDCone) = SDConstraint(x)

_constrain_symmetry(m::Model, c::SDConstraint{T, MT}, sdref) where {T, MT<:Symmetric} = sdref

function _constrain_symmetry(m::Model, c::SDConstraint{T}, sdref) where T
    expr = Dict{Variable, Float64}()
    symaff = AffExpr[]
    symidx = IntSet()
    n = Base.LinAlg.checksquare(c.terms)
    for i in 1:n, j in 1:i
        empty!(expr)
        fill_expr!(m, expr, c.terms[i, j])
        fill_expr!(m, expr, c.terms[j, i], -1)
        constant = c.terms[i, j].constant - c.terms[j, i].constant
        # if the symmetry-enforcing row is empty or has only tiny coefficients due to unintended numerical asymmetry, drop it
        largestabs = abs(constant)
        for (var, coeff) in expr
            if iszero(coeff)
                delete!(expr, var)
            else
                largestabs = max(largestabs, abs(coeff))
            end
        end
        if largestabs < 1e-10
            continue
        end
        push!(symidx, trimap(i, j))
        aff = AffExpr(collect(keys(expr)), collect(values(expr)), constant)
        push!(symaff, aff)
    end
    if isempty(symaff)
        symref = Nullable{MOICON{MOI.VectorAffineFunction{Float64}, MOI.Zeros}}()
    else
        symcon = VectorAffExprConstraint(symaff, MOI.Zeros(length(symaff)))
        symref = Nullable{MOICON{MOI.VectorAffineFunction{Float64}, MOI.Zeros}}(addconstraint(m, symcon).instanceref)
    end
    SDSymInstanceRef(sdref, symref, symidx)
end

"""
    addconstraint(m::Model, c::SDConstraint)

Adds the semidefinite constraint `c` and returns a reference `cref` to the constraint.
If `c.terms` is a matrix of type `Symmetric`, then `cref.instanceref` will be an instance constraint reference
to the positive semidefiniteness of the matrix.
Otherwise, `cref.instanceref` is a `SDSymInstanceRef` containing the PSD constraint and the symmetry constraints.
"""
function addconstraint(m::Model, c::SDConstraint)
    @assert !m.solverinstanceattached # TODO
    n = Base.LinAlg.checksquare(c.terms)
    sdref = MOI.addconstraint!(m.instance, MOI.VectorAffineFunction([c.terms[i, j] for j in 1:n for i in 1:j]), MOI.PositiveSemidefiniteConeTriangle(n))
    cref = _constrain_symmetry(m, c, sdref)
    return ConstraintRef(m, cref)
end
