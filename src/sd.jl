# Used in @constraint model x in PSDCone
struct PSDCone end

# Used by the @variable macro. It can also be used with the @constraint macro,
# this allows to get the constraint reference, e.g.
# @variable model x[1:2,1:2] Symmetric # x is Symmetric{VariableRef,Matrix{VariableRef}}
# varpsd = @constraint model x in PSDCone()
function buildconstraint(_error::Function, Q::Symmetric{V, Matrix{V}}, ::PSDCone) where V<:AbstractVariableRef
    n = Base.LinAlg.checksquare(Q)
    VectorOfVariablesConstraint([Q[i, j] for j in 1:n for i in 1:j], MOI.PositiveSemidefiniteConeTriangle(n))
end
# @variable model x[1:2,1:2] # x is Matrix{VariableRef}
# varpsd = @constraint model x in PSDCone()
function buildconstraint(_error::Function, Q::Matrix{<:AbstractVariableRef}, ::PSDCone)
    n = Base.LinAlg.checksquare(Q)
    VectorOfVariablesConstraint(vec(Q), MOI.PositiveSemidefiniteConeSquare(n))
end

function buildconstraint(_error::Function, x::AbstractMatrix, ::PSDCone)
    n = Base.LinAlg.checksquare(x)
    # Support for non-symmetric matrices as done prior to JuMP v0.19
    # will be added once the appropriate cone has been added in MathOptInterface
    # as discussed in the following PR:
    # https://github.com/JuliaOpt/JuMP.jl/pull/1122#issuecomment-344980944
    @assert issymmetric(x)
    aff = [x[i, j] for j in 1:n for i in 1:j]
    s = MOI.PositiveSemidefiniteConeTriangle(n)
    return VectorAffExprConstraint(aff, s)
end
