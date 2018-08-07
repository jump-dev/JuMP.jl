# Used in @constraint model x in PSDCone
struct PSDCone end

"""
    SymmetricMatrixShape

Shape object for a symmetric square matrix of `side_dimension` rows and columns.
The vectorized form contains the entries of the upper-right triangular part of
the matrix given column by column (or equivalently, the entries of the
lower-left triangular part given row by row).
"""
struct SymmetricMatrixShape <: AbstractShape
    side_dimension::Int
end
function reshape(vectorized_form::Vector{T}, shape::SymmetricMatrixShape) where T
    matrix = Matrix{T}(undef, shape.side_dimension, shape.side_dimension)
    k = 0
    for j in 1:shape.side_dimension
        for i in 1:j
            k += 1
            matrix[j, i] = matrix[i, j] = vectorized_form[k]
        end
    end
    return Symmetric(matrix)
end

"""
    SquareMatrixShape

Shape object for a square matrix of `side_dimension` rows and columns. The
vectorized form contains the entries of the the matrix given column by column
(or equivalently, the entries of the lower-left triangular part given row by
row).
"""
struct SquareMatrixShape <: AbstractShape
    side_dimension::Int
end
function reshape(vectorized_form::Vector{T}, shape::SquareMatrixShape) where T
    return Base.reshape(vectorized_form,
                        shape.side_dimension,
                        shape.side_dimension)
end

# Used by the @variable macro. It can also be used with the @constraint macro,
# this allows to get the constraint reference, e.g.
# @variable model x[1:2,1:2] Symmetric # x is Symmetric{VariableRef,Matrix{VariableRef}}
# varpsd = @constraint model x in PSDCone()
function buildconstraint(_error::Function, Q::Symmetric{V, Matrix{V}}, ::PSDCone) where V<:AbstractVariableRef
    n = Base.LinAlg.checksquare(Q)
    VectorOfVariablesConstraint([Q[i, j] for j in 1:n for i in 1:j],
                                MOI.PositiveSemidefiniteConeTriangle(n),
                                SymmetricMatrixShape(n))
end
# @variable model x[1:2,1:2] # x is Matrix{VariableRef}
# varpsd = @constraint model x in PSDCone()
function buildconstraint(_error::Function, Q::Matrix{<:AbstractVariableRef}, ::PSDCone)
    n = Base.LinAlg.checksquare(Q)
    VectorOfVariablesConstraint(vec(Q),
                                MOI.PositiveSemidefiniteConeSquare(n),
                                SquareMatrixShape(n))
end

function buildconstraint(_error::Function, x::AbstractMatrix, ::PSDCone)
    n = Base.LinAlg.checksquare(x)
    # Support for non-symmetric matrices as done prior to JuMP v0.19
    # will be added once the appropriate cone has been added in MathOptInterface
    # as discussed in the following PR:
    # https://github.com/JuliaOpt/JuMP.jl/pull/1122#issuecomment-344980944
    @assert issymmetric(x)
    aff = [x[i, j] for j in 1:n for i in 1:j]
    return VectorAffExprConstraint(aff,
                                   MOI.PositiveSemidefiniteConeTriangle(n),
                                   SymmetricMatrixShape(n))
end
