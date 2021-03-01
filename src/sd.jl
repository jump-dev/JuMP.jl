"""
    SymMatrixSpace()

Use in the [`@variable`](@ref) macro to constrain a matrix of variables to be
symmetric.

## Examples

```jldoctest; setup=:(model = Model())
julia> @variable(model, Q[1:2, 1:2] in SymMatrixSpace())
2×2 LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}:
 Q[1,1]  Q[1,2]
 Q[1,2]  Q[2,2]
```
"""
struct SymMatrixSpace end

"""
    SkewSymmetricMatrixSpace()

Use in the [`@variable`](@ref) macro to constrain a matrix of variables to be
skew-symmetric.

## Examples

```jldoctest; setup=:(model = Model())
@variable(model, Q[1:2, 1:2] in SkewSymmetricMatrixSpace())
```
"""
struct SkewSymmetricMatrixSpace end

"""
    PSDCone

Positive semidefinite cone object that can be used to constrain a square matrix
to be positive semidefinite in the [`@constraint`](@ref) macro. If the matrix
has type `Symmetric` then the columns vectorization (the vector obtained by
concatenating the columns) of its upper triangular part is constrained to belong
to the `MOI.PositiveSemidefiniteConeTriangle` set, otherwise its column
vectorization is constrained to belong to the
`MOI.PositiveSemidefiniteConeSquare` set.

## Examples

Consider the following example:
```jldoctest PSDCone; setup = :(using JuMP)
julia> model = Model();

julia> @variable(model, x)
x

julia> a = [ x 2x
            2x  x];

julia> b = [1 2
            2 4];

julia> cref = @SDconstraint(model, a ⪰ b)
[x - 1    2 x - 2;
 2 x - 2  x - 4  ] ∈ PSDCone()

julia> jump_function(constraint_object(cref))
4-element Array{GenericAffExpr{Float64,VariableRef},1}:
 x - 1
 2 x - 2
 2 x - 2
 x - 4

julia> moi_set(constraint_object(cref))
MathOptInterface.PositiveSemidefiniteConeSquare(2)
```
We see in the output of the last command that the matrix the vectorization of the
matrix is constrained to belong to the `PositiveSemidefiniteConeSquare`.

```jldoctest PSDCone
julia> using LinearAlgebra # For Symmetric

julia> cref = @constraint(model, Symmetric(a - b) in PSDCone())
[x - 1    2 x - 2;
 2 x - 2  x - 4  ] ∈ PSDCone()

julia> jump_function(constraint_object(cref))
3-element Array{GenericAffExpr{Float64,VariableRef},1}:
 x - 1
 2 x - 2
 x - 4

julia> moi_set(constraint_object(cref))
MathOptInterface.PositiveSemidefiniteConeTriangle(2)
```
As we see in the output of the last command, the vectorization of only the upper
triangular part of the matrix is constrained to belong to the
`PositiveSemidefiniteConeSquare`.
"""
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
function reshape_vector(
    vectorized_form::Vector{T},
    shape::SymmetricMatrixShape,
) where {T}
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
function reshape_set(
    ::MOI.PositiveSemidefiniteConeTriangle,
    ::SymmetricMatrixShape,
)
    return PSDCone()
end
function vectorize(matrix::Matrix, ::SymmetricMatrixShape)
    n = LinearAlgebra.checksquare(matrix)
    return [matrix[i, j] for j in 1:n for i in 1:j]
end

"""
    SkewSymmetricMatrixShape

Shape object for a skew symmetric square matrix of `side_dimension` rows and
columns. The vectorized form contains the entries of the upper-right triangular
part of the matrix (without the diagonal) given column by column (or
equivalently, the entries of the lower-left triangular part given row by row).
The diagonal is zero.
"""
struct SkewSymmetricMatrixShape <: AbstractShape
    side_dimension::Int
end
function reshape_vector(
    vectorized_form::Vector{T},
    shape::SkewSymmetricMatrixShape,
) where {T}
    NewType = Base.promote_type(
        T,
        _MA.promote_operation(-, T),
        _MA.promote_operation(zero, T),
    )
    matrix = Matrix{NewType}(undef, shape.side_dimension, shape.side_dimension)
    k = 0
    for j in 1:shape.side_dimension
        for i in 1:(j-1)
            k += 1
            matrix[i, j] = vectorized_form[k]
            matrix[j, i] = -vectorized_form[k]
        end
        matrix[j, j] = zero(NewType)
    end
    return matrix
end

function vectorize(matrix::Matrix, ::SkewSymmetricMatrixShape)
    n = LinearAlgebra.checksquare(matrix)
    return [matrix[i, j] for j in 1:n for i in 1:j-1]
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
function reshape_vector(
    vectorized_form::Vector{T},
    shape::SquareMatrixShape,
) where {T}
    return reshape(vectorized_form, shape.side_dimension, shape.side_dimension)
end
function reshape_set(::MOI.PositiveSemidefiniteConeSquare, ::SquareMatrixShape)
    return PSDCone()
end
vectorize(matrix::Matrix, ::SquareMatrixShape) = vec(matrix)

function vectorize(matrix, shape::Union{SymmetricMatrixShape,SquareMatrixShape})
    return vectorize(Matrix(matrix), shape)
end

function _square_side(_error::Function, ::Containers.SparseAxisArray)
    return _error("Cannot have index dependencies in symmetric variables.")
end
function _square_side(_error::Function, ::Containers.DenseAxisArray)
    return _error(
        "Index sets for symmetric variables must be ranges of the form 1:N.",
    )
end
function _square_side(_error::Function, ::Array)
    return _error("Symmetric variables must be 2-dimensional.")
end
function _square_side(_error::Function, variables::Matrix)
    n, m = size(variables)
    if n != m
        _error("Symmetric variables must be square. Got size ($n, $m).")
    end
    return n
end

function _vectorize_variables(_error::Function, matrix::Matrix)
    n = LinearAlgebra.checksquare(matrix)
    for j in 1:n
        for i in 1:j
            if matrix[i, j] != matrix[j, i]
                _error(
                    "Non-symmetric bounds, integrality or starting values for symmetric variable.",
                )
            end
        end
    end
    return vectorize(matrix, SymmetricMatrixShape(n))
end

"""
    build_variable(_error::Function, variables, ::SymMatrixSpace)

Return a `VariablesConstrainedOnCreation` of shape [`SymmetricMatrixShape`](@ref)
creating variables in `MOI.Reals`, i.e. "free" variables unless they are
constrained after their creation.

This function is used by the [`@variable`](@ref) macro as follows:
```julia
@variable(model, Q[1:2, 1:2], Symmetric)
```
"""
function build_variable(
    _error::Function,
    variables::Matrix{<:ScalarVariable},
    ::SymMatrixSpace,
)
    n = _square_side(_error, variables)
    set = MOI.Reals(MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n)))
    shape = SymmetricMatrixShape(n)
    return VariablesConstrainedOnCreation(
        _vectorize_variables(_error, variables),
        set,
        shape,
    )
end

"""
    build_variable(_error::Function, variables, ::SkewSymmetricMatrixSpace)

Return a `VariablesConstrainedOnCreation` of shape [`SkewSymmetricMatrixShape`](@ref)
creating variables in `MOI.Reals`, i.e. "free" variables unless they are
constrained after their creation.

This function is used by the [`@variable`](@ref) macro as follows:
```julia
@variable(model, Q[1:2, 1:2] in SkewSymmetricMatrixSpace())
```
"""
function build_variable(
    _error::Function,
    variables::Matrix{<:ScalarVariable},
    ::SkewSymmetricMatrixSpace,
)
    n = _square_side(_error, variables)
    set = MOI.Reals(div(n^2 - n, 2))
    shape = SkewSymmetricMatrixShape(n)
    return VariablesConstrainedOnCreation(
        vectorize(variables, SkewSymmetricMatrixShape(n)),
        set,
        shape,
    )
end

"""
    build_variable(_error::Function, variables, ::PSDCone)

Return a `VariablesConstrainedOnCreation` of shape [`SymmetricMatrixShape`](@ref)
constraining the variables to be positive semidefinite.

This function is used by the [`@variable`](@ref) macro as follows:
```julia
@variable(model, Q[1:2, 1:2], PSD)
```
"""
function build_variable(
    _error::Function,
    variables::Matrix{<:ScalarVariable},
    ::PSDCone,
)
    n = _square_side(_error, variables)
    set = MOI.PositiveSemidefiniteConeTriangle(n)
    shape = SymmetricMatrixShape(n)
    return VariablesConstrainedOnCreation(
        _vectorize_variables(_error, variables),
        set,
        shape,
    )
end

"""
    build_constraint(_error::Function, Q::Symmetric{V, M},
                     ::PSDCone) where {V <: AbstractJuMPScalar,
                                       M <: AbstractMatrix{V}}

Return a `VectorConstraint` of shape [`SymmetricMatrixShape`](@ref) constraining
the matrix `Q` to be positive semidefinite.

This function is used by the [`@constraint`](@ref) macros as follows:
```julia
@constraint(model, Symmetric(Q) in PSDCone())
```
The form above is usually used when the entries of `Q` are affine or quadratic
expressions, but it can also be used when the entries are variables to get the
reference of the semidefinite constraint, e.g.,
```julia
@variable model Q[1:2,1:2] Symmetric
# The type of `Q` is `Symmetric{VariableRef, Matrix{VariableRef}}`
var_psd = @constraint model Q in PSDCone()
# The `var_psd` variable contains a reference to the constraint
```
"""
function build_constraint(
    _error::Function,
    Q::Symmetric{V,M},
    ::PSDCone,
) where {V<:AbstractJuMPScalar,M<:AbstractMatrix{V}}
    n = LinearAlgebra.checksquare(Q)
    shape = SymmetricMatrixShape(n)
    return VectorConstraint(
        vectorize(Q, shape),
        MOI.PositiveSemidefiniteConeTriangle(n),
        shape,
    )
end

"""
    build_constraint(_error::Function,
                     Q::AbstractMatrix{<:AbstractJuMPScalar},
                     ::PSDCone)

Return a `VectorConstraint` of shape [`SquareMatrixShape`](@ref) constraining
the matrix `Q` to be symmetric and positive semidefinite.

This function is used by the [`@constraint`](@ref) and [`@SDconstraint`](@ref)
macros as follows:
```julia
@constraint(model, Q in PSDCone())
@SDconstraint(model, P ⪰ Q)
```
The [`@constraint`](@ref) call above is usually used when the entries of `Q` are
affine or quadratic expressions, but it can also be used when the entries are
variables to get the reference of the semidefinite constraint, e.g.,
```julia
@variable model Q[1:2,1:2]
# The type of `Q` is `Matrix{VariableRef}`
var_psd = @constraint model Q in PSDCone()
# The `var_psd` variable contains a reference to the constraint
```
"""
function build_constraint(
    _error::Function,
    Q::AbstractMatrix{<:AbstractJuMPScalar},
    ::PSDCone,
)
    n = LinearAlgebra.checksquare(Q)
    shape = SquareMatrixShape(n)
    return VectorConstraint(
        vectorize(Q, shape),
        MOI.PositiveSemidefiniteConeSquare(n),
        shape,
    )
end
