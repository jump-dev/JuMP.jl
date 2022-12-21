#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

"""
    SymmetricMatrixSpace()

Use in the [`@variable`](@ref) macro to constrain a matrix of variables to be
symmetric.

## Examples

```jldoctest; setup=:(model = Model())
julia> @variable(model, Q[1:2, 1:2] in SymmetricMatrixSpace())
2×2 LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}:
 Q[1,1]  Q[1,2]
 Q[1,2]  Q[2,2]
```
"""
struct SymmetricMatrixSpace end

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

julia> cref = @constraint(model, a >= b, PSDCone())
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
We see in the output of the last command that the vectorization of the matrix
is constrained to belong to the `PositiveSemidefiniteConeSquare`.

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

function build_constraint(
    _error::Function,
    f::AbstractMatrix{<:AbstractJuMPScalar},
    s::MOI.GreaterThan,
    extra::PSDCone,
)
    @assert iszero(s.lower)
    return build_constraint(_error, f, extra)
end

function build_constraint(
    _error::Function,
    f::AbstractMatrix{<:AbstractJuMPScalar},
    s::MOI.LessThan,
    extra::PSDCone,
)
    @assert iszero(s.upper)
    new_f = _MA.operate!!(*, -1, f)
    return build_constraint(_error, new_f, extra)
end

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

# This is a special method because calling `Matrix(matrix)` accesses an undef
# reference.
function vectorize(matrix::UpperTriangular, ::SquareMatrixShape)
    n = LinearAlgebra.checksquare(matrix)
    return [matrix[i, j] for j in 1:n for i in 1:n]
end

# This is a special method because calling `Matrix(matrix)` accesses an undef
# reference.
function vectorize(matrix::LowerTriangular, ::SquareMatrixShape)
    n = LinearAlgebra.checksquare(matrix)
    return [matrix[i, j] for j in 1:n for i in 1:n]
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
    build_variable(_error::Function, variables, ::SymmetricMatrixSpace)

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
    variables::Matrix{<:AbstractVariable},
    ::SymmetricMatrixSpace,
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
    variables::Matrix{<:AbstractVariable},
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
    variables::Matrix{<:AbstractVariable},
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

function value(
    Q::LinearAlgebra.Symmetric{V,Matrix{V}},
) where {V<:AbstractVariableRef}
    return LinearAlgebra.Symmetric(
        value.(LinearAlgebra.parent(Q)),
        LinearAlgebra.sym_uplo(Q.uplo),
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
    build_constraint(
        _error::Function,
        Q::AbstractMatrix{<:AbstractJuMPScalar},
        ::PSDCone,
    )

Return a `VectorConstraint` of shape [`SquareMatrixShape`](@ref) constraining
the matrix `Q` to be symmetric and positive semidefinite.

This function is used by the [`@constraint`](@ref) macro as follows:
```julia
@constraint(model, Q in PSDCone())
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

"""
    HermitianPSDCone

Hermitian positive semidefinite cone object that can be used to create a
Hermitian positive semidefinite square matrix in the [`@variable`](@ref)
and [`@constraint`](@ref) macros.

## Examples

Consider the following example:
```jldoctest; setup = :(using JuMP)
julia> model = Model();

julia> @variable(model, H[1:3, 1:3] in HermitianPSDCone())
3×3 Matrix{GenericAffExpr{ComplexF64, VariableRef}}:
 real(H[1,1])                                real(H[1,2]) + (0.0 + 1.0im) imag(H[1,2])   real(H[1,3]) + (0.0 + 1.0im) imag(H[1,3])
 real(H[1,2]) + (-0.0 - 1.0im) imag(H[1,2])  real(H[2,2])                                real(H[2,3]) + (0.0 + 1.0im) imag(H[2,3])
 real(H[1,3]) + (-0.0 - 1.0im) imag(H[1,3])  real(H[2,3]) + (-0.0 - 1.0im) imag(H[2,3])  real(H[3,3])

 julia> v = all_variables(model)
 9-element Vector{VariableRef}:
  real(H[1,1])
  real(H[1,2])
  real(H[2,2])
  real(H[1,3])
  real(H[2,3])
  real(H[3,3])
  imag(H[1,2])
  imag(H[1,3])
  imag(H[2,3])

julia> all_constraints(model, Vector{VariableRef}, MOI.HermitianPositiveSemidefiniteConeTriangle)
1-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.VectorOfVariables, MathOptInterface.HermitianPositiveSemidefiniteConeTriangle}}}:
 [real(H[1,1]), real(H[1,2]), real(H[2,2]), real(H[1,3]), real(H[2,3]), real(H[3,3]), imag(H[1,2]), imag(H[1,3]), imag(H[2,3])] in MathOptInterface.HermitianPositiveSemidefiniteConeTriangle(3)
```
We see in the output of the last commands that 9 real variables were created.
The matrix `H` contrains affine expressions in terms of these 9 variables that
parametrize a Hermitian matrix.
"""
struct HermitianPSDCone end

"""
    HermitianMatrixShape

Shape object for a Hermitian square matrix of `side_dimension` rows and
columns. The vectorized form corresponds to
[`MOI.HermitianPositiveSemidefiniteConeTriangle`](@ref).
"""
struct HermitianMatrixShape <: AbstractShape
    side_dimension::Int
end

function vectorize(matrix, shape::HermitianMatrixShape)
    return vectorize(Matrix(matrix), shape)
end

function vectorize(matrix::Matrix, ::HermitianMatrixShape)
    n = LinearAlgebra.checksquare(matrix)
    return vcat(
        vectorize(_real.(matrix), SymmetricMatrixShape(n)),
        vectorize(
            _imag.(matrix[1:(end-1), 2:end]),
            SymmetricMatrixShape(n - 1),
        ),
    )
end

function reshape_set(
    ::MOI.HermitianPositiveSemidefiniteConeTriangle,
    ::HermitianMatrixShape,
)
    return HermitianPSDCone()
end

function reshape_vector(v::Vector{T}, shape::HermitianMatrixShape) where {T}
    NewType = _MA.promote_operation(_MA.add_mul, T, Complex{Bool}, T)
    n = shape.side_dimension
    matrix = Matrix{NewType}(undef, n, n)
    real_k = 0
    imag_k = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    for j in 1:n
        for i in 1:(j-1)
            real_k += 1
            imag_k += 1
            matrix[i, j] = v[real_k] + im * v[imag_k]
            matrix[j, i] = v[real_k] - im * v[imag_k]
        end
        real_k += 1
        matrix[j, j] = v[real_k]
    end
    return matrix
end

function _vectorize_complex_variables(_error::Function, matrix::Matrix)
    n = LinearAlgebra.checksquare(matrix)
    for j in 1:n
        if !_isreal(matrix[j, j])
            _error(
                "Non-real bounds or starting values for diagonal of Hermitian variable.",
            )
        end
        for i in 1:j
            if matrix[i, j] != _conj(matrix[j, i])
                _error(
                    "Non-conjugate bounds or starting values for Hermitian variable.",
                )
            end
        end
    end
    return vectorize(matrix, HermitianMatrixShape(n))
end

function build_variable(
    _error::Function,
    variables::Matrix{<:AbstractVariable},
    ::HermitianPSDCone,
)
    n = _square_side(_error, variables)
    set = MOI.HermitianPositiveSemidefiniteConeTriangle(n)
    shape = HermitianMatrixShape(n)
    if any(_is_binary, variables) || any(_is_integer, variables)
        # We would then need to fix the imaginary value to zero. Let's wait to
        # see if there is need for such complication first.
        _error(
            "Binary or integer variables in a Hermitian matrix is not supported.",
        )
    end
    return VariablesConstrainedOnCreation(
        _vectorize_complex_variables(_error, variables),
        set,
        shape,
    )
end

"""
    build_constraint(
        _error::Function,
        Q::Hermitian{V,M},
        ::HermitianPSDCone,
    ) where {V<:AbstractJuMPScalar,M<:AbstractMatrix{V}}

Return a `VectorConstraint` of shape [`HermitianMatrixShape`](@ref) constraining
the matrix `Q` to be Hermitian positive semidefinite.

This function is used by the [`@constraint`](@ref) macros as follows:
```julia
@constraint(model, Hermitian(Q) in HermitianPSDCone())
```
"""
function build_constraint(
    ::Function,
    Q::Hermitian{V,M},
    ::HermitianPSDCone,
) where {V<:AbstractJuMPScalar,M<:AbstractMatrix{V}}
    n = LinearAlgebra.checksquare(Q)
    shape = HermitianMatrixShape(n)
    return VectorConstraint(
        vectorize(Q, shape),
        MOI.HermitianPositiveSemidefiniteConeTriangle(n),
        shape,
    )
end
