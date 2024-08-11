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
    vectorize(matrix::AbstractMatrix, ::Shape)

Convert the `matrix` into a vector according to `Shape`.
"""
function vectorize end

"""
    SymmetricMatrixSpace()

Use in the [`@variable`](@ref) macro to constrain a matrix of variables to be
symmetric.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, Q[1:2, 1:2] in SymmetricMatrixSpace())
2×2 LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}:
 Q[1,1]  Q[1,2]
 Q[1,2]  Q[2,2]
```
"""
struct SymmetricMatrixSpace end

"""
    SkewSymmetricMatrixSpace()

Use in the [`@variable`](@ref) macro to constrain a matrix of variables to be
skew-symmetric.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, Q[1:2, 1:2] in SkewSymmetricMatrixSpace())
2×2 Matrix{AffExpr}:
 0        Q[1,2]
 -Q[1,2]  0
```
"""
struct SkewSymmetricMatrixSpace end

"""
    HermitianMatrixSpace()

Use in the [`@variable`](@ref) macro to constrain a matrix of variables to be
hermitian.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, Q[1:2, 1:2] in HermitianMatrixSpace())
2×2 LinearAlgebra.Hermitian{GenericAffExpr{ComplexF64, VariableRef}, Matrix{GenericAffExpr{ComplexF64, VariableRef}}}:
 real(Q[1,1])                    real(Q[1,2]) + imag(Q[1,2]) im
 real(Q[1,2]) - imag(Q[1,2]) im  real(Q[2,2])
```
"""
struct HermitianMatrixSpace end

"""
    PSDCone

Positive semidefinite cone object that can be used to constrain a square matrix
to be positive semidefinite in the [`@constraint`](@ref) macro.

If the matrix has type `Symmetric` then the columns vectorization (the vector
obtained by concatenating the columns) of its upper triangular part is
constrained to belong to the [`MOI.PositiveSemidefiniteConeTriangle`](@ref) set,
otherwise its column vectorization is constrained to belong to the
[`MOI.PositiveSemidefiniteConeSquare`](@ref) set.

## Example

Non-symmetric case:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> a = [x 2x; 2x x];

julia> b = [1 2; 2 4];

julia> cref = @constraint(model, a >= b, PSDCone())
[x - 1    2 x - 2
 2 x - 2  x - 4] ∈ PSDCone()

julia> jump_function(constraint_object(cref))
4-element Vector{AffExpr}:
 x - 1
 2 x - 2
 2 x - 2
 x - 4

julia> moi_set(constraint_object(cref))
MathOptInterface.PositiveSemidefiniteConeSquare(2)
```

Symmetric case:

```jldoctest PSDCone
julia> using LinearAlgebra # For Symmetric

julia> model = Model();

julia> @variable(model, x);

julia> a = [x 2x; 2x x];

julia> b = [1 2; 2 4];

julia> cref = @constraint(model, Symmetric(a - b) in PSDCone())
[x - 1  2 x - 2
 ⋯      x - 4] ∈ PSDCone()

julia> jump_function(constraint_object(cref))
3-element Vector{AffExpr}:
 x - 1
 2 x - 2
 x - 4

julia> moi_set(constraint_object(cref))
MathOptInterface.PositiveSemidefiniteConeTriangle(2)
```
"""
struct PSDCone end

"""
    SymmetricMatrixShape(
        side_dimension::Int;
        needs_adjoint_dual::Bool = false,
    )

The shape object for a symmetric square matrix of `side_dimension` rows and
columns.

The vectorized form contains the entries of the upper-right triangular part of
the matrix given column by column (or equivalently, the entries of the
lower-left triangular part given row by row).

## `needs_adjoint_dual`

By default, the [`dual_shape`](@ref) of [`SymmetricMatrixShape`](@ref) is also
[`SymmetricMatrixShape`](@ref). This is true for cases such as a
`LinearAlgebra.Symmetric` matrix in [`PSDCone`](@ref).

However, JuMP also supports `LinearAlgebra.Symmetric` matrix in `Zeros`, which
is interpreted as an element-wise equality constraint. By exploiting symmetry,
we pass only the upper triangle of the equality constraints. This works for the
primal, but it leads to a factor of 2 difference in the off-diagonal dual
elements. (The dual value of the `(i, j)` element in the triangle formulation
should be divided by 2 when spread across the `(i, j)` and `(j, i)` elements in
the square matrix formulation.) If the constraint has this dual inconsistency,
set `needs_adjoint_dual = true`.
"""
struct SymmetricMatrixShape <: AbstractShape
    side_dimension::Int
    needs_adjoint_dual::Bool

    function SymmetricMatrixShape(
        side_dimension::Int;
        needs_adjoint_dual::Bool = false,
    )
        return new(side_dimension, needs_adjoint_dual)
    end
end

function dual_shape(s::SymmetricMatrixShape)
    if s.needs_adjoint_dual
        return SymmetricMatrixAdjointShape(s.side_dimension)
    end
    return s
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
    return LinearAlgebra.Symmetric(matrix)
end

function reshape_set(
    ::MOI.PositiveSemidefiniteConeTriangle,
    ::SymmetricMatrixShape,
)
    return PSDCone()
end

"""
    SymmetricMatrixAdjointShape(side_dimension)

The [`dual_shape`](@ref) of [`SymmetricMatrixShape`](@ref).

This shape is not intended for regular use.
"""
struct SymmetricMatrixAdjointShape <: AbstractShape
    side_dimension::Int
end

function vectorize(matrix::AbstractMatrix, s::SymmetricMatrixAdjointShape)
    n = LinearAlgebra.checksquare(matrix)
    return [((i == j) ? 1 : 2) * matrix[i, j] for j in 1:n for i in 1:j]
end

function reshape_vector(
    v::Vector{T},
    shape::SymmetricMatrixAdjointShape,
) where {T}
    matrix = Matrix{T}(undef, shape.side_dimension, shape.side_dimension)
    k = 0
    for j in 1:shape.side_dimension
        for i in 1:j-1
            k += 1
            matrix[j, i] = matrix[i, j] = 0.5 * v[k]
        end
        k += 1
        matrix[j, j] = v[k]
    end
    return LinearAlgebra.Symmetric(matrix)
end

"""
    triangle_vec(matrix::Matrix)

Return the upper triangle of a matrix concatenated into a vector in the order
required by JuMP and MathOptInterface for `Triangle` sets.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, X[1:3, 1:3], Symmetric);

julia> @variable(model, t)
t

julia> @constraint(model, [t; triangle_vec(X)] in MOI.RootDetConeTriangle(3))
[t, X[1,1], X[1,2], X[2,2], X[1,3], X[2,3], X[3,3]] ∈ MathOptInterface.RootDetConeTriangle(3)
```
"""
function triangle_vec(matrix::AbstractMatrix)
    n = LinearAlgebra.checksquare(matrix)
    return [matrix[i, j] for j in 1:n for i in 1:j]
end

vectorize(matrix, ::SymmetricMatrixShape) = triangle_vec(matrix)

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

function vectorize(matrix, ::SkewSymmetricMatrixShape)
    n = LinearAlgebra.checksquare(matrix)
    return [matrix[i, j] for j in 1:n for i in 1:j-1]
end

"""
    SquareMatrixShape

Shape object for a square matrix of `side_dimension` rows and columns. The
vectorized form contains the entries of the matrix given column by column
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

function vectorize(matrix, shape::SquareMatrixShape)
    return vectorize(Matrix(matrix), shape)
end

# This is a special method because calling `Matrix(matrix)` accesses an undef
# reference.
function vectorize(matrix::LinearAlgebra.UpperTriangular, ::SquareMatrixShape)
    n = LinearAlgebra.checksquare(matrix)
    return [matrix[i, j] for j in 1:n for i in 1:n]
end

# This is a special method because calling `Matrix(matrix)` accesses an undef
# reference.
function vectorize(matrix::LinearAlgebra.LowerTriangular, ::SquareMatrixShape)
    n = LinearAlgebra.checksquare(matrix)
    return [matrix[i, j] for j in 1:n for i in 1:n]
end

function _square_side(error_fn::Function, variables::Matrix)
    n, m = size(variables)
    if n != m
        error_fn("Symmetric variables must be square. Got size ($n, $m).")
    end
    return n
end

function _vectorize_variables(error_fn::Function, matrix::Matrix)
    n = LinearAlgebra.checksquare(matrix)
    for j in 1:n
        for i in 1:j
            if matrix[i, j] != matrix[j, i]
                error_fn(
                    "Non-symmetric bounds, integrality or starting values for symmetric variable.",
                )
            end
        end
    end
    return vectorize(matrix, SymmetricMatrixShape(n))
end

function build_variable(
    error_fn::Function,
    variables::Matrix{<:AbstractVariable},
    ::SymmetricMatrixSpace,
)
    n = _square_side(error_fn, variables)
    set = MOI.Reals(MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n)))
    shape = SymmetricMatrixShape(n)
    return VariablesConstrainedOnCreation(
        _vectorize_variables(error_fn, variables),
        set,
        shape,
    )
end

function build_variable(
    error_fn::Function,
    variables::Matrix{<:AbstractVariable},
    ::SkewSymmetricMatrixSpace,
)
    n = _square_side(error_fn, variables)
    set = MOI.Reals(div(n^2 - n, 2))
    shape = SkewSymmetricMatrixShape(n)
    return VariablesConstrainedOnCreation(
        vectorize(variables, SkewSymmetricMatrixShape(n)),
        set,
        shape,
    )
end

function build_variable(
    error_fn::Function,
    variables::Matrix{<:AbstractVariable},
    ::HermitianMatrixSpace,
)
    n = _square_side(error_fn, variables)
    set = MOI.Reals(
        MOI.dimension(MOI.HermitianPositiveSemidefiniteConeTriangle(n)),
    )
    shape = HermitianMatrixShape(n)
    return VariablesConstrainedOnCreation(
        _vectorize_complex_variables(error_fn, variables),
        set,
        shape,
    )
end

function build_variable(
    error_fn::Function,
    variables::Matrix{<:AbstractVariable},
    ::PSDCone,
)
    n = _square_side(error_fn, variables)
    set = MOI.PositiveSemidefiniteConeTriangle(n)
    return build_variable(error_fn, variables, set)
end

function value(
    Q::LinearAlgebra.Symmetric{V,Matrix{V}},
) where {V<:AbstractVariableRef}
    return LinearAlgebra.Symmetric(
        value.(LinearAlgebra.parent(Q)),
        LinearAlgebra.sym_uplo(Q.uplo),
    )
end

function build_constraint(
    error_fn::Function,
    Q::LinearAlgebra.Symmetric{V,M},
    ::PSDCone,
) where {V<:AbstractJuMPScalar,M<:AbstractMatrix{V}}
    n = LinearAlgebra.checksquare(Q)
    return build_constraint(
        error_fn,
        Q,
        MOI.PositiveSemidefiniteConeTriangle(n),
    )
end

function build_constraint(
    error_fn::Function,
    Q::AbstractMatrix{<:AbstractJuMPScalar},
    ::PSDCone,
)
    n = LinearAlgebra.checksquare(Q)
    return build_constraint(error_fn, Q, MOI.PositiveSemidefiniteConeSquare(n))
end

"""
    HermitianPSDCone

Hermitian positive semidefinite cone object that can be used to create a
Hermitian positive semidefinite square matrix in the [`@variable`](@ref)
and [`@constraint`](@ref) macros.

## Example

Consider the following example:
```jldoctest
julia> model = Model();

julia> @variable(model, H[1:3, 1:3] in HermitianPSDCone())
3×3 LinearAlgebra.Hermitian{GenericAffExpr{ComplexF64, VariableRef}, Matrix{GenericAffExpr{ComplexF64, VariableRef}}}:
 real(H[1,1])                    …  real(H[1,3]) + imag(H[1,3]) im
 real(H[1,2]) - imag(H[1,2]) im     real(H[2,3]) + imag(H[2,3]) im
 real(H[1,3]) - imag(H[1,3]) im     real(H[3,3])

julia> all_variables(model)
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
 [real(H[1,1]), real(H[1,2]), real(H[2,2]), real(H[1,3]), real(H[2,3]), real(H[3,3]), imag(H[1,2]), imag(H[1,3]), imag(H[2,3])] ∈ MathOptInterface.HermitianPositiveSemidefiniteConeTriangle(3)
```
We see in the output of the last commands that 9 real variables were created.
The matrix `H` constrains affine expressions in terms of these 9 variables that
parametrize a Hermitian matrix.
"""
struct HermitianPSDCone end

"""
    HermitianMatrixShape(
        side_dimension::Int;
        needs_adjoint_dual::Bool = false,
    )

The shape object for a Hermitian square matrix of `side_dimension` rows and
columns.

The vectorized form corresponds to
[`MOI.HermitianPositiveSemidefiniteConeTriangle`](@ref).

## `needs_adjoint_dual`

By default, the [`dual_shape`](@ref) of [`HermitianMatrixShape`](@ref) is also
[`HermitianMatrixShape`](@ref). This is true for cases such as a
`LinearAlgebra.Hermitian` matrix in [`HermitianPSDCone`](@ref).

However, JuMP also supports `LinearAlgebra.Hermitian` matrix in `Zeros`, which
is interpreted as an element-wise equality constraint. By exploiting symmetry,
we pass only the upper triangle of the equality constraints. This works for the
primal, but it leads to a factor of 2 difference in the off-diagonal dual
elements. (The dual value of the `(i, j)` element in the triangle formulation
should be divided by 2 when spread across the `(i, j)` and `(j, i)` elements in
the square matrix formulation.) If the constraint has this dual inconsistency,
set `needs_adjoint_dual = true`.
"""
struct HermitianMatrixShape <: AbstractShape
    side_dimension::Int
    needs_adjoint_dual::Bool

    function HermitianMatrixShape(
        side_dimension::Int;
        needs_adjoint_dual::Bool = false,
    )
        return new(side_dimension, needs_adjoint_dual)
    end
end

function dual_shape(s::HermitianMatrixShape)
    if s.needs_adjoint_dual
        return HermitianMatrixAdjointShape(s.side_dimension)
    end
    return s
end

function vectorize(matrix, ::HermitianMatrixShape)
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
    return LinearAlgebra.Hermitian(matrix)
end

"""
    HermitianMatrixAdjointShape(side_dimension)

The [`dual_shape`](@ref) of [`HermitianMatrixShape`](@ref).

This shape is not intended for regular use.
"""
struct HermitianMatrixAdjointShape <: AbstractShape
    side_dimension::Int
end

function vectorize(matrix, ::HermitianMatrixAdjointShape)
    n = LinearAlgebra.checksquare(matrix)
    real_shape = SymmetricMatrixAdjointShape(n)
    imag_shape = SymmetricMatrixShape(n - 1)
    return vcat(
        vectorize(_real.(matrix), real_shape),
        vectorize(2 * _imag.(matrix[1:(end-1), 2:end]), imag_shape),
    )
end

function reshape_vector(
    v::Vector{T},
    shape::HermitianMatrixAdjointShape,
) where {T}
    NewType = _MA.promote_operation(_MA.add_mul, T, Complex{Bool}, T)
    n = shape.side_dimension
    matrix = Matrix{NewType}(undef, n, n)
    real_k = 0
    imag_k = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    for j in 1:n
        for i in 1:(j-1)
            real_k += 1
            imag_k += 1
            matrix[i, j] = (v[real_k] + im * v[imag_k]) / 2
            matrix[j, i] = (v[real_k] - im * v[imag_k]) / 2
        end
        real_k += 1
        matrix[j, j] = v[real_k]
    end
    return LinearAlgebra.Hermitian(matrix)
end

function _vectorize_complex_variables(error_fn::Function, matrix::Matrix)
    if any(_is_binary, matrix) || any(_is_integer, matrix)
        # We would then need to fix the imaginary value to zero. Let's wait to
        # see if there is need for such complication first.
        error_fn(
            "Binary or integer variables in a Hermitian matrix are not supported.",
        )
    end
    n = LinearAlgebra.checksquare(matrix)
    for j in 1:n
        if !_isreal(matrix[j, j])
            error_fn(
                "Non-real bounds or starting values for diagonal of Hermitian variable.",
            )
        end
        for i in 1:j
            if matrix[i, j] != _conj(matrix[j, i])
                error_fn(
                    "Non-conjugate bounds or starting values for Hermitian variable.",
                )
            end
        end
    end
    return vectorize(matrix, HermitianMatrixShape(n))
end

function build_variable(
    error_fn::Function,
    variables::Matrix{<:AbstractVariable},
    ::HermitianPSDCone,
)
    n = _square_side(error_fn, variables)
    set = MOI.HermitianPositiveSemidefiniteConeTriangle(n)
    shape = HermitianMatrixShape(n)
    return VariablesConstrainedOnCreation(
        _vectorize_complex_variables(error_fn, variables),
        set,
        shape,
    )
end

function build_constraint(
    ::Function,
    Q::LinearAlgebra.Hermitian{V,M},
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

function build_constraint(
    error_fn::Function,
    Q::AbstractMatrix{<:AbstractJuMPScalar},
    cone::HermitianPSDCone,
)
    return error_fn(
        "Unable to add matrix in HermitianPSDCone because the matrix is " *
        "not a subtype of `LinearAlgebra.Hermitian`. To fix, wrap the matrix " *
        "`H` in `LinearAlgebra.Hermitian(H)`.",
    )
end

function build_constraint(
    error_fn::Function,
    H::LinearAlgebra.Hermitian,
    ::Zeros,
)
    n = LinearAlgebra.checksquare(H)
    shape = HermitianMatrixShape(n; needs_adjoint_dual = true)
    x = vectorize(H, shape)
    return VectorConstraint(x, MOI.Zeros(length(x)), shape)
end

# If we have a real-valued Hermitian matrix, then it is actually Symmetric, and
# not Complex-valued Hermitian.
function build_constraint(
    error_fn::Function,
    H::LinearAlgebra.Hermitian{V},
    set::Zeros,
) where {
    V<:Union{
        GenericVariableRef{<:Real},
        GenericAffExpr{<:Real},
        GenericQuadExpr{<:Real},
    },
}
    return build_constraint(error_fn, LinearAlgebra.Symmetric(H), set)
end

reshape_set(s::MOI.Zeros, ::HermitianMatrixShape) = Zeros()

function build_constraint(error_fn::Function, ::AbstractMatrix, ::Nonnegatives)
    return error_fn(
        "Unsupported matrix in vector-valued set. Did you mean to use the " *
        "broadcasting syntax `.>=` instead? Alternatively, perhaps you are " *
        "missing a set argument like `@constraint(model, X >= 0, PSDCone())` " *
        "or `@constraint(model, X >= 0, HermitianPSDCone())`.",
    )
end

function build_constraint(error_fn::Function, ::AbstractMatrix, ::Nonpositives)
    return error_fn(
        "Unsupported matrix in vector-valued set. Did you mean to use the " *
        "broadcasting syntax `.<=` instead? Alternatively, perhaps you are " *
        "missing a set argument like `@constraint(model, X <= 0, PSDCone())` " *
        "or `@constraint(model, X <= 0, HermitianPSDCone())`.",
    )
end

function build_constraint(error_fn::Function, ::AbstractMatrix, ::Zeros)
    return error_fn(
        "Unsupported matrix in vector-valued set. Did you mean to use the " *
        "broadcasting syntax `.==` for element-wise equality? Alternatively, " *
        "this syntax is supported in the special case that the matrices are " *
        "`Array`, `LinearAlgebra.Symmetric`, or `LinearAlgebra.Hermitian`.",
    )
end

function build_constraint(
    error_fn::Function,
    Q::LinearAlgebra.Symmetric{V,M},
    set::MOI.AbstractSymmetricMatrixSetTriangle,
) where {V<:AbstractJuMPScalar,M<:AbstractMatrix{V}}
    n = LinearAlgebra.checksquare(Q)
    shape = SymmetricMatrixShape(n)
    return VectorConstraint(vectorize(Q, shape), set, shape)
end

function build_constraint(
    error_fn::Function,
    Q::AbstractMatrix{<:AbstractJuMPScalar},
    set::MOI.AbstractSymmetricMatrixSetSquare,
)
    n = LinearAlgebra.checksquare(Q)
    shape = SquareMatrixShape(n)
    return VectorConstraint(vectorize(Q, shape), set, shape)
end

function build_constraint(
    error_fn::Function,
    f::AbstractMatrix{<:AbstractJuMPScalar},
    ::Nonnegatives,
    extra::Union{
        MOI.AbstractSymmetricMatrixSetTriangle,
        MOI.AbstractSymmetricMatrixSetSquare,
        PSDCone,
        HermitianPSDCone,
    },
)
    return build_constraint(error_fn, f, extra)
end

function build_constraint(
    error_fn::Function,
    f::AbstractMatrix{<:AbstractJuMPScalar},
    ::Nonpositives,
    extra::Union{
        MOI.AbstractSymmetricMatrixSetTriangle,
        MOI.AbstractSymmetricMatrixSetSquare,
        PSDCone,
        HermitianPSDCone,
    },
)
    new_f = _MA.operate!!(*, -1, f)
    return build_constraint(error_fn, new_f, extra)
end

function build_variable(
    error_fn::Function,
    variables::Matrix{<:AbstractVariable},
    set::MOI.AbstractSymmetricMatrixSetTriangle,
)
    n = _square_side(error_fn, variables)
    x = _vectorize_variables(error_fn, variables)
    return VariablesConstrainedOnCreation(x, set, SymmetricMatrixShape(n))
end

moi_set(::Nonnegatives, dim::Int) = MOI.Nonnegatives(dim)
moi_set(::Nonpositives, dim::Int) = MOI.Nonpositives(dim)
moi_set(::Zeros, dim::Int) = MOI.Zeros(dim)

function _shape_for_orthants(f::LinearAlgebra.Symmetric)
    n = LinearAlgebra.checksquare(f)
    return SymmetricMatrixShape(n; needs_adjoint_dual = true)
end

reshape_set(::MOI.Nonnegatives, ::SymmetricMatrixShape) = Nonnegatives()
reshape_set(::MOI.Nonpositives, ::SymmetricMatrixShape) = Nonpositives()
reshape_set(::MOI.Zeros, ::SymmetricMatrixShape) = Zeros()

_shape_for_orthants(f::Array) = ArrayShape(size(f))

reshape_set(::MOI.Nonnegatives, ::ArrayShape) = Nonnegatives()
reshape_set(::MOI.Nonpositives, ::ArrayShape) = Nonpositives()
reshape_set(::MOI.Zeros, ::ArrayShape) = Zeros()

# We use an @eval loop because a `Union` introduces ambiguities.
for S in (Nonnegatives, Nonpositives, Zeros)
    for F in (Array, LinearAlgebra.Symmetric)
        @eval begin
            function build_constraint(error_fn::Function, f::$F, set::$S)
                s = _shape_for_orthants(f)
                x = vectorize(f, s)
                return VectorConstraint(x, moi_set(set, length(x)), s)
            end

            function build_constraint(
                error_fn::Function,
                f::$F,
                ::Nonnegatives,
                set::$S,
            )
                return build_constraint(error_fn, f, set)
            end

            function build_constraint(
                error_fn::Function,
                ::$F,
                ::Nonpositives,
                set::$S,
            )
                return error_fn(
                    "The syntax `x <= y, $set` not supported. Use `y >= x, $set` instead.",
                )
            end
        end
    end
end
