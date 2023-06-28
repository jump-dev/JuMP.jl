```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
    import HiGHS
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Complex number support

This page explains the complex-valued variables and constraints that JuMP
supports. For a worked-example using these features, read the
[Quantum state discrimination](@ref) tutorial.

## Complex-valued variables

Create a complex-valued variable using [`ComplexPlane`](@ref):

```jldoctest complex_variables
julia> model = Model();

julia> @variable(model, x in ComplexPlane())
real(x) + imag(x) im
```

Note that `x` is not a [`VariableRef`](@ref); instead, it is an affine
expression with `Complex{Float64}`-valued coefficients:

```jldoctest complex_variables
julia> typeof(x)
GenericAffExpr{ComplexF64, VariableRef}
```

Behind the scenes, JuMP has created two real-valued variables, with names
`"real(x)"` and `"imag(x)"`:

```jldoctest complex_variables
julia> all_variables(model)
2-element Vector{VariableRef}:
 real(x)
 imag(x)

julia> name.(all_variables(model))
2-element Vector{String}:
 "real(x)"
 "imag(x)"
```

Use the `real` and `imag` functions on `x` to return a real-valued affine
expression representing each variable:

```jldoctest complex_variables
julia> typeof(real(x))
AffExpr (alias for GenericAffExpr{Float64, GenericVariableRef{Float64}})

julia> typeof(imag(x))
AffExpr (alias for GenericAffExpr{Float64, GenericVariableRef{Float64}})
```

To create an anonymous variable, use the `set` keyword argument:

```jldoctest
julia> model = Model();

julia> x = @variable(model, set = ComplexPlane())
_[1] + _[2] im
```

## Complex-valued variable bounds

Because complex-valued variables lack a total ordering, the definition of a
variable bound for a complex-valued variable is ambiguous. If you pass a real-
or complex-valued argument to keywords such as `lower_bound`, `upper_bound`,
and `start_value`, JuMP will apply the real and imaginary parts to the
associated real-valued variables.

```jldoctest complex_variables
julia> model = Model();

julia> @variable(
           model,
           x in ComplexPlane(),
           lower_bound = 1.0,
           upper_bound = 2.0 + 3.0im,
           start = 4im,
       )
real(x) + imag(x) im

julia> vars = all_variables(model)
2-element Vector{VariableRef}:
 real(x)
 imag(x)

julia> lower_bound.(vars)
2-element Vector{Float64}:
 1.0
 0.0

julia> upper_bound.(vars)
2-element Vector{Float64}:
 2.0
 3.0

julia> start_value.(vars)
2-element Vector{Float64}:
 0.0
 4.0
```

## Complex-valued equality constraints

JuMP reformulates complex-valued equality constraints into two real-valued
constraints: one representing the real part, and one representing the imaginary
part. Thus, complex-valued equality constraints can be solved any solver that
supports the real-valued constraint type.

For example:

```jldoctest
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x[1:2]);

julia> @constraint(model, (1 + 2im) * x[1] + 3 * x[2] == 4 + 5im)
(1 + 2im) x[1] + 3 x[2] = (4 + 5im)

julia> optimize!(model)

julia> value.(x)
2-element Vector{Float64}:
 2.5
 0.5
```

is equivalent to

```jldoctest
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x[1:2]);

julia> @constraint(model, 1 * x[1] + 3 * x[2] == 4)  # real component
x[1] + 3 x[2] = 4

julia> @constraint(model, 2 * x[1] == 5)  # imag component
2 x[1] = 5

julia> optimize!(model)

julia> value.(x)
2-element Vector{Float64}:
 2.5
 0.5
```

This also applies if the variables are complex-valued:

```jldoctest
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x in ComplexPlane());

julia> @constraint(model, (1 + 2im) * x + 3 * x == 4 + 5im)
(4 + 2im) real(x) + (-2 + 4im) imag(x) = (4 + 5im)

julia> optimize!(model)

julia> value(x)
1.3 + 0.6000000000000001im
```

which is equivalent to

```jldoctest
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x_real);

julia> @variable(model, x_imag);

julia> @constraint(model, x_real - 2 * x_imag + 3 * x_real == 4)
4 x_real - 2 x_imag = 4

julia> @constraint(model, x_imag + 2 * x_real + 3 * x_imag == 5)
2 x_real + 4 x_imag = 5

julia> optimize!(model)

julia> value(x_real) + value(x_imag) * im
1.3 + 0.6000000000000001im
```

## Hermitian PSD Cones

JuMP supports creating matrices where are Hermitian.
```jldoctest hermitian_psd_cone
julia> model = Model();

julia> @variable(model, H[1:3, 1:3] in HermitianPSDCone())
3×3 LinearAlgebra.Hermitian{GenericAffExpr{ComplexF64, VariableRef}, Matrix{GenericAffExpr{ComplexF64, VariableRef}}}:
 real(H[1,1])                    …  real(H[1,3]) + imag(H[1,3]) im
 real(H[1,2]) - imag(H[1,2]) im     real(H[2,3]) + imag(H[2,3]) im
 real(H[1,3]) - imag(H[1,3]) im     real(H[3,3])
```

Behind the scenes, JuMP has created nine real-valued decision variables:

```jldoctest hermitian_psd_cone
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
```

and a `Vector{VariableRef}-in-MOI.HermitianPositiveSemidefiniteConeTriangle`
constraint:

```jldoctest hermitian_psd_cone
julia> num_constraints(model, Vector{VariableRef}, MOI.HermitianPositiveSemidefiniteConeTriangle)
1
```

The [`MOI.HermitianPositiveSemidefiniteConeTriangle`](@ref) set can be
efficiently bridged to [`MOI.PositiveSemidefiniteConeTriangle`](@ref), so it can
be solved by any solver that supports PSD constraints.

Each element of `H` is an affine expression with `Complex{Float64}`-valued
coefficients:

```jldoctest hermitian_psd_cone
julia> typeof(H[1, 1])
GenericAffExpr{ComplexF64, VariableRef}

julia> typeof(H[2, 1])
GenericAffExpr{ComplexF64, VariableRef}
```

## Hermitian PSD constraints

The [`HermitianPSDCone`](@ref) can also be used in the [`@constraint`](@ref)
macro:
```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> import LinearAlgebra

julia> H = LinearAlgebra.Hermitian([x[1] 1im; -1im -x[2]])
2×2 LinearAlgebra.Hermitian{GenericAffExpr{ComplexF64, VariableRef}, Matrix{GenericAffExpr{ComplexF64, VariableRef}}}:
 x[1]  im
 -im   -x[2]

julia> @constraint(model, H in HermitianPSDCone())
[x[1]  im;
 -im   -x[2]] ∈ HermitianPSDCone()
```

!!! note
    The matrix `H` in `H in HermitianPSDCone()` must be a `LinearAlgebra.Hermitian`
    matrix type. A `build_constraint` error will be thrown if the matrix is
    a different matrix type.
