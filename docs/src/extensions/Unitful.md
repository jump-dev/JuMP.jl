# Unitful.jl

[Unitful.jl](https://github.com/JuliaPhysics/Unitful.jl) supports physical
units and dimensional analysis.

!!! compat
    Using the Unitful extension with JuMP requires Julia v1.9 or later.

The Unitful extension in JuMP lets you create variables with units, and it
automatically converts the units of the constraints and objective.

## License

Unitful.jl is licensed under the [MIT license](https://github.com/JuliaPhysics/Unitful.jl/blob/master/LICENSE.md).

## Installation

Install Unitful using `Pkg.add`:

```julia
import Pkg
Pkg.add("Unitful")
```

## Use with JuMP

Activate the extension by loading both JuMP and Unitful:

```jldoctest ext_unitful
julia> using JuMP, Unitful

julia> import HiGHS
```

### Variables

Create a variable with units by passing the unit as a positional argument to
[`@variable`](@ref):

```jldoctest ext_unitful
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 0, u"m/s")
(x) m s^-1

julia> @variable(model, y >= 0, u"ft/s")
(y) ft s^-1
```

Bounds without units, such as `0` above, are interpreted in the unit of the
variable.

Similarly to [Complex number support](@ref), `x` is not a [`VariableRef`](@ref);
instead, it is an affine expression with `Unitful.Quantity`-valued
coefficients:

```jldoctest ext_unitful
julia> typeof(x)
GenericAffExpr{Quantity{Float64, 𝐋 𝐓^-1, Unitful.FreeUnits{(m, s^-1), 𝐋 𝐓^-1, nothing}}, VariableRef}
```

Behind the scenes, JuMP has created real-valued variables without units:

```jldoctest ext_unitful
julia> all_variables(model)
2-element Vector{VariableRef}:
 x
 y
```

### Constraints

When adding a constraint, the units are converted to the unit of the
expression, and then removed:

```jldoctest ext_unitful
julia> @constraint(model, c1, x + y <= 72u"km/hr")
c1 : x + 0.3048 y ≤ 20
```

Pass a unit as a positional argument to [`@constraint`](@ref) to choose the unit
of the constraint, for example, for better scaling:

```jldoctest ext_unitful
julia> @constraint(model, c2, x + y <= 72u"km/hr", u"km/hr")
c2 : 3.6 x + 1.09728 y ≤ 72
```

Note that `0.5` has no unit, so the units of both sides of the following
constraint are compatible:

```jldoctest ext_unitful
julia> @constraint(model, c3, x <= 0.5 * y)
c3 : x - 0.1524 y ≤ 0
```

Constraints with incompatible units throw an error:

```jldoctest ext_unitful
julia> @constraint(model, x <= 1u"m")
ERROR: DimensionError: 0.0 m s^-1 and 0 m are not dimensionally compatible.
[...]
```

### Objective

When setting the objective, the units are converted to the unit of the
expression and then removed:

```jldoctest ext_unitful
julia> @objective(model, Max, 2 * x + y)
(2 x + 0.3048 y) m s^-1

julia> objective_function(model)
2 x + 0.3048 y
```

### Solutions

Querying the value of an expression with units returns a value with units:

```jldoctest ext_unitful
julia> optimize!(model)

julia> value(x)
6.666666666666667 m s^-1

julia> value(y)
43.744531933508306 ft s^-1

julia> value(x + y)
20.0 m s^-1
```

However, because the constraints and objective are passed to the solver
without units, solution values such as [`objective_value`](@ref) and
[`dual`](@ref) do not have units:

```jldoctest ext_unitful
julia> objective_value(model)
26.666666666666664
```
