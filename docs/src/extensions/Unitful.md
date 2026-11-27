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

### Variables with units

Create a variable with units by passing the unit as a positional argument to
[`@variable`](@ref):

```jldoctest ext_unitful
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 0u"m/s", u"m/s")
(x) m s^-1

julia> @variable(model, y >= 0u"m/s", u"ft/s")
(y) ft s^-1
```

The bounds and start values must have units, even if they are zero:

```jldoctest ext_unitful
julia> @variable(model, z >= 0, u"m/s")
ERROR: At none:1: `@variable(model, z >= 0, u"m/s")`:

Unable to use `0` as a bound or start value because it has no units.

Use a quantity with the same dimension as `m s^-1`.

Stacktrace:
[...]
```

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

To improve the scaling of the problem, pass a quantity instead of a unit:

```jldoctest ext_unitful_scale
julia> using JuMP, Unitful

julia> model = Model();

julia> @variable(model, 0u"m/s" <= w <= 100u"m/s", 10u"m/s")
(10 w) m s^-1

julia> upper_bound(only(all_variables(model)))
10.0
```

Here, `w` is the expression `10.0u"m/s" * w`, where the variable `w` without
units is between `0` and `10`.

### Expressions with units

The unit of an expression is the unit of its coefficient type. If all its terms
have the same unit, the expression keeps this unit:

```jldoctest ext_unitful
julia> y + 1u"ft/s"
(y + 1) ft s^-1
```

If the terms have different units of the same dimension, Unitful converts them
to the preferred unit of this dimension, regardless of the order of the terms.
See `Unitful.upreferred` and `Unitful.preferunits` for details. By default, the
preferred units are the SI units:

```jldoctest ext_unitful
julia> y + x
(0.3048 y + x) m s^-1

julia> Unitful.upreferred(u"ft/s")
m s^-1
```

### Constraints with units

When adding a constraint, the units are converted to the unit of the
expression `lhs - rhs`, as defined in [Expressions with units](@ref), and then
removed before the constraint is passed to the solver. The unit is stored in
the [`AbstractShape`](@ref) of the constraint:

```jldoctest ext_unitful
julia> @constraint(model, c1, x + y <= 72u"km/hr")
c1 : (x + 0.3048 y) m s^-1 ≤ 20 m s^-1
```

Pass a unit as a positional argument to [`@constraint`](@ref) to choose the unit
of the constraint, for example, for better scaling:

```jldoctest ext_unitful
julia> @constraint(model, c2, x + y <= 72u"km/hr", u"km/hr")
c2 : (3.6 x + 1.09728 y) km hr^-1 ≤ 72 km hr^-1
```

Note that a number like `0.5` has no unit, so the units of both sides of the following
constraint are compatible:

```jldoctest ext_unitful
julia> @constraint(model, c3, x <= 0.5 * y)
c3 : (x - 0.1524 y) m s^-1 ≤ 0 m s^-1
```

Constraints with incompatible units throw an error:

```jldoctest ext_unitful
julia> @constraint(model, x <= 1u"m")
ERROR: DimensionError: 0.0 m s^-1 and 0 m are not dimensionally compatible.
[...]
```

### Objective with units

When setting the objective, the units are converted to the unit of the
expression and then removed:

```jldoctest ext_unitful
julia> @objective(model, Max, 2 * x + y)
(2 x + 0.3048 y) m s^-1

julia> objective_function(model)
2 x + 0.3048 y
```

### Solutions with units

Querying the value of an expression or a constraint with units returns a value
with units:

```jldoctest ext_unitful
julia> optimize!(model)

julia> value(x)
6.666666666666667 m s^-1

julia> value(y)
43.744531933508306 ft s^-1

julia> value(x + y)
20.0 m s^-1

julia> value(c2)
72.0 km hr^-1
```

The objective is passed to the solver without units, so the
[`objective_value`](@ref) does not have units:

```jldoctest ext_unitful
julia> objective_value(model)
26.666666666666664
```

Because the objective has no units, the [`dual`](@ref) of a constraint has the
inverse unit of the constraint:

```jldoctest ext_unitful
julia> dual(c1)
-1.3333333333333333 s m^-1

julia> dual(c3)
-0.6666666666666666 s m^-1
```
