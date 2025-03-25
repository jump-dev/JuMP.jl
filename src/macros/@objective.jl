#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    @objective(model::GenericModel, sense, func)

Set the objective sense to `sense` and objective function to `func`.

## `sense`

The objective sense must be either be the literals `Min` or `Max`, or one of the
three [`MOI.OptimizationSense`](@ref) enum values ([`MIN_SENSE`](@ref),
[`MAX_SENSE`](@ref), or [`FEASIBILITY_SENSE`](@ref)).

In order to set the sense programmatically, that is, when `sense` is a Julia
variable whose value is the sense, you must use a [`MOI.OptimizationSense`](@ref).

## FEASIBILITY_SENSE

[`FEASIBILITY_SENSE`](@ref) implies that there is no objective function.
Therefore, you should not set `sense` to [`FEASIBILITY_SENSE`](@ref) with a
non-zero `func`.

To reset the model to [`FEASIBILITY_SENSE`](@ref), do
`@objective(model, FEASIBILITY_SENSE, 0)`, or use [`set_objective_sense`](@ref):
`set_objective_sense(model, FEASIBILITY_SENSE)`.

## Example

Minimize the value of the variable `x`, do:
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @objective(model, Min, x)
x
```

Maximize the value of the affine expression `2x - 1`:
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @objective(model, Max, 2x - 1)
2 x - 1
```

Set the objective sense programmatically:
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> sense = MIN_SENSE
MIN_SENSE::OptimizationSense = 0

julia> @objective(model, sense, x^2 - 2x + 1)
x² - 2 x + 1
```

Remove an objective function:
```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0);

julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> print(model)
Min 2 x + 1
Subject to
 x ≥ 0

julia> @objective(model, FEASIBILITY_SENSE, 0)
0

julia> print(model)
Feasibility
Subject to
 x ≥ 0
```
"""
macro objective(input_args...)
    error_fn = Containers.build_error_fn(:objective, input_args, __source__)
    args, kwargs = Containers.parse_macro_arguments(
        error_fn,
        input_args;
        num_positional_args = 3,
        valid_kwargs = Symbol[],
    )
    esc_model = esc(args[1])
    sense = _parse_moi_sense(error_fn, args[2])
    expr, parse_code = _rewrite_expression(args[3])
    code = quote
        $parse_code
        # Don't leak a `_MA.Zero` if the objective expression is an empty
        # summation, or other structure that returns `_MA.Zero()`.
        $expr = _replace_zero($esc_model, $expr)
        set_objective($esc_model, $sense, $expr)
        $expr
    end
    return _finalize_macro(esc_model, code, __source__)
end

function _parse_moi_sense(error_fn::Function, sense)
    if sense == :Min
        return MIN_SENSE
    elseif sense == :Max
        return MAX_SENSE
    end
    return :(_moi_sense($error_fn, $(esc(sense))))
end

function _moi_sense(error_fn::Function, sense)
    return error_fn(
        "unexpected sense `$sense`. The sense must be an " *
        "`::MOI.OptimizatonSense`, or the symbol `:Min` or `:Max`.",
    )
end

_moi_sense(::Function, sense::MOI.OptimizationSense) = sense
