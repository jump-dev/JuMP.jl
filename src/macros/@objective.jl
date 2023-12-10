#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    @objective(model::GenericModel, sense, func)

Set the objective sense to `sense` and objective function to `func`. The
objective sense can be either `Min`, `Max`, `MOI.MIN_SENSE`, `MOI.MAX_SENSE` or
`MOI.FEASIBILITY_SENSE`; see [`MOI.ObjectiveSense`](@ref).

In order to set the sense programmatically, i.e., when `sense` is a Julia
variable whose value is the sense, one of the three `MOI.ObjectiveSense` values
should be used.

## Example

To minimize the value of the variable `x`, do as follows:
```jldoctest @objective
julia> model = Model();

julia> @variable(model, x)
x

julia> @objective(model, Min, x)
x
```

To maximize the value of the affine expression `2x - 1`, do as follows:
```jldoctest @objective
julia> @objective(model, Max, 2x - 1)
2 x - 1
```

To set a quadratic objective and set the objective sense programmatically, do
as follows:
```jldoctest @objective
julia> sense = MIN_SENSE
MIN_SENSE::OptimizationSense = 0

julia> @objective(model, sense, x^2 - 2x + 1)
x² - 2 x + 1
```
"""
macro objective(input_args...)
    error_fn(str...) = _macro_error(:objective, input_args, __source__, str...)
    args, kwargs = Containers.parse_macro_arguments(error_fn, input_args)
    if length(args) != 3
        error_fn("expected 3 positional arguments, got $(length(args)).")
    elseif !isempty(kwargs)
        for key in keys(kwargs)
            error_fn("unsupported keyword argument `$key`.")
        end
    end
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
