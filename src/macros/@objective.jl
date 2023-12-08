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
macro objective(model, args...)
    function error_fn(str...)
        return _macro_error(:objective, (model, args...), __source__, str...)
    end

    # We don't overwrite `model` as it is used in `error_fn`
    esc_model = esc(model)
    if length(args) != 2
        # Either just an objective sense, or just an expression.
        error_fn(
            "needs three arguments: model, objective sense (Max or Min) and expression.",
        )
    end
    sense, x = args
    sense_expr = _moi_sense(error_fn, sense)
    newaff, parsecode = _rewrite_expression(x)
    code = quote
        $parsecode
        # Don't leak a `_MA.Zero` if the objective expression is an empty
        # summation, or other structure that returns `_MA.Zero()`.
        $newaff = _replace_zero($esc_model, $newaff)
        set_objective($esc_model, $sense_expr, $newaff)
        $newaff
    end
    return _finalize_macro(esc_model, code, __source__)
end

"""
    _moi_sense(error_fn::Function, sense)

Return an expression whose value is an `MOI.OptimizationSense` corresponding
to `sense`. Sense is either the symbol `:Min` or `:Max`, corresponding
respectively to `MIN_SENSE` and `MAX_SENSE` or it is another symbol,
which should be the name of a variable or expression whose value is an
`MOI.OptimizationSense`.
In the last case, the expression throws an error using the `error_fn`
function in case the value is not an `MOI.OptimizationSense`.
"""
function _moi_sense(error_fn::Function, sense)
    expr = if sense == :Min
        MIN_SENSE
    elseif sense == :Max
        MAX_SENSE
    else
        esc(sense)
    end
    return :(_throw_error_for_invalid_sense($error_fn, $expr))
end

function _throw_error_for_invalid_sense(error_fn::Function, sense)
    return error_fn(
        "unexpected sense `$sense`. The sense must be an " *
        "`::MOI.OptimizatonSense`, or the symbol `:Min` or `:Max`.",
    )
end

_throw_error_for_invalid_sense(::Function, sense::MOI.OptimizationSense) = sense
