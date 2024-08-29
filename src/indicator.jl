#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

# This file extends JuMP to indicator constraints. It is a good example of how
# JuMP can be extended.

function _build_indicator_constraint(
    error_fn::Function,
    variable::AbstractVariableRef,
    constraint::ScalarConstraint,
    ::Type{MOI.Indicator{A}},
) where {A}
    set = MOI.Indicator{A}(moi_set(constraint))
    return VectorConstraint([variable, jump_function(constraint)], set)
end

function _build_indicator_constraint(
    error_fn::Function,
    lhs::F,
    ::ScalarConstraint,
    ::Type{<:MOI.Indicator},
) where {F}
    return error_fn(
        "unable to build indicator constraint with the left-hand side term " *
        "`($lhs)::$F`. The left-hand side must be a binary decision variable.",
    )
end

function _indicator_variable_set(::Function, variable::Symbol)
    return variable, MOI.Indicator{MOI.ACTIVATE_ON_ONE}
end

function _indicator_variable_set(error_fn::Function, expr::Expr)
    if expr.args[1] == :¬ || expr.args[1] == :!
        if length(expr.args) != 2
            error_fn(
                "Invalid binary variable expression `$(expr)` for indicator constraint.",
            )
        end
        return expr.args[2], MOI.Indicator{MOI.ACTIVATE_ON_ZERO}
    else
        return expr, MOI.Indicator{MOI.ACTIVATE_ON_ONE}
    end
end

function parse_constraint_head(error_fn::Function, ::Val{:(-->)}, lhs, rhs)
    code, call = parse_constraint_call(error_fn, false, Val(:(=>)), lhs, rhs)
    return false, code, call
end

function parse_constraint_call(
    error_fn::Function,
    vectorized::Bool,
    ::Union{Val{:(=>)},Val{:⇒}},
    lhs,
    rhs,
)
    variable, S = _indicator_variable_set(error_fn, lhs)
    if !Meta.isexpr(rhs, :braces) || length(rhs.args) != 1
        error_fn(
            "Invalid right-hand side `$(rhs)` of indicator constraint. Expected constraint surrounded by `{` and `}`.",
        )
    end
    rhs_vectorized, rhs_parsecode, rhs_build_call =
        parse_constraint(error_fn, rhs.args[1])
    if vectorized != rhs_vectorized
        error_fn(
            "Inconsistent use of `.` in symbols to indicate vectorization.",
        )
    end
    f, lhs_parse_code = _rewrite_expression(variable)
    parsecode = quote
        $rhs_parsecode
        $lhs_parse_code
    end
    build_call = if vectorized
        :(_build_indicator_constraint.($error_fn, $f, $rhs_build_call, $S))
    else
        :(_build_indicator_constraint($error_fn, $f, $rhs_build_call, $S))
    end
    return parsecode, build_call
end

function constraint_string(
    print_mode,
    constraint::VectorConstraint{F,<:MOI.Indicator{A}},
) where {F,A}
    # TODO Implement pretty IJulia printing
    var_str = function_string(print_mode, constraint.func[1])
    if A == MOI.ACTIVATE_ON_ZERO
        var_str = "!" * var_str
    end
    con = ScalarConstraint(constraint.func[2], constraint.set.set)
    con_str = constraint_string(print_mode, con)
    return var_str * " --> {" * con_str * "}"
end
