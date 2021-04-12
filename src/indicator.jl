#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

# This file extends JuMP to indicator constraints. It is a good example of how
# JuMP can be extended.

function _build_indicator_constraint(
    _error::Function,
    variable::AbstractVariableRef,
    constraint::ScalarConstraint,
    ::Type{MOI.IndicatorSet{A}},
) where {A}
    set = MOI.IndicatorSet{A}(moi_set(constraint))
    return VectorConstraint([variable, jump_function(constraint)], set)
end
function _indicator_variable_set(::Function, variable::Symbol)
    return variable, MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}
end
function _indicator_variable_set(_error::Function, expr::Expr)
    if expr.args[1] == :¬ || expr.args[1] == :!
        if length(expr.args) != 2
            _error(
                "Invalid binary variable expression `$(expr)` for indicator constraint.",
            )
        end
        return expr.args[2], MOI.IndicatorSet{MOI.ACTIVATE_ON_ZERO}
    else
        return expr, MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}
    end
end
function parse_one_operator_constraint(
    _error::Function,
    vectorized::Bool,
    ::Union{Val{:(=>)},Val{:⇒}},
    lhs,
    rhs,
)
    variable, S = _indicator_variable_set(_error, lhs)
    if !isexpr(rhs, :braces) || length(rhs.args) != 1
        _error(
            "Invalid right-hand side `$(rhs)` of indicator constraint. Expected constraint surrounded by `{` and `}`.",
        )
    end
    rhs_con = rhs.args[1]
    rhs_vectorized, rhs_parsecode, rhs_buildcall =
        parse_constraint_expr(_error, rhs_con)
    if vectorized != rhs_vectorized
        _error("Inconsistent use of `.` in symbols to indicate vectorization.")
    end
    if vectorized
        buildcall = :(
            _build_indicator_constraint.(
                $_error,
                $(esc(variable)),
                $rhs_buildcall,
                $S,
            )
        )
    else
        buildcall = :(_build_indicator_constraint(
            $_error,
            $(esc(variable)),
            $rhs_buildcall,
            $S,
        ))
    end
    return rhs_parsecode, buildcall
end

function constraint_string(
    print_mode,
    constraint::VectorConstraint{F,<:MOI.IndicatorSet{A}},
) where {F,A}
    # TODO Implement pretty IJulia printing
    var_str = function_string(print_mode, constraint.func[1])
    if A == MOI.ACTIVATE_ON_ZERO
        var_str = "!" * var_str
    end
    con = ScalarConstraint(constraint.func[2], constraint.set.set)
    con_str = constraint_string(print_mode, con)
    return var_str * " => {" * con_str * "}"
end
