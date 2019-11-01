function build_indicator_constraint(
    _error::Function, variable::JuMP.AbstractVariableRef, constraint::JuMP.ScalarConstraint, ::Type{MOI.IndicatorSet{A}}) where A

    set = MOI.IndicatorSet{A}(moi_set(constraint))
    return VectorConstraint([variable, jump_function(constraint)], set)
end
function _indicator_variable_set(::Function, variable::Symbol)
    return variable, MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}
end
function _indicator_variable_set(_error::Function, expr::Expr)
    if expr.args[1] == :¬ || expr.args[1] == :!
        if length(expr.args) != 2
            _error("Invalid binary variable expression `$(expr)` for indicator constraint.")
        end
        return expr.args[2], MOI.IndicatorSet{MOI.ACTIVATE_ON_ZERO}
    else
        return expr, MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}
    end
end
function parse_one_operator_constraint(
    _error::Function, vectorized::Bool, ::Union{Val{:(=>)}, Val{:⇒}}, lhs, rhs)

    variable, S = _indicator_variable_set(_error, lhs)
    if !(rhs isa Expr)
        _error("Invalid right-hand side `$(rhs)`.")
    end
    rhs_vectorized, rhs_parsecode, rhs_buildcall = parse_constraint(_error, rhs.args...)
    if vectorized != rhs_vectorized
        _error("Inconsistent use of `.` in symbols to indicate vectorization.")
    end
    if vectorized
        buildcall = :(build_indicator_constraint.($_error, $(esc(variable)), $rhs_buildcall, $S))
    else
        buildcall = :(build_indicator_constraint($_error, $(esc(variable)), $rhs_buildcall, $S))
    end
    return rhs_parsecode, buildcall
end
