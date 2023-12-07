#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function _build_reified_constraint(error_fn::Function, lhs, rhs)
    set = MOI.Reified(moi_set(rhs))
    return VectorConstraint([lhs; jump_function(rhs)], set)
end

function parse_constraint_call(
    error_fn::Function,
    ::Bool,
    ::Union{Val{:(<-->)},Val{:⟺}},
    lhs,
    rhs,
)
    if !Meta.isexpr(rhs, :braces) || length(rhs.args) != 1
        error_fn(
            "Invalid right-hand side `$(rhs)` of reified constraint. " *
            "Expected constraint surrounded by `{` and `}`.",
        )
    end
    vectorized, parse_code, build_call = parse_constraint(error_fn, rhs.args[1])
    if vectorized
        error_fn("vectorized constraints cannot be used with reification.")
    end
    new_build_call =
        Expr(:call, :_build_reified_constraint, error_fn, esc(lhs), build_call)
    return parse_code, new_build_call
end

function constraint_string(
    print_mode,
    constraint::VectorConstraint{F,<:MOI.Reified},
) where {F}
    lhs = function_string(print_mode, constraint.func[1])
    set = constraint.set.set
    con = if set isa MOI.AbstractVectorSet
        VectorConstraint(constraint.func[2:end], set)
    else
        ScalarConstraint(constraint.func[2], set)
    end
    return string(lhs, " <--> {", constraint_string(print_mode, con), "}")
end
