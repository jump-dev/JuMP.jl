#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function _build_inequality_constraint(
    error_fn::Function,
    vectorized::Bool,
    lhs::VariableRef,
    rhs::VariableRef,
)
    @assert !vectorized
    set = MOI.AllDifferent(2)
    return VectorConstraint([lhs; rhs], set)
end

function _build_inequality_constraint(
    error_fn::Function,
    vectorized::Bool,
    lhs::Vector{VariableRef},
    rhs::Vector{VariableRef},
)
    if !vectorized
        error_fn(
            "Ineqality operator with vector operands must be explicitly " *
            "vectorized, use `.!=` instead of `!=`.",
        )
    end
    if length(lhs) != length(rhs)
        error_fn("Operand length mismatch, $(length(lhs)) vs $(length(rhs)).")
    end
    lhs = _desparsify(lhs)
    rhs = _desparsify(rhs)
    return _build_inequality_constraint.(error_fn, false, lhs, rhs)
end

function _build_inequality_constraint(error_fn::Function, ::Bool, lhs, rhs)
    return error_fn(
        "Unsupported form of inequality constraint. The left- and right-hand " *
        "sides must both be decision variables.",
    )
end

function parse_constraint_call(
    error_fn::Function,
    vectorized::Bool,
    ::Val{:(!=)},
    lhs,
    rhs,
)
    build_call = Expr(
        :call,
        :_build_inequality_constraint,
        error_fn,
        vectorized,
        esc(lhs),
        esc(rhs),
    )
    return nothing, build_call
end

function constraint_string(
    print_mode,
    constraint::VectorConstraint{F,<:MOI.AllDifferent},
) where {F}
    set = constraint.set
    if set.dimension == 2
        ineq_sym = JuMP._math_symbol(print_mode, :(!=))
        lhs = function_string(print_mode, constraint.func[1])
        rhs = function_string(print_mode, constraint.func[2])
        return string(lhs, " $ineq_sym ", rhs)
    end

    # FIXME: can we just fallback to the generic handling here?
    ops = [function_string(print_mode, op) for op in constraint.func[1:end]]
    in_sym = JuMP._math_symbol(print_mode, :in)
    return string("[", join(ops, ", "), "] $in_sym $set")
end
