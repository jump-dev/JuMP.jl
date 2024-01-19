#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestInequality

using JuMP
using Test

include(joinpath(@__DIR__, "utilities.jl"))

function test_inequality_two_int_scalars()
    model = Model()
    @variable(model, -4 <= x <= 4, Int)
    @variable(model, -4 <= y <= 4, Int)
    c = @constraint(model, x != y)
    ineq_sym = JuMP._math_symbol(MIME("text/plain"), :(!=))
    set = MOI.AllDifferent(2)
    @test sprint(show, c) == "x $ineq_sym y"
    obj = constraint_object(c)
    @test obj.func == [x; y]
    @test obj.set == set
    return
end

function test_inequality_latex()
    model = Model()
    @variable(model, -4 <= x <= 4, Int)
    @variable(model, -4 <= y <= 4, Int)
    c = @constraint(model, x != y)
    set = MOI.AllDifferent(2)
    @test sprint(io -> show(io, MIME("text/latex"), c)) == "\$\$ x \\neq y \$\$"
    obj = constraint_object(c)
    @test obj.func == [x; y]
    @test obj.set == set
    return
end

function test_inequality_two_bin_scalars()
    model = Model()
    @variable(model, x, Bin)
    @variable(model, y, Bin)
    c = @constraint(model, x != y)
    ineq_sym = JuMP._math_symbol(MIME("text/plain"), :(!=))
    set = MOI.AllDifferent(2)
    @test sprint(show, c) == "x $ineq_sym y"
    obj = constraint_object(c)
    @test obj.func == [x; y]
    @test obj.set == set
    return
end

# FIXME: should this fail?
function test_inequality_two_scalars_only_one_being_int()
    model = Model()
    @variable(model, -4 <= x <= 4)
    @variable(model, -4 <= y <= 4, Int)
    c = @constraint(model, x != y)
    ineq_sym = JuMP._math_symbol(MIME("text/plain"), :(!=))
    set = MOI.AllDifferent(2)
    @test sprint(show, c) == "x $ineq_sym y"
    obj = constraint_object(c)
    @test obj.func == [x; y]
    @test obj.set == set
    return
end

# FIXME: should this fail?
function test_inequality_two_scalars_only_one_being_bin()
    model = Model()
    @variable(model, -4 <= x <= 4)
    @variable(model, y, Bin)
    c = @constraint(model, x != y)
    ineq_sym = JuMP._math_symbol(MIME("text/plain"), :(!=))
    set = MOI.AllDifferent(2)
    @test sprint(show, c) == "x $ineq_sym y"
    obj = constraint_object(c)
    @test obj.func == [x; y]
    @test obj.set == set
    return
end

# FIXME: should this fail?
function test_inequality_two_scalars_int_vs_bin()
    model = Model()
    @variable(model, -4 <= x <= 4, Int)
    @variable(model, y, Bin)
    c = @constraint(model, x != y)
    ineq_sym = JuMP._math_symbol(MIME("text/plain"), :(!=))
    set = MOI.AllDifferent(2)
    @test sprint(show, c) == "x $ineq_sym y"
    obj = constraint_object(c)
    @test obj.func == [x; y]
    @test obj.set == set
    return
end

# FIXME: should this fail?
function test_inequality_two_scalars_real_scalars()
    model = Model()
    @variable(model, -4 <= x <= 4)
    @variable(model, -4 <= y <= 4)
    c = @constraint(model, x != y)
    ineq_sym = JuMP._math_symbol(MIME("text/plain"), :(!=))
    set = MOI.AllDifferent(2)
    @test sprint(show, c) == "x $ineq_sym y"
    obj = constraint_object(c)
    @test obj.func == [x; y]
    @test obj.set == set
    return
end

function test_inequality_two_vectors_vectorized()
    model = Model()
    @variable(model, -4 <= x[1:3] <= 4, Int)
    @variable(model, -4 <= y[1:3] <= 4, Int)
    c = @constraint(model, x .!= y)
    ineq_sym = JuMP._math_symbol(MIME("text/plain"), :(!=))
    set = MOI.AllDifferent(2)
    for (i, ci) in enumerate(c)
        @test sprint(show, ci) == "x[$i] $ineq_sym y[$i]"
        obj = constraint_object(ci)
        @test obj.func == [x[i]; y[i]]
        @test obj.set == set
    end
    return
end

function test_inequality_two_vectors_nonvectorized()
    model = Model()
    @variable(model, -4 <= x[1:3] <= 4, Int)
    @variable(model, -4 <= y[1:3] <= 4, Int)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x != y)`: Ineqality operator with " *
            "vector operands must be explicitly vectorized, " *
            "use `.!=` instead of `!=`.",
        ),
        @constraint(model, x != y)
    )
    return
end

function test_inequality_two_vectors_nonvectorized_len_mismatch()
    model = Model()
    @variable(model, -4 <= x[1:3] <= 4, Int)
    @variable(model, -4 <= y[1:2] <= 4, Int)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x != y)`: Ineqality operator with " *
            "vector operands must be explicitly vectorized, " *
            "use `.!=` instead of `!=`.",
        ),
        @constraint(model, x != y)
    )
    return
end

function test_inequality_two_vectors_vectorized_len_mismatch()
    model = Model()
    @variable(model, -4 <= x[1:3] <= 4, Int)
    @variable(model, -4 <= y[1:2] <= 4, Int)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x .!= y)`: " *
            "Operand length mismatch, 3 vs 2.",
        ),
        @constraint(model, x .!= y)
    )
    return
end

function test_inequality_non_variables()
    model = Model()
    @variable(model, -4 <= x <= 4, Int)
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x != 0)`: Unsupported form of " *
            "inequality constraint. The left- and right-hand sides must both " *
            "be decision variables.",
        ),
        @constraint(model, x != 0)
    )
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, 0 != x)`: Unsupported form of " *
            "inequality constraint. The left- and right-hand sides must both " *
            "be decision variables.",
        ),
        @constraint(model, 0 != x)
    )
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, 2x != 0)`: Unsupported form of " *
            "inequality constraint. The left- and right-hand sides must both " *
            "be decision variables.",
        ),
        @constraint(model, 2 * x != 0)
    )
    @test_throws_runtime(
        ErrorException(
            "In `@constraint(model, x != 2x)`: Unsupported form of " *
            "inequality constraint. The left- and right-hand sides must both " *
            "be decision variables.",
        ),
        @constraint(model, x != 2 * x)
    )
    return
end

end  # module
