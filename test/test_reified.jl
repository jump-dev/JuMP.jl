#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestReified

using JuMP
using Test

include(joinpath(@__DIR__, "utilities.jl"))

function test_reified_all_different()
    model = Model()
    @variable(model, 1 <= x[1:4] <= 4, Int)
    @variable(model, z, Bin)
    c = @constraint(model, z <--> {x in MOI.AllDifferent(4)})
    in_sym = JuMP._math_symbol(MIME("text/plain"), :in)
    set = MOI.AllDifferent(4)
    @test sprint(show, c) == "z <--> {[x[1], x[2], x[3], x[4]] $in_sym $set}"
    obj = constraint_object(c)
    @test obj.func == [z; x]
    @test obj.set == MOI.Reified(set)
    return
end

function test_reified_iff()
    model = Model()
    @variable(model, 1 <= x[1:4] <= 4, Int)
    @variable(model, z, Bin)
    c = @constraint(model, z ⟺ {x in MOI.AllDifferent(4)})
    in_sym = JuMP._math_symbol(MIME("text/plain"), :in)
    set = MOI.AllDifferent(4)
    @test sprint(show, c) == "z <--> {[x[1], x[2], x[3], x[4]] $in_sym $set}"
    obj = constraint_object(c)
    @test obj.func == [z; x]
    @test obj.set == MOI.Reified(set)
    return
end

function test_reified_equal_to()
    model = Model()
    @variable(model, x)
    @variable(model, z, Bin)
    c = @constraint(model, z <--> {x == 1})
    eq_sym = JuMP._math_symbol(MIME("text/plain"), :eq)
    @test sprint(show, c) == "z <--> {x $eq_sym 1}"
    obj = constraint_object(c)
    @test obj.func == AffExpr[z, x]
    @test obj.set == MOI.Reified(MOI.EqualTo(1.0))
    return
end

function test_reified_inconsistent()
    model = Model()
    @variable(model, x)
    @variable(model, z, Bin)
    expr = :(z <--> {x .== 1})
    @test_throws_parsetime(
        ErrorException(
            "In `@constraint(model, $expr)`: " *
            "vectorized constraints cannot be used with reification.",
        ),
        @constraint(model, z <--> {x .== 1})
    )
    return
end

function test_reified_no_curly_bracket()
    model = Model()
    @variable(model, x)
    @variable(model, z, Bin)
    expr = :(z <--> x == 1)
    @test_throws_parsetime(
        ErrorException(
            "In `@constraint(model, $expr)`: " *
            "Invalid right-hand side `x == 1` of reified constraint. " *
            "Expected constraint surrounded by `{` and `}`.",
        ),
        @constraint(model, z <--> x == 1)
    )
    return
end

end  # module
