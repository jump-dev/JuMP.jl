#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/print.jl
# Testing $fa pretty-printing-related functionality
#############################################################################
using JuMP
using Compat
using Compat.Test
import JuMP.REPLMode, JuMP.IJuliaMode

# Helper function to test IO methods work correctly
function io_test(mode, obj, exp_str; repl=:both)
    if mode == REPLMode
        repl != :show  && @test sprint(print, obj) == exp_str
        repl != :print && @test sprint(show,  obj) == exp_str
    else
        @test sprint(show, "text/latex", obj) == string("\$\$ ",exp_str," \$\$")
    end
end

# Used to test that JuMP printing works correctly for types for which
# oneunit is not convertible to Float64
struct UnitNumber <: Number
    α::Float64
end
Base.zero(::Union{UnitNumber, Type{UnitNumber}}) = UnitNumber(0.0)
Base.oneunit(::Union{UnitNumber, Type{UnitNumber}}) = UnitNumber(1.0)
Base.:(+)(u::UnitNumber, v::UnitNumber) = UnitNumber(u.α + v.α)
Base.:(-)(u::UnitNumber, v::UnitNumber) = UnitNumber(u.α - v.α)
Base.:(*)(α::Float64, u::UnitNumber) = UnitNumber(α * u.α)
Base.abs(u::UnitNumber) = UnitNumber(abs(u.α))
Base.isless(u::UnitNumber, v::UnitNumber) = isless(u.α, v.α)

@testset "Printing" begin

    @testset "expressions" begin
        # Most of the expression logic is well covered by test/operator.jl
        # This is really just to check IJulia printing for expressions
        le = JuMP.math_symbol(REPLMode, :leq)
        ge = JuMP.math_symbol(REPLMode, :geq)

        #------------------------------------------------------------------
        mod = Model()
        @variable(mod, x[1:5])
        @variable(mod, y[i=2:4,j=i:5])
        @variable(mod, z)

        ex = @expression(mod, x[1] + 2*y[2,3])
        io_test(REPLMode, ex, "x[1] + 2 y[2,3]")
        io_test(IJuliaMode, ex, "x_{1} + 2 y_{2,3}")

        ex = @expression(mod, x[1] + 2*y[2,3] + x[1])
        io_test(REPLMode, ex, "2 x[1] + 2 y[2,3]")
        io_test(IJuliaMode, ex, "2 x_{1} + 2 y_{2,3}")

        # TODO: These tests shouldn't depend on order of the two variables in
        # quadratic terms, i.e., x*y vs y*x.
        ex = @expression(mod, (x[1]+x[2])*(y[2,2]+3.0))
        io_test(REPLMode, ex, "x[1]*y[2,2] + x[2]*y[2,2] + 3 x[1] + 3 x[2]")
        io_test(IJuliaMode, ex, "x_{1}\\times y_{2,2} + x_{2}\\times y_{2,2} + 3 x_{1} + 3 x_{2}")

        ex = @expression(mod, (y[2,2]+3.0)*(x[1]+x[2]))
        io_test(REPLMode, ex, "y[2,2]*x[1] + y[2,2]*x[2] + 3 x[1] + 3 x[2]")
        io_test(IJuliaMode, ex, "y_{2,2}\\times x_{1} + y_{2,2}\\times x_{2} + 3 x_{1} + 3 x_{2}")

        ex = @expression(mod, (x[1]+x[2])*(y[2,2]+3.0) + z^2 - 1)
        repl_sq = JuMP.math_symbol(REPLMode, :sq)
        io_test(REPLMode, ex, "x[1]*y[2,2] + x[2]*y[2,2] + z$repl_sq + 3 x[1] + 3 x[2] - 1")
        ijulia_sq = JuMP.math_symbol(IJuliaMode, :sq)
        io_test(IJuliaMode, ex, "x_{1}\\times y_{2,2} + x_{2}\\times y_{2,2} + z$ijulia_sq + 3 x_{1} + 3 x_{2} - 1")

        ex = @expression(mod, -z*x[1] - x[1]*z + x[1]*x[2] + 0*z^2)
        io_test(REPLMode, ex, "-2 x[1]*z + x[1]*x[2]")
        io_test(IJuliaMode, ex, "-2 x_{1}\\times z + x_{1}\\times x_{2}")
    end

    @testset "VariableRef" begin
        m = Model()
        @variable(m, 0 <= x <= 2)

        @test    JuMP.name(x) == "x"
        io_test(REPLMode,   x, "x")
        io_test(IJuliaMode, x, "x")

        JuMP.setname(x, "x2")
        @test    JuMP.name(x) == "x2"
        io_test(REPLMode,   x, "x2")
        io_test(IJuliaMode, x, "x2")

        JuMP.setname(x, "")
        @test    JuMP.name(x) == ""
        io_test(REPLMode,   x, "noname")
        io_test(IJuliaMode, x, "noname")

        @variable(m, z[1:2,3:5])
        @test       JuMP.name(z[1,3]) == "z[1,3]"
        io_test(REPLMode,   z[1,3],    "z[1,3]")
        io_test(IJuliaMode, z[1,3],    "z_{1,3}")
        @test       JuMP.name(z[2,4]) == "z[2,4]"
        io_test(REPLMode,   z[2,4],    "z[2,4]")
        io_test(IJuliaMode, z[2,4],    "z_{2,4}")
        @test       JuMP.name(z[2,5]) == "z[2,5]"
        io_test(REPLMode,   z[2,5],    "z[2,5]")
        io_test(IJuliaMode, z[2,5],    "z_{2,5}")

        @variable(m, w[3:9,["red","blue","green"]])
        @test    JuMP.name(w[7,"green"]) == "w[7,green]"
        io_test(REPLMode,   w[7,"green"], "w[7,green]")
        io_test(IJuliaMode, w[7,"green"], "w_{7,green}")

        rng = 2:5
        @variable(m, v[rng,rng,rng,rng,rng,rng,rng])
        a_v = v[4,5,2,3,2,2,4]
        @test    JuMP.name(a_v) == "v[4,5,2,3,2,2,4]"
        io_test(REPLMode,   a_v, "v[4,5,2,3,2,2,4]")
        io_test(IJuliaMode, a_v, "v_{4,5,2,3,2,2,4}")
    end

    @testset "User-created Array{VariableRef}" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)

        v = [x,y,x]
        A = [x y; y x]
        io_test(REPLMode,   v, "JuMP.VariableRef[x, y, x]")
        #io_test(IJuliaMode, v, "JuMP.VariableRef[x, y, x]")

        io_test(REPLMode,   A, "JuMP.VariableRef[x y; y x]")
        #io_test(IJuliaMode, A, "JuMP.VariableRef[x y; y x]")
    end

    @testset "basename keyword argument" begin
        m = Model()
        @variable(m, x, basename="foo")
        @variable(m, y[1:3], basename="bar")
        num = 123
        @variable(m, z[[:red,:blue]], basename="color_$num")
        @variable(m, v[1:2,1:2], PSD, basename=string("i","$num",num))
        @variable(m, w[1:3,1:3], Symmetric, basename="symm")

        io_test(REPLMode,   x, "foo")
        io_test(IJuliaMode, x, "foo")
        io_test(REPLMode,   y[2], "bar[2]")
        io_test(IJuliaMode, y[2], "bar_{2}")
        io_test(REPLMode,   z[:red], "color_123[red]")
        io_test(IJuliaMode, z[:red], "color_123_{red}")
        io_test(REPLMode,   v[2,1], "i123123[1,2]")
        io_test(IJuliaMode, v[2,1], "i123123_{1,2}")
        io_test(REPLMode,   w[1,3], "symm[1,3]")
        io_test(IJuliaMode, w[1,3], "symm_{1,3}")
    end

    # See https://github.com/JuliaOpt/JuMP.jl/pull/1352
    @testset "Expression of coefficient type with unit" begin
        m = Model()
        @variable m x
        @variable m y
        u = UnitNumber(2.0)
        aff = JuMP.GenericAffExpr(zero(u), x => u, y => zero(u))
        io_test(REPLMode,   aff, "UnitNumber(2.0) x")
        io_test(IJuliaMode, aff, "UnitNumber(2.0) x")
        quad = aff * x
        io_test(REPLMode,   quad, "UnitNumber(2.0) x² + UnitNumber(0.0)")
        io_test(IJuliaMode, quad, "UnitNumber(2.0) x^2 + UnitNumber(0.0)")
    end

    @testset "SingleVariable constraints" begin
        ge = JuMP.math_symbol(REPLMode, :geq)
        in_sym = JuMP.math_symbol(REPLMode, :in)
        model = Model()
        @variable(model, x >= 10)
        zero_one = @constraint(model, x in MathOptInterface.ZeroOne())

        io_test(REPLMode, JuMP.LowerBoundRef(x), "x $ge 10.0")
        io_test(REPLMode, zero_one, "x $in_sym MathOptInterface.ZeroOne()")
        # TODO: Test in IJulia mode and do nice printing for {0, 1}.
    end

    @testset "VectorOfVariable constraints" begin
        ge = JuMP.math_symbol(REPLMode, :geq)
        in_sym = JuMP.math_symbol(REPLMode, :in)
        model = Model()
        @variable(model, x)
        @variable(model, y)
        zero_constr = @constraint(model, [x, y] in MathOptInterface.Zeros(2))

        io_test(REPLMode, zero_constr,
                "[x, y] $in_sym MathOptInterface.Zeros(2)")
        # TODO: Test in IJulia mode and do nice printing for Zeros().
    end

    @testset "Scalar AffExpr constraints" begin
        le = JuMP.math_symbol(REPLMode, :leq)
        ge = JuMP.math_symbol(REPLMode, :geq)
        eq = JuMP.math_symbol(REPLMode, :eq)
        in_sym = JuMP.math_symbol(REPLMode, :in)

        model = Model()
        @variable(model, x)
        @constraint(model, linear_le, x + 0 <= 1)
        @constraint(model, linear_ge, x + 0 >= 1)
        @constraint(model, linear_eq, x + 0 == 1)
        @constraint(model, linear_range, -1 <= x + 0 <= 1)
        linear_noname = @constraint(model, x + 0 <= 1)

        io_test(REPLMode, linear_le, "linear_le : x $le 1.0")
        io_test(REPLMode, linear_eq, "linear_eq : x $eq 1.0")
        io_test(REPLMode, linear_range, "linear_range : x $in_sym [-1.0, 1.0]")
        io_test(REPLMode, linear_noname, "x $le 1.0")

        # io_test doesn't work here because constraints print with a mix of math
        # and non-math.
        @test sprint(show, "text/latex", linear_le) ==
              "linear_le : \$ x \\leq 1.0 \$"
        @test sprint(show, "text/latex", linear_ge) ==
              "linear_ge : \$ x \\geq 1.0 \$"
        @test sprint(show, "text/latex", linear_eq) ==
              "linear_eq : \$ x = 1.0 \$"
        @test sprint(show, "text/latex", linear_range) ==
              "linear_range : \$ x \\in \\[-1.0, 1.0\\] \$"
        @test sprint(show, "text/latex", linear_noname) == "\$ x \\leq 1.0 \$"
    end

    @testset "Vector AffExpr constraints" begin
        in_sym = JuMP.math_symbol(REPLMode, :in)

        model = Model()
        @variable(model, x)
        @constraint(model, soc_constr, [x - 1, x + 1] in SecondOrderCone())

        io_test(REPLMode, soc_constr, "soc_constr : " *
                   "[x - 1, x + 1] $in_sym MathOptInterface.SecondOrderCone(2)")

        # TODO: Test in IJulia mode.
    end

    @testset "Scalar QuadExpr constraints" begin
        in_sym = JuMP.math_symbol(REPLMode, :in)
        le = JuMP.math_symbol(REPLMode, :leq)
        sq = JuMP.math_symbol(REPLMode, :sq)

        model = Model()
        @variable(model, x)
        quad_constr = @constraint(model, 2x^2 <= 1)

        io_test(REPLMode, quad_constr, "2 x$sq $le 1.0")
        # TODO: Test in IJulia mode.
    end

    @testset "Nonlinear expressions" begin
        model = Model()
        @variable(model, x)
        expr = @NLexpression(model, x + 1)
        io_test(REPLMode, expr, "\"Reference to nonlinear expression #1\"")
    end

    @testset "Nonlinear parameters" begin
        model = Model()
        @NLparameter(model, param == 1.0)
        io_test(REPLMode, param, "\"Reference to nonlinear parameter #1\"")
    end

    @testset "Nonlinear constraints" begin
        le = JuMP.math_symbol(REPLMode, :leq)
        ge = JuMP.math_symbol(REPLMode, :geq)
        eq = JuMP.math_symbol(REPLMode, :eq)

        model = Model()
        @variable(model, x)
        constr_le = @NLconstraint(model, sin(x) <= 1)
        constr_ge = @NLconstraint(model, sin(x) >= 1)
        constr_eq = @NLconstraint(model, sin(x) == 1)
        constr_range = @NLconstraint(model, 0 <= sin(x) <= 1)

        io_test(REPLMode, constr_le, "sin(x) - 1.0 $le 0")
        io_test(REPLMode, constr_ge, "sin(x) - 1.0 $ge 0")
        io_test(REPLMode, constr_eq, "sin(x) - 1.0 $eq 0")
        # Note: This is inconsistent with the "x in [-1, 1]" printing for
        # regular constraints.
        io_test(REPLMode, constr_range, "0 $le sin(x) $le 1")

        io_test(IJuliaMode, constr_le, "sin(x) - 1.0 \\leq 0")
        io_test(IJuliaMode, constr_ge, "sin(x) - 1.0 \\geq 0")
        io_test(IJuliaMode, constr_eq, "sin(x) - 1.0 = 0")
        io_test(IJuliaMode, constr_range, "0 \\leq sin(x) \\leq 1")
    end

    @testset "Nonlinear constraints with embedded parameters/expressions" begin
        le = JuMP.math_symbol(REPLMode, :leq)

        model = Model()
        @variable(model, x)
        expr = @NLexpression(model, x + 1)
        @NLparameter(model, param == 1.0)

        constr = @NLconstraint(model, expr - param <= 0)
        io_test(REPLMode, constr, "(subexpression[1] - parameter[1]) - 0.0 $le 0")
        io_test(IJuliaMode, constr, "(subexpression_{1} - parameter_{1}) - 0.0 \\leq 0")
    end
end
