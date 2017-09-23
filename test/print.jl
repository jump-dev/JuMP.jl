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
using Base.Test
import JuMP.REPLMode, JuMP.IJuliaMode
import JuMP.repl, JuMP.ijulia

# Helper function to test IO methods work correctly
function io_test(mode, obj, exp_str; repl=:both)
    if mode == REPLMode
        repl != :show  && @test sprint(print, obj) == exp_str
        repl != :print && @test sprint(show,  obj) == exp_str
    else
        @test sprint(show, "text/latex", obj) == string("\$\$ ",exp_str," \$\$")
    end
end


@testset "Printing" begin

    @testset "expressions" begin
        # Most of the expression logic is well covered by test/operator.jl
        # This is really just to check IJulia printing for expressions
        le, ge = repl[:leq], repl[:geq]

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

        ex = @expression(mod, (x[1]+x[2])*(y[2,2]+3.0))
        io_test(REPLMode, ex, "x[1]*y[2,2] + x[2]*y[2,2] + 3 x[1] + 3 x[2]")
        io_test(IJuliaMode, ex, "x_{1}\\times y_{2,2} + x_{2}\\times y_{2,2} + 3 x_{1} + 3 x_{2}")

        ex = @expression(mod, (y[2,2]+3.0)*(x[1]+x[2]))
        io_test(REPLMode, ex, "y[2,2]*x[1] + y[2,2]*x[2] + 3 x[1] + 3 x[2]")
        io_test(IJuliaMode, ex, "y_{2,2}\\times x_{1} + y_{2,2}\\times x_{2} + 3 x_{1} + 3 x_{2}")

        ex = @expression(mod, (x[1]+x[2])*(y[2,2]+3.0) + z^2 - 1)
        io_test(REPLMode, ex, "x[1]*y[2,2] + x[2]*y[2,2] + z$(repl[:sq]) + 3 x[1] + 3 x[2] - 1")
        io_test(IJuliaMode, ex, "x_{1}\\times y_{2,2} + x_{2}\\times y_{2,2} + z$(ijulia[:sq]) + 3 x_{1} + 3 x_{2} - 1")

        ex = @expression(mod, -z*x[1] - x[1]*z + x[1]*x[2] + 0*z^2)
        io_test(REPLMode, ex, "-2 z*x[1] + x[1]*x[2]")
        io_test(IJuliaMode, ex, "-2 z\\times x_{1} + x_{1}\\times x_{2}")
    end

    @testset "Variable" begin
        m = Model()
        @variable(m, 0 <= x <= 2)

        @test    JuMP.name(x) == "x"
        io_test(REPLMode,   x, "x")
        io_test(IJuliaMode, x, "x")

        JuMP.setname(x, "x2")
        @test    JuMP.name(x) == "x2"
        io_test(REPLMode,   x, "x2")
        io_test(IJuliaMode, x, "x2")

        JuMP.deletename(x)
        @test    JuMP.name(x) == "noname"
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

    @testset "User-created Array{Variable}" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)

        v = [x,y,x]
        A = [x y; y x]
        io_test(REPLMode,   v, "JuMP.Variable[x, y, x]")
        #io_test(IJuliaMode, v, "JuMP.Variable[x, y, x]")

        io_test(REPLMode,   A, "JuMP.Variable[x y; y x]")
        #io_test(IJuliaMode, A, "JuMP.Variable[x y; y x]")
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

end
