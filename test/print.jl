#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/print.jl
# Testing $fa pretty-printing-related functionality
#############################################################################

using JuMP
using LinearAlgebra, Test
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

# Used to test extensibility of JuMP printing for `JuMP.AbstractConstraint`
struct CustomConstraint{S <: JuMP.AbstractShape} <: JuMP.AbstractConstraint
    function_str::String
    in_set_str::String
    shape::S
end
function JuMP.function_string(print_mode,
                              constraint::CustomConstraint)
    return constraint.function_str
end
function JuMP.in_set_string(print_mode,
                            constraint::CustomConstraint)
    return constraint.in_set_str
end
struct CustomIndex
    value::Int
end
function JuMP.add_constraint(model::JuMP.Model, constraint::CustomConstraint,
                             name::String)
    if !haskey(model.ext, :custom)
        model.ext[:custom_constraints] = CustomConstraint[]
        model.ext[:custom_names] = String[]
    end
    constraints = model.ext[:custom_constraints]
    push!(constraints, constraint)
    push!(model.ext[:custom_names], name)
    return JuMP.ConstraintRef(model, CustomIndex(length(constraints)),
                              constraint.shape)
end
function JuMP.constraint_object(cref::JuMP.ConstraintRef{JuMP.Model,
                                                         CustomIndex})
    return cref.model.ext[:custom_constraints][cref.index.value]
end
function JuMP.name(cref::JuMP.ConstraintRef{JuMP.Model, CustomIndex})
    return cref.model.ext[:custom_names][cref.index.value]
end

@testset "Printing" begin

    @testset "expressions" begin
        # Most of the expression logic is well covered by test/operator.jl
        # This is really just to check IJulia printing for expressions
        le = JuMP._math_symbol(REPLMode, :leq)
        ge = JuMP._math_symbol(REPLMode, :geq)

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
        repl_sq = JuMP._math_symbol(REPLMode, :sq)
        io_test(REPLMode, ex, "x[1]*y[2,2] + x[2]*y[2,2] + z$repl_sq + 3 x[1] + 3 x[2] - 1")
        ijulia_sq = JuMP._math_symbol(IJuliaMode, :sq)
        io_test(IJuliaMode, ex, "x_{1}\\times y_{2,2} + x_{2}\\times y_{2,2} + z$ijulia_sq + 3 x_{1} + 3 x_{2} - 1")

        ex = @expression(mod, -z*x[1] - x[1]*z + x[1]*x[2] + 0*z^2)
        io_test(REPLMode, ex, "-2 x[1]*z + x[1]*x[2]")
        io_test(IJuliaMode, ex, "-2 x_{1}\\times z + x_{1}\\times x_{2}")

        ex = 0 * z^2 + 0 * x[1]
        io_test(REPLMode, ex, "0 z² + 0 x[1]")
        io_test(IJuliaMode, ex, "0 z$ijulia_sq + 0 x_{1}")

    end

    # See https://github.com/JuliaOpt/JuMP.jl/pull/1352
    @testset "Expression of coefficient type with unit" begin
        m = Model()
        @variable m x
        @variable m y
        u = UnitNumber(2.0)
        aff = JuMP.GenericAffExpr(zero(u), x => u, y => zero(u))
        io_test(REPLMode,   aff, "UnitNumber(2.0) x + UnitNumber(0.0) y")
        io_test(IJuliaMode, aff, "UnitNumber(2.0) x + UnitNumber(0.0) y")
        drop_zeros!(aff)
        io_test(REPLMode,   aff, "UnitNumber(2.0) x")
        io_test(IJuliaMode, aff, "UnitNumber(2.0) x")
        quad = aff * x
        io_test(REPLMode,   quad, "UnitNumber(2.0) x² + UnitNumber(0.0)")
        io_test(IJuliaMode, quad, "UnitNumber(2.0) x^2 + UnitNumber(0.0)")
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

    @testset "NLPEvaluator" begin
        model = Model()
        evaluator = JuMP.NLPEvaluator(model)
        io_test(REPLMode, evaluator, "\"A JuMP.NLPEvaluator\"")
    end

    @testset "Nonlinear constraints" begin
        le = JuMP._math_symbol(REPLMode, :leq)
        ge = JuMP._math_symbol(REPLMode, :geq)
        eq = JuMP._math_symbol(REPLMode, :eq)

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
        le = JuMP._math_symbol(REPLMode, :leq)

        model = Model()
        @variable(model, x)
        expr = @NLexpression(model, x + 1)
        @NLparameter(model, param == 1.0)

        constr = @NLconstraint(model, expr - param <= 0)
        io_test(REPLMode, constr, "(subexpression[1] - parameter[1]) - 0.0 $le 0")
        io_test(IJuliaMode, constr, "(subexpression_{1} - parameter_{1}) - 0.0 \\leq 0")
    end

    @testset "Custom constraint" begin
        model = Model()
        function test_constraint(function_str, in_set_str, name)
            constraint = CustomConstraint(function_str, in_set_str,
                                          JuMP.ScalarShape())
            cref = JuMP.add_constraint(model, constraint, name)
            @test string(cref) == "$name : $function_str $in_set_str"
        end
        test_constraint("fun", "set", "name")
        test_constraint("a", "b", "c")
    end
end

function printing_test(ModelType::Type{<:JuMP.AbstractModel})
    @testset "VariableRef" begin
        m = ModelType()
        @variable(m, 0 <= x <= 2)

        @test    JuMP.name(x) == "x"
        io_test(REPLMode,   x, "x")
        io_test(IJuliaMode, x, "x")

        JuMP.set_name(x, "x2")
        @test    JuMP.name(x) == "x2"
        io_test(REPLMode,   x, "x2")
        io_test(IJuliaMode, x, "x2")

        JuMP.set_name(x, "")
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

    @testset "base_name keyword argument" begin
        m = ModelType()
        @variable(m, x, base_name="foo")
        @variable(m, y[1:3], base_name="bar")
        num = 123
        @variable(m, z[[:red,:blue]], base_name="color_$num")
        @variable(m, v[1:2,1:2], PSD, base_name=string("i","$num",num))
        @variable(m, w[1:3,1:3], Symmetric, base_name="symm")

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

    @testset "VectorOfVariable constraints" begin
        ge = JuMP._math_symbol(REPLMode, :geq)
        in_sym = JuMP._math_symbol(REPLMode, :in)
        model = ModelType()
        @variable(model, x)
        @variable(model, y)
        zero_constr = @constraint(model, [x, y] in MOI.Zeros(2))

        io_test(REPLMode, zero_constr,
                "[x, y] $in_sym MathOptInterface.Zeros(2)")
        # TODO: Test in IJulia mode and do nice printing for Zeros().
    end

    @testset "Scalar AffExpr constraints" begin
        le = JuMP._math_symbol(REPLMode, :leq)
        ge = JuMP._math_symbol(REPLMode, :geq)
        eq = JuMP._math_symbol(REPLMode, :eq)
        in_sym = JuMP._math_symbol(REPLMode, :in)

        model = ModelType()
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
        in_sym = JuMP._math_symbol(REPLMode, :in)

        model = ModelType()
        @variable(model, x)
        @constraint(model, soc_constr, [x - 1, x + 1] in SecondOrderCone())

        io_test(REPLMode, soc_constr, "soc_constr : " *
                   "[x - 1, x + 1] $in_sym MathOptInterface.SecondOrderCone(2)")

        # TODO: Test in IJulia mode.
    end

    @testset "Scalar QuadExpr constraints" begin
        in_sym = JuMP._math_symbol(REPLMode, :in)
        le = JuMP._math_symbol(REPLMode, :leq)
        sq = JuMP._math_symbol(REPLMode, :sq)

        model = ModelType()
        @variable(model, x)
        quad_constr = @constraint(model, 2x^2 <= 1)

        io_test(REPLMode, quad_constr, "2 x$sq $le 1.0")
        # TODO: Test in IJulia mode.
    end
end

# Test printing of models of type `ModelType` for which the model is stored in
# an MOI backend
function model_printing_test(ModelType::Type{<:JuMP.AbstractModel})
    @testset "Model" begin
        repl(s) = JuMP._math_symbol(REPLMode, s)
        le, ge, eq, fa = repl(:leq), repl(:geq), repl(:eq), repl(:for_all)
        inset, dots = repl(:in), repl(:dots)
        infty, union = repl(:infty), repl(:union)
        Vert, sub2 = repl(:Vert), repl(:sub2)
        for_all = repl(:for_all)

        #------------------------------------------------------------------

        model_1 = ModelType()
        @variable(model_1, a>=1)
        @variable(model_1, b<=1)
        @variable(model_1, -1<=c<=1)
        @variable(model_1, a1>=1,Int)
        @variable(model_1, b1<=1,Int)
        @variable(model_1, -1<=c1<=1,Int)
        @variable(model_1, x, Bin)
        @variable(model_1, y)
        @variable(model_1, z, Int)
        @variable(model_1, u[1:3], Bin)
        @variable(model_1, fi == 9)
        @objective(model_1, Max, a - b + 2a1 - 10x)
        @constraint(model_1, a + b - 10c - 2x + c1 <= 1)
        @constraint(model_1, a*b <= 2)
        @constraint(model_1, [1 - a; u] in SecondOrderCone())
        @constraint(model_1, [a b; c x] in PSDCone())
        @constraint(model_1, Symmetric([a b; b x]) in PSDCone())
        @constraint(model_1, [a, b, c] in MOI.PositiveSemidefiniteConeTriangle(2))
        @constraint(model_1, [a, b, c, x] in MOI.PositiveSemidefiniteConeSquare(2))

        VariableType = typeof(a)

        io_test(REPLMode, model_1, """
    Max a - b + 2 a1 - 10 x
    Subject to
     a + b - 10 c - 2 x + c1 $le 1.0
     a*b $le 2.0
     [a  b;
      b  x] $inset PSDCone()
     [a, b, c] $inset MathOptInterface.PositiveSemidefiniteConeTriangle(2)
     [a  b;
      c  x] $inset PSDCone()
     [a, b, c, x] $inset MathOptInterface.PositiveSemidefiniteConeSquare(2)
     [-a + 1, u[1], u[2], u[3]] $inset MathOptInterface.SecondOrderCone(4)
     fi $eq 9.0
     a $ge 1.0
     c $ge -1.0
     a1 $ge 1.0
     c1 $ge -1.0
     b $le 1.0
     c $le 1.0
     b1 $le 1.0
     c1 $le 1.0
     a1 integer
     b1 integer
     c1 integer
     z integer
     x binary
     u[1] binary
     u[2] binary
     u[3] binary
    """, repl=:print)

        io_test(REPLMode, model_1, """
    A JuMP Model
    Maximization problem with:
    Variables: 13
    Objective function type: GenericAffExpr{Float64,$VariableType}
    `GenericAffExpr{Float64,$VariableType}`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
    `GenericQuadExpr{Float64,$VariableType}`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
    `Array{$VariableType,1}`-in-`MathOptInterface.PositiveSemidefiniteConeTriangle`: 2 constraints
    `Array{$VariableType,1}`-in-`MathOptInterface.PositiveSemidefiniteConeSquare`: 2 constraints
    `Array{GenericAffExpr{Float64,$VariableType},1}`-in-`MathOptInterface.SecondOrderCone`: 1 constraint
    `$VariableType`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
    `$VariableType`-in-`MathOptInterface.GreaterThan{Float64}`: 4 constraints
    `$VariableType`-in-`MathOptInterface.LessThan{Float64}`: 4 constraints
    `$VariableType`-in-`MathOptInterface.Integer`: 4 constraints
    `$VariableType`-in-`MathOptInterface.ZeroOne`: 4 constraints
    Model mode: AUTOMATIC
    CachingOptimizer state: NO_OPTIMIZER
    Solver name: No optimizer attached.
    Names registered in the model: a, a1, b, b1, c, c1, fi, u, x, y, z""", repl=:show)

        io_test(IJuliaMode, model_1, """
    \\begin{alignat*}{1}\\max\\quad & a - b + 2 a1 - 10 x\\\\
    \\text{Subject to} \\quad & a + b - 10 c - 2 x + c1 \\leq 1.0\\\\
     & a\\times b \\leq 2.0\\\\
     & \\begin{bmatrix}
    a & b\\\\
    \\cdot & x\\\\
    \\end{bmatrix} \\in PSDCone()\\\\
     & [a, b, c] \\in MathOptInterface.PositiveSemidefiniteConeTriangle(2)\\\\
     & \\begin{bmatrix}
    a & b\\\\
    c & x\\\\
    \\end{bmatrix} \\in PSDCone()\\\\
     & [a, b, c, x] \\in MathOptInterface.PositiveSemidefiniteConeSquare(2)\\\\
     & [-a + 1, u_{1}, u_{2}, u_{3}] \\in MathOptInterface.SecondOrderCone(4)\\\\
     & fi = 9.0\\\\
     & a \\geq 1.0\\\\
     & c \\geq -1.0\\\\
     & a1 \\geq 1.0\\\\
     & c1 \\geq -1.0\\\\
     & b \\leq 1.0\\\\
     & c \\leq 1.0\\\\
     & b1 \\leq 1.0\\\\
     & c1 \\leq 1.0\\\\
     & a1 integer\\\\
     & b1 integer\\\\
     & c1 integer\\\\
     & z integer\\\\
     & x binary\\\\
     & u_{1} binary\\\\
     & u_{2} binary\\\\
     & u_{3} binary\\\\
    \\end{alignat*}
    """)

        #------------------------------------------------------------------

        model_2 = ModelType()
        @variable(model_2, x, Bin)
        @variable(model_2, y, Int)
        @constraint(model_2, x*y <= 1)

        io_test(REPLMode, model_2, """
    A JuMP Model
    Feasibility problem with:
    Variables: 2
    `GenericQuadExpr{Float64,$VariableType}`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
    `$VariableType`-in-`MathOptInterface.Integer`: 1 constraint
    `$VariableType`-in-`MathOptInterface.ZeroOne`: 1 constraint
    Model mode: AUTOMATIC
    CachingOptimizer state: NO_OPTIMIZER
    Solver name: No optimizer attached.
    Names registered in the model: x, y""", repl=:show)

        model_3 = ModelType()
        @variable(model_3, x)
        @constraint(model_3, x <= 3)

        io_test(REPLMode, model_3, """
    A JuMP Model
    Feasibility problem with:
    Variable: 1
    `GenericAffExpr{Float64,$VariableType}`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
    Model mode: AUTOMATIC
    CachingOptimizer state: NO_OPTIMIZER
    Solver name: No optimizer attached.
    Names registered in the model: x""", repl=:show)
    end
end

# Test printing of models of type `ModelType` for which the model is stored in
# its JuMP form, e.g., as `AbstractVariable`s and `AbstractConstraint`s.
# This is used by `JuMPExtension` but can also be used by external packages such
# as `StructJuMP`, see https://github.com/JuliaOpt/JuMP.jl/issues/1711
function model_extension_printing_test(ModelType::Type{<:JuMP.AbstractModel})
    @testset "Model" begin
        repl(s) = JuMP._math_symbol(REPLMode, s)
        le, ge, eq, fa = repl(:leq), repl(:geq), repl(:eq), repl(:for_all)
        inset, dots = repl(:in), repl(:dots)
        infty, union = repl(:infty), repl(:union)
        Vert, sub2 = repl(:Vert), repl(:sub2)
        for_all = repl(:for_all)

        #------------------------------------------------------------------

        model_1 = ModelType()
        @variable(model_1, a>=1)
        @variable(model_1, b<=1)
        @variable(model_1, -1<=c<=1)
        @variable(model_1, a1>=1,Int)
        @variable(model_1, b1<=1,Int)
        @variable(model_1, -1<=c1<=1,Int)
        @variable(model_1, x, Bin)
        @variable(model_1, y)
        @variable(model_1, z, Int)
        @variable(model_1, u[1:3], Bin)
        @variable(model_1, fi == 9)
        @objective(model_1, Max, a - b + 2a1 - 10x)
        @constraint(model_1, a + b - 10c - 2x + c1 <= 1)
        @constraint(model_1, a*b <= 2)
        @constraint(model_1, [1 - a; u] in SecondOrderCone())

        VariableType = typeof(a)

        # TODO variable constraints
        io_test(REPLMode, model_1, """
    Max a - b + 2 a1 - 10 x
    Subject to
     a + b - 10 c - 2 x + c1 $le 1.0
     a*b $le 2.0
     [-a + 1, u[1], u[2], u[3]] $inset MathOptInterface.SecondOrderCone(4)
    """, repl=:print)

        io_test(REPLMode, model_1, """
    A JuMP Model
    Maximization problem with:
    Variables: 13
    Objective function type: GenericAffExpr{Float64,$VariableType}
    Constraints: 3
    Names registered in the model: a, a1, b, b1, c, c1, fi, u, x, y, z""", repl=:show)

        io_test(IJuliaMode, model_1, """
    \\begin{alignat*}{1}\\max\\quad & a - b + 2 a1 - 10 x\\\\
    \\text{Subject to} \\quad & a + b - 10 c - 2 x + c1 \\leq 1.0\\\\
     & a\\times b \\leq 2.0\\\\
     & [-a + 1, u_{1}, u_{2}, u_{3}] \\in MathOptInterface.SecondOrderCone(4)\\\\
    \\end{alignat*}
    """)

        #------------------------------------------------------------------

        model_2 = ModelType()
        @variable(model_2, x, Bin)
        @variable(model_2, y, Int)
        @constraint(model_2, x*y <= 1)

        io_test(REPLMode, model_2, """
    A JuMP Model
    Feasibility problem with:
    Variables: 2
    Constraint: 1
    Names registered in the model: x, y""", repl=:show)

        model_3 = ModelType()
        @variable(model_3, x)
        @constraint(model_3, x <= 3)

        io_test(REPLMode, model_3, """
    A JuMP Model
    Feasibility problem with:
    Variable: 1
    Constraint: 1
    Names registered in the model: x""", repl=:show)
    end
end

@testset "Printing for JuMP.Model" begin
    printing_test(Model)
    model_printing_test(Model)
    @testset "Model with nonlinear terms" begin
        eq = JuMP._math_symbol(REPLMode, :eq)
        model = Model()
        @variable(model, x)
        @NLobjective(model, Max, sin(x))
        @NLconstraint(model, cos(x) == 0)

        io_test(REPLMode, model, """
    A JuMP Model
    Maximization problem with:
    Variable: 1
    Objective function type: Nonlinear
    Nonlinear: 1 constraint
    Model mode: AUTOMATIC
    CachingOptimizer state: NO_OPTIMIZER
    Solver name: No optimizer attached.
    Names registered in the model: x""", repl=:show)

        io_test(REPLMode, model, """
    Max sin(x)
    Subject to
     cos(x) - 0.0 $eq 0
    """, repl=:print)

        io_test(IJuliaMode, model, """
    \\begin{alignat*}{1}\\max\\quad & sin(x)\\\\
    \\text{Subject to} \\quad & cos(x) - 0.0 = 0\\\\
    \\end{alignat*}
    """)
    end
    @testset "SingleVariable constraints" begin
        ge = JuMP._math_symbol(REPLMode, :geq)
        in_sym = JuMP._math_symbol(REPLMode, :in)
        model = Model()
        @variable(model, x >= 10)
        zero_one = @constraint(model, x in MOI.ZeroOne())

        io_test(REPLMode, JuMP.LowerBoundRef(x), "x $ge 10.0")
        io_test(REPLMode, zero_one, "x binary")
        # TODO: Test in IJulia mode
    end
end

@testset "Printing for JuMPExtension.MyModel" begin
    printing_test(JuMPExtension.MyModel)
    model_extension_printing_test(JuMPExtension.MyModel)
end
