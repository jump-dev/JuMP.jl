#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestJuMPExtension

import JuMP
import Test

include(joinpath(@__DIR__, "JuMPExtension.jl"))

# Test printing of models of type `ModelType` for which the model is stored in
# its JuMP form, for example, as `AbstractVariable`s and `AbstractConstraint`s.
# This is used by `JuMPExtension` but can also be used by external packages such
# as `StructJuMP`, see https://github.com/jump-dev/JuMP.jl/issues/1711

function test_model_extension_printing()
    repl(s) = JuMP._math_symbol(MIME("text/plain"), s)

    model_1 = JuMPExtension.MyModel()
    JuMP.@variable(model_1, a >= 1)
    JuMP.@variable(model_1, b <= 1)
    JuMP.@variable(model_1, -1 <= c <= 1)
    JuMP.@variable(model_1, a1 >= 1, Int)
    JuMP.@variable(model_1, b1 <= 1, Int)
    JuMP.@variable(model_1, -1 <= c1 <= 1, Int)
    JuMP.@variable(model_1, x, Bin)
    JuMP.@variable(model_1, y)
    JuMP.@variable(model_1, z, Int)
    JuMP.@variable(model_1, u[1:3], Bin)
    JuMP.@variable(model_1, fi == 9)
    JuMP.@objective(model_1, Max, a - b + 2a1 - 10x)
    JuMP.@constraint(model_1, a + b - 10c - 2x + c1 <= 1)
    JuMP.@constraint(model_1, a * b <= 2)
    JuMP.@constraint(model_1, [1 - a; u] in JuMP.SecondOrderCone())

    model_2 = JuMPExtension.MyModel()
    JuMP.@variable(model_2, x, Bin)
    JuMP.@variable(model_2, y, Int)
    JuMP.@constraint(model_2, x * y <= 1)

    model_3 = JuMPExtension.MyModel()
    JuMP.@variable(model_3, x)
    JuMP.@constraint(model_3, x <= 3)

    Test.@test JuMP.model_string(MIME("text/plain"), model_1) == """
Max a - b + 2 a1 - 10 x
Subject to
 a + b - 10 c - 2 x + c1 $(repl(:leq)) 1
 a*b $(repl(:leq)) 2
 [-a + 1, u[1], u[2], u[3]] $(repl(:in)) MathOptInterface.SecondOrderCone(4)
"""

    Test.@test JuMP.model_string(MIME("text/latex"), model_1) == """
\$\$ \\begin{aligned}
\\max\\quad & a - b + 2 a1 - 10 x\\\\
\\text{Subject to} \\quad & a + b - 10 c - 2 x + c1 \\leq 1\\\\
 & a\\times b \\leq 2\\\\
 & [-a + 1, u_{1}, u_{2}, u_{3}] \\in \\text{MathOptInterface.SecondOrderCone(4)}\\\\
\\end{aligned} \$\$"""

    Test.@test sprint(show, model_1) == """
An Abstract JuMP Model
Maximization problem with:
Variables: 13
Objective function type: $(JuMP.GenericAffExpr{Float64,JuMPExtension.MyVariableRef})
Constraints: 3
Names registered in the model: a, a1, b, b1, c, c1, fi, u, x, y, z"""

    Test.@test sprint(show, model_2) == """
An Abstract JuMP Model
Feasibility problem with:
Variables: 2
Constraint: 1
Names registered in the model: x, y"""

    Test.@test sprint(show, model_3) == """
An Abstract JuMP Model
Feasibility problem with:
Variable: 1
Constraint: 1
Names registered in the model: x"""
    return
end

end
