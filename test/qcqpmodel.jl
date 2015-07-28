#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/quadmodel.jl
# Testing quadratic models (objective and constraints)
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using JuMP, FactCheck

facts("[qcqpmodel] Test quad objective (discrete)") do
for solver in quad_mip_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @defVar(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    @setObjective(modQ, Min, 10*x[1]*x[1] + 3*x[1]*x[2] + 5*x[2]*x[2] + 9*x[3]*x[3])
    @addConstraint(modQ, x[2] <= 1.7*x[3])
    @addConstraint(modQ, x[2] >= 0.5*x[1])
    @fact solve(modQ) => :Optimal
    @fact modQ.objVal => roughly( 247.0, 1e-5)
    @fact getValue(x)[:] => roughly([2.0, 3.0, 4.0], 1e-6)

    # vectorized version
    modV = Model(solver=solver)
    @defVar(modV, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    obj = x'*[10 1.5 0; 1.5 5 0; 0 0 9]*x
    @setObjective(modV, Min, obj[1])
    A = [  0  1 -1.7
         0.5 -1    0]
    @addConstraint(modV, A*x .<= zeros(2))
    @fact solve(modV) => :Optimal
    @fact modV.objVal => roughly( 247.0, 1e-5)
    @fact getValue(x)[:] => roughly([2.0, 3.0, 4.0], 1e-6)

    modQ = Model(solver=solver)
    @defVar(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    @setObjective(modQ, :Max, -10*x[1]*x[1] - 3*x[1]*x[2] - 5*x[2]*x[2] - 9*x[3]*x[3])
    @addConstraint(modQ, x[2] <= 1.7*x[3])
    @addConstraint(modQ, x[2] >= 0.5*x[1])
    @fact solve(modQ) => :Optimal
    @fact modQ.objVal => roughly(-247.0, 1e-5)
    @fact getValue(x)[:] => roughly([2.0, 3.0, 4.0], 1e-6)

    # vectorized version
    modV = Model(solver=solver)
    @defVar(modV, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    Q = [10 3 0; 0 5 0; 0 0 9]
    obj = x'Q*x
    @setObjective(modV, Min, obj[1])
    A = [   0 -1 1.7
         -0.5  1   0]
    @addConstraint(modV, A*x .>= zeros(2))
    @fact solve(modV) => :Optimal
    @fact modV.objVal => roughly( 247.0, 1e-5)
    @fact getValue(x)[:] => roughly([2.0, 3.0, 4.0], 1e-6)
end; end; end

facts("[qcqpmodel] Test quad objective (continuous)") do
for solver in quad_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @defVar(modQ, 0.5 <= x <= 2 )
    @defVar(modQ, 0 <= y <= 30 )
    @setObjective(modQ, :Min, (x+y)*(x+y) )
    @addConstraint(modQ, x + y >= 1 )
    @fact solve(modQ) => :Optimal
    @fact modQ.objVal => roughly(1.0, 1e-6)
    @fact (getValue(x) + getValue(y)) => roughly(1.0, 1e-6)

    # Vectorized version
    modV = Model(solver=solver)
    @defVar(modV, 0.5 <= x <= 2 )
    @defVar(modV, 0 <= y <= 30 )
    obj = [x y]*ones(2,2)*[x,y]
    @setObjective(modV, Min, obj[1])
    @addConstraint(modV, ones(1,2)*[x,y] .>= 1)
    @fact solve(modV) => :Optimal
    @fact modV.objVal => roughly(1.0, 1e-6)
    @fact (getValue(x) + getValue(y)) => roughly(1.0, 1e-6)
end; end; end



facts("[qcqpmodel] Test quad constraints (continuous)") do
for solver in quad_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @defVar(modQ, -2 <= x <= 2 )
    @defVar(modQ, -2 <= y <= 2 )
    @setObjective(modQ, Min, x - y )
    @addConstraint(modQ, x + x*x + x*y + y*y <= 1 )
    @fact MathProgBase.numquadconstr(modQ) => 1
    @fact MathProgBase.numlinconstr(modQ) => 0
    @fact MathProgBase.numconstr(modQ) => 1

    @fact solve(modQ) => :Optimal
    @fact modQ.objVal => roughly(-1-4/sqrt(3), 1e-6)
    @fact (getValue(x) + getValue(y)) => roughly(-1/3, 1e-3)

    # Vectorized version
    modV = Model(solver=solver)
    @defVar(modV, -2 <= x <= 2 )
    @defVar(modV, -2 <= y <= 2 )
    obj = [1 -1]*[x,y]
    @setObjective(modV, Min, obj[1])
    A = [1 0.5; 0.5 1]
    @addConstraint(modV, [x y]*A*[x,y] + [x] .<= [1])
    @fact MathProgBase.numquadconstr(modV) => 1
    @fact MathProgBase.numlinconstr(modV) => 0
    @fact MathProgBase.numconstr(modV) => 1

    @fact solve(modV) => :Optimal
    @fact modV.objVal => roughly(-1-4/sqrt(3), 1e-6)
    @fact (getValue(x) + getValue(y)) => roughly(-1/3, 1e-3)

end; end; end

facts("[qcqpmodel] Test SOC constraints (continuous)") do
for solver in soc_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @defVar(modQ, x)
    @defVar(modQ, y)
    @defVar(modQ, t >= 0)
    @setObjective(modQ, Min, t)
    @addConstraint(modQ, x+y >= 1)
    @addConstraint(modQ, x^2 + y^2 <= t^2)

    @fact solve(modQ) => :Optimal
    @fact modQ.objVal => roughly(sqrt(1/2), 1e-6)
    @fact norm([getValue(x), getValue(y), getValue(t)] - [0.5,0.5,sqrt(1/2)]) => roughly(0.0,1e-3)

    # Vectorized version
    modV = Model(solver=solver)
    @defVar(modV, x)
    @defVar(modV, y)
    @defVar(modV, t >= 0)
    @setObjective(modV, Min, t)
    @addConstraint(modV, [1 1 0]*[x,y,t] .>= 1)
    Q = [1 0  0
         0 1  0
         0 0 -1]
    @addConstraint(modV, [x y t]*Q*[x,y,t] .<= 0)

    @fact solve(modV) => :Optimal
    @fact modV.objVal => roughly(sqrt(1/2), 1e-6)
    @fact norm([getValue(x), getValue(y), getValue(t)] - [0.5,0.5,sqrt(1/2)]) => roughly(0.0,1e-3)

    # SOC version 1
    modN = Model(solver=solver)
    @defVar(modN, x)
    @defVar(modN, y)
    @defVar(modN, t >= 0)
    @setObjective(modN, Min, t)
    @addConstraint(modN, x + y >= 1)
    @addConstraint(modN, norm([x,y]) <= t)

    @fact solve(modN) => :Optimal
    @fact modN.objVal => roughly(sqrt(1/2), 1e-6)
    @fact norm([getValue(x), getValue(y), getValue(t)] - [0.5,0.5,sqrt(1/2)]) => roughly(0.0,1e-3)

    # SOC version 2
    modN = Model(solver=solver)
    @defVar(modN, x)
    @defVar(modN, y)
    @defVar(modN, t >= 0)
    @setObjective(modN, Min, t)
    @addConstraint(modN, x + y >= 1)
    tmp = [x,y]
    @addConstraint(modN, norm2{tmp[i], i=1:2} <= t)

    @fact solve(modN) => :Optimal
    @fact modN.objVal => roughly(sqrt(1/2), 1e-6)
    @fact norm([getValue(x), getValue(y), getValue(t)] - [0.5,0.5,sqrt(1/2)]) => roughly(0.0,1e-3)
end; end; end

facts("[qcqpmodel] Test SOC duals") do
for solver in soc_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @defVar(modQ, x >= 0)
    @defVar(modQ, y)
    @defVar(modQ, z)
    @setObjective(modQ, Min, -y-z)
    @addConstraint(modQ, eq, x <= 1)
    @addConstraint(modQ, y^2 + z^2 <= x^2)

    @fact solve(modQ) => :Optimal
    @fact modQ.objVal => roughly(-sqrt(2), 1e-6)
    @fact getValue(y) => roughly(1/sqrt(2), 1e-6)
    @fact getValue(z) => roughly(1/sqrt(2), 1e-6)
    @fact getDual(eq) => roughly(-sqrt(2), 1e-6)

    @setObjective(modQ, Max, y+z)
    @fact solve(modQ) => :Optimal
    @fact modQ.objVal => roughly(sqrt(2), 1e-6)
    @fact getValue(y) => roughly(1/sqrt(2), 1e-6)
    @fact getValue(z) => roughly(1/sqrt(2), 1e-6)
    @fact getDual(eq) => roughly(sqrt(2), 1e-6)

    # # SOC syntax version
    # modN = Model(solver=solver)
    # @defVar(modN, x >= 0)
    # @defVar(modN, y)
    # @defVar(modN, z)
    # @setObjective(modN, Min, -y-z)
    # @addConstraint(modN, eq, x <= 1)
    # # @addConstraint(modN, y^2 + z^2 <= x^2)
    # @addConstraint(modN, (1/4)*norm([4.0y,4.0z]) <= x)
    #
    # @fact solve(modN) => :Optimal
    # @fact modN.objVal => roughly(-sqrt(2), 1e-6)
    # @fact getValue(y) => roughly(1/sqrt(2), 1e-6)
    # @fact getValue(z) => roughly(1/sqrt(2), 1e-6)
    # @fact getDual(eq) => roughly(-sqrt(2), 1e-6)
    #
    # @setObjective(modN, Max, y+z)
    # @fact solve(modN) => :Optimal
    # @fact modN.objVal => roughly(sqrt(2), 1e-6)
    # @fact getValue(y) => roughly(1/sqrt(2), 1e-6)
    # @fact getValue(z) => roughly(1/sqrt(2), 1e-6)
    # @fact getDual(eq) => roughly(sqrt(2), 1e-6)

end; end; end

facts("[qcqpmodel] Test quad constraints (discrete)") do
for solver in quad_mip_solvers
context("With solver $(typeof(solver))") do
    modQ = Model(solver=solver)
    @defVar(modQ, -2 <= x <= 2, Int )
    @defVar(modQ, -2 <= y <= 2, Int )
    @setObjective(modQ, Min, x - y )
    @addConstraint(modQ, x + x*x + x*y + y*y <= 1 )

    @fact solve(modQ) => :Optimal
    @fact modQ.objVal => roughly(-3, 1e-6)
    @fact (getValue(x) + getValue(y)) => roughly(-1, 1e-6)

end; end; end

facts("[qcqpmodel] Test simple normed problem") do
for solver in soc_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver);
    @defVar(m, x[1:3]);
    @addConstraint(m, 2norm2{x[i]-1, i=1:3} <= 2)
    @setObjective(m, Max, x[1]+x[2])

    @fact solve(m) => :Optimal
    @fact getObjectiveValue(m) => roughly(2+sqrt(2), 1e-5)
    @fact norm(getValue(x)-[1+sqrt(1/2),1+sqrt(1/2),1]) => roughly(0, 1e-6)
end; end; end

facts("[qcqpmodel] Test quad problem modification") do
for solver in quad_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @defVar(modQ, x >= 0)
    @addConstraint(modQ, x*x <= 1)
    @setObjective(modQ, Max, x)
    @fact solve(modQ) => :Optimal
    @fact getObjectiveValue(modQ) => roughly(1.0, 1e-6)

    @addConstraint(modQ, 2x*x <= 1)
    @fact modQ.internalModelLoaded => true
    @fact solve(modQ) => :Optimal
    @fact getObjectiveValue(modQ) => roughly(sqrt(0.5), 1e-6)

    modQ = Model(solver=solver)
    @defVar(modQ,   0 <= x <= 1)
    @defVar(modQ, 1/2 <= y <= 1)
    setObjective(modQ, :Min, x*x - y)
    @fact solve(modQ) => :Optimal
    @fact getObjectiveValue(modQ) => roughly(-1.0, 1e-6)

    setObjective(modQ, :Min, y*y - x)
    @fact modQ.internalModelLoaded => true
    @fact solve(modQ) => :Optimal
    @fact getObjectiveValue(modQ) => roughly(-0.75, 1e-6)
end; end; end
