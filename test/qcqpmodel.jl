#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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

# If solvers not loaded, load them (i.e running just these tests)
!isdefined(:lp_solvers) && include("solvers.jl")

facts("[qcqpmodel] Test quad objective (discrete)") do
for solver in quad_mip_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @variable(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    @objective(modQ, Min, 10*x[1]*x[1] + 3*x[1]*x[2] + 5*x[2]*x[2] + 9*x[3]*x[3])
    @constraint(modQ, x[2] <= 1.7*x[3])
    @constraint(modQ, x[2] >= 0.5*x[1])
    @fact solve(modQ) --> :Optimal
    @fact modQ.objVal --> roughly( 247.0, 1e-5)
    @fact getvalue(x)[:] --> roughly([2.0, 3.0, 4.0], 1e-6)

    # vectorized version
    modV = Model(solver=solver)
    @variable(modV, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    obj = x'*[10 1.5 0; 1.5 5 0; 0 0 9]*x
    @objective(modV, Min, obj[1])
    A = [  0  1 -1.7
         0.5 -1    0]
    @constraint(modV, A*x .<= zeros(2))
    @fact solve(modV) --> :Optimal
    @fact modV.objVal --> roughly( 247.0, 1e-5)
    @fact getvalue(x)[:] --> roughly([2.0, 3.0, 4.0], 1e-6)

    modQ = Model(solver=solver)
    @variable(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    @objective(modQ, :Max, -10*x[1]*x[1] - 3*x[1]*x[2] - 5*x[2]*x[2] - 9*x[3]*x[3])
    @constraint(modQ, x[2] <= 1.7*x[3])
    @constraint(modQ, x[2] >= 0.5*x[1])
    @fact solve(modQ) --> :Optimal
    @fact modQ.objVal --> roughly(-247.0, 1e-5)
    @fact getvalue(x)[:] --> roughly([2.0, 3.0, 4.0], 1e-6)

    # sparse vectorized version
    modV = Model(solver=solver)
    @variable(modV, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    obj = x'*sparse([10 1.5 0; 1.5 5 0; 0 0 9])*x
    @objective(modV, Min, obj[1])
    A = sparse([  0  1 -1.7
                0.5 -1    0])
    @constraint(modV, A*x .<= zeros(2))
    @fact solve(modV) --> :Optimal
    @fact modV.objVal --> roughly( 247.0, 1e-5)
    @fact getvalue(x)[:] --> roughly([2.0, 3.0, 4.0], 1e-6)

    modQ = Model(solver=solver)
    @variable(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    @objective(modQ, :Max, -10*x[1]*x[1] - 3*x[1]*x[2] - 5*x[2]*x[2] - 9*x[3]*x[3])
    @constraint(modQ, x[2] <= 1.7*x[3])
    @constraint(modQ, x[2] >= 0.5*x[1])
    @fact solve(modQ) --> :Optimal
    @fact modQ.objVal --> roughly(-247.0, 1e-5)
    @fact getvalue(x)[:] --> roughly([2.0, 3.0, 4.0], 1e-6)

    # vectorized version
    modV = Model(solver=solver)
    @variable(modV, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
    Q = [10 3 0; 0 5 0; 0 0 9]
    obj = x'Q*x
    @objective(modV, Min, obj[1])
    A = [   0 -1 1.7
         -0.5  1   0]
    @constraint(modV, A*x .>= zeros(2))
    @fact solve(modV) --> :Optimal
    @fact modV.objVal --> roughly( 247.0, 1e-5)
    @fact getvalue(x)[:] --> roughly([2.0, 3.0, 4.0], 1e-6)
end; end; end

facts("[qcqpmodel] Test quad objective (continuous)") do
for solver in quad_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @variable(modQ, 0.5 <= x <= 2 )
    @variable(modQ, 0 <= y <= 30 )
    @objective(modQ, :Min, (x+y)*(x+y) )
    @constraint(modQ, x + y >= 1 )
    @fact solve(modQ) --> :Optimal
    @fact modQ.objVal --> roughly(1.0, 1e-6)
    @fact (getvalue(x) + getvalue(y)) --> roughly(1.0, 1e-6)

    # Vectorized version
    modV = Model(solver=solver)
    @variable(modV, 0.5 <= x <= 2 )
    @variable(modV, 0 <= y <= 30 )
    obj = [x y]*ones(2,2)*[x,y]
    @objective(modV, Min, obj[1])
    @constraint(modV, ones(1,2)*[x,y] .>= 1)
    @fact solve(modV) --> :Optimal
    @fact modV.objVal --> roughly(1.0, 1e-6)
    @fact (getvalue(x) + getvalue(y)) --> roughly(1.0, 1e-6)

    # Sparse vectorized version
    modV = Model(solver=solver)
    @variable(modV, 0.5 <= x <= 2 )
    @variable(modV, 0 <= y <= 30 )
    obj = [x y]*sparse(ones(2,2))*[x,y]
    @objective(modV, Min, obj[1])
    @constraint(modV, sparse(ones(1,2))*[x,y] .>= 1)
    @fact solve(modV) --> :Optimal
    @fact modV.objVal --> roughly(1.0, 1e-6)
    @fact (getvalue(x) + getvalue(y)) --> roughly(1.0, 1e-6)
end; end; end



facts("[qcqpmodel] Test quad constraints (continuous)") do
for solver in quad_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @variable(modQ, -2 <= x <= 2 )
    @variable(modQ, -2 <= y <= 2 )
    @objective(modQ, Min, x - y )
    @constraint(modQ, x + x*x + x*y + y*y <= 1 )
    @fact MathProgBase.numquadconstr(modQ) --> 1
    @fact MathProgBase.numlinconstr(modQ) --> 0
    @fact MathProgBase.numconstr(modQ) --> 1
    @fact JuMP.isquadsoc(modQ) --> false

    @fact solve(modQ) --> :Optimal
    @fact modQ.objVal --> roughly(-1-4/sqrt(3), 1e-6)
    @fact (getvalue(x) + getvalue(y)) --> roughly(-1/3, 1e-3)

    # Vectorized version
    modV = Model(solver=solver)
    @variable(modV, -2 <= x <= 2 )
    @variable(modV, -2 <= y <= 2 )
    obj = [1 -1]*[x,y]
    @objective(modV, Min, obj[1])
    A = [1 0.5; 0.5 1]
    @constraint(modV, [x y]*A*[x,y] + [x] .<= [1])
    @fact MathProgBase.numquadconstr(modV) --> 1
    @fact MathProgBase.numlinconstr(modV) --> 0
    @fact MathProgBase.numconstr(modV) --> 1

    @fact solve(modV) --> :Optimal
    @fact modV.objVal --> roughly(-1-4/sqrt(3), 1e-6)
    @fact (getvalue(x) + getvalue(y)) --> roughly(-1/3, 1e-3)

    # Sparse vectorized version
    modV = Model(solver=solver)
    @variable(modV, -2 <= x <= 2 )
    @variable(modV, -2 <= y <= 2 )
    obj = [1 -1]*[x,y]
    @objective(modV, Min, obj[1])
    A = sparse([1 0.5; 0.5 1])
    @constraint(modV, [x y]*A*[x,y] + [x] .<= [1])
    @fact MathProgBase.numquadconstr(modV) --> 1
    @fact MathProgBase.numlinconstr(modV) --> 0
    @fact MathProgBase.numconstr(modV) --> 1

    @fact solve(modV) --> :Optimal
    @fact modV.objVal --> roughly(-1-4/sqrt(3), 1e-6)
    @fact (getvalue(x) + getvalue(y)) --> roughly(-1/3, 1e-3)
end; end; end

facts("[qcqpmodel] Test SOC constraints (continuous)") do
for solver in soc_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @variable(modQ, x)
    @variable(modQ, y)
    @variable(modQ, t >= 0)
    @objective(modQ, Min, t)
    @constraint(modQ, x+y >= 1)
    @constraint(modQ, x^2 + y^2 <= t^2)
    @fact JuMP.isquadsoc(modQ) --> true

    @fact solve(modQ) --> :Optimal
    @fact modQ.objVal --> roughly(sqrt(1/2), 1e-6)
    @fact norm([getvalue(x), getvalue(y), getvalue(t)] - [0.5,0.5,sqrt(1/2)]) --> roughly(0.0,1e-3)

    # Vectorized version
    modV = Model(solver=solver)
    @variable(modV, x)
    @variable(modV, y)
    @variable(modV, t >= 0)
    @objective(modV, Min, t)
    @constraint(modV, [1 1 0]*[x,y,t] .>= 1)
    Q = [1 0  0
         0 1  0
         0 0 -1]
    @constraint(modV, [x y t]*Q*[x,y,t] .<= 0)

    @fact solve(modV) --> :Optimal
    @fact modV.objVal --> roughly(sqrt(1/2), 1e-6)
    @fact norm([getvalue(x), getvalue(y), getvalue(t)] - [0.5,0.5,sqrt(1/2)]) --> roughly(0.0,1e-3)

    # Sparse vectorized version
    modV = Model(solver=solver)
    @variable(modV, x)
    @variable(modV, y)
    @variable(modV, t >= 0)
    @objective(modV, Min, t)
    @constraint(modV, sparse([1 1 0])*[x,y,t] .>= 1)
    Q = sparse([1 0  0
                0 1  0
                0 0 -1])
    @constraint(modV, [x y t]*Q*[x,y,t] .<= 0)

    @fact solve(modV) --> :Optimal
    @fact modV.objVal --> roughly(sqrt(1/2), 1e-6)
    @fact norm([getvalue(x), getvalue(y), getvalue(t)] - [0.5,0.5,sqrt(1/2)]) --> roughly(0.0,1e-3)

    # SOC version 1
    modN = Model(solver=solver)
    @variable(modN, x)
    @variable(modN, y)
    @variable(modN, t >= 0)
    @objective(modN, Min, t)
    @constraint(modN, x + y >= 1)
    @constraint(modN, norm([x,y]) <= t)

    @fact solve(modN) --> :Optimal
    @fact modN.objVal --> roughly(sqrt(1/2), 1e-6)
    @fact norm([getvalue(x), getvalue(y), getvalue(t)] - [0.5,0.5,sqrt(1/2)]) --> roughly(0.0,1e-3)

    # SOC version 2
    modN = Model(solver=solver)
    @variable(modN, x)
    @variable(modN, y)
    @variable(modN, t >= 0)
    @objective(modN, Min, t)
    @constraint(modN, x + y >= 1)
    tmp = [x,y]
    @constraint(modN, norm(tmp[i] for i=1:2) <= t)

    @fact solve(modN) --> :Optimal
    @fact modN.objVal --> roughly(sqrt(1/2), 1e-6)
    @fact norm([getvalue(x), getvalue(y), getvalue(t)] - [0.5,0.5,sqrt(1/2)]) --> roughly(0.0,1e-3)
end; end; end

facts("[qcqpmodel] Test SOC duals") do
for solver in soc_solvers
contains("$(typeof(solver))", "MosekSolver") && continue # Mosek doesn't support duals with conic-through-quadratic
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @variable(modQ, x >= 0)
    @variable(modQ, y)
    @variable(modQ, z)
    @objective(modQ, Min, -y-z)
    @constraint(modQ, eq, x <= 1)
    @constraint(modQ, y^2 + z^2 <= x^2)

    @fact solve(modQ) --> :Optimal
    @fact modQ.objVal --> roughly(-sqrt(2), 1e-6)
    @fact getvalue(y) --> roughly(1/sqrt(2), 1e-6)
    @fact getvalue(z) --> roughly(1/sqrt(2), 1e-6)
    @fact getdual(eq) --> roughly(-sqrt(2), 1e-6)

    @objective(modQ, Max, y+z)
    @fact solve(modQ) --> :Optimal
    @fact modQ.objVal --> roughly(sqrt(2), 1e-6)
    @fact getvalue(y) --> roughly(1/sqrt(2), 1e-6)
    @fact getvalue(z) --> roughly(1/sqrt(2), 1e-6)
    @fact getdual(eq) --> roughly(sqrt(2), 1e-6)

    # # SOC syntax version
    # modN = Model(solver=solver)
    # @variable(modN, x >= 0)
    # @variable(modN, y)
    # @variable(modN, z)
    # @objective(modN, Min, -y-z)
    # @constraint(modN, eq, x <= 1)
    # # @constraint(modN, y^2 + z^2 <= x^2)
    # @constraint(modN, (1/4)*norm([4.0y,4.0z]) <= x)
    #
    # @fact solve(modN) --> :Optimal
    # @fact modN.objVal --> roughly(-sqrt(2), 1e-6)
    # @fact getvalue(y) --> roughly(1/sqrt(2), 1e-6)
    # @fact getvalue(z) --> roughly(1/sqrt(2), 1e-6)
    # @fact getdual(eq) --> roughly(-sqrt(2), 1e-6)
    #
    # @objective(modN, Max, y+z)
    # @fact solve(modN) --> :Optimal
    # @fact modN.objVal --> roughly(sqrt(2), 1e-6)
    # @fact getvalue(y) --> roughly(1/sqrt(2), 1e-6)
    # @fact getvalue(z) --> roughly(1/sqrt(2), 1e-6)
    # @fact getdual(eq) --> roughly(sqrt(2), 1e-6)

end; end; end

facts("[qcqpmodel] Test quad constraints (discrete)") do
for solver in quad_mip_solvers
context("With solver $(typeof(solver))") do
    modQ = Model(solver=solver)
    @variable(modQ, -2 <= x <= 2, Int )
    @variable(modQ, -2 <= y <= 2, Int )
    @objective(modQ, Min, x - y )
    @constraint(modQ, x + x*x + x*y + y*y <= 1 )

    @fact solve(modQ) --> :Optimal
    @fact modQ.objVal --> roughly(-3, 1e-6)
    @fact (getvalue(x) + getvalue(y)) --> roughly(-1, 1e-6)

end; end; end

facts("[qcqpmodel] Test simple normed problem") do
for solver in soc_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver);
    @variable(m, x[1:3]);
    @constraint(m, 2norm(x[i]-1 for i=1:3) <= 2)
    @objective(m, Max, x[1]+x[2])

    @fact solve(m) --> :Optimal
    @fact getobjectivevalue(m) --> roughly(2+sqrt(2), 1e-5)
    @fact norm(getvalue(x)-[1+sqrt(1/2),1+sqrt(1/2),1]) --> roughly(0, 1e-6)
    @fact getvalue(norm(x-1)) --> roughly(1, 1e-5)
    @fact getvalue(norm(x-1)-2) --> roughly(-1, 1e-5)
end; end; end

facts("[qcqpmodel] Test quad problem modification") do
for solver in quad_solvers
context("With solver $(typeof(solver))") do

    modQ = Model(solver=solver)
    @variable(modQ, x >= 0)
    @constraint(modQ, x*x <= 1)
    @objective(modQ, Max, x)
    @fact solve(modQ) --> :Optimal
    @fact getobjectivevalue(modQ) --> roughly(1.0, 1e-6)

    @constraint(modQ, 2x*x <= 1)
    @fact modQ.internalModelLoaded --> true
    @fact solve(modQ) --> :Optimal
    @fact getobjectivevalue(modQ) --> roughly(sqrt(0.5), 1e-6)

    modQ = Model(solver=solver)
    @variable(modQ,   0 <= x <= 1)
    @variable(modQ, 1/2 <= y <= 1)
    @objective(modQ, Min, x*x - y)
    @fact solve(modQ) --> :Optimal
    @fact getobjectivevalue(modQ) --> roughly(-1.0, 1e-6)

    @objective(modQ, Min, y*y - x)
    @fact modQ.internalModelLoaded --> true
    @fact solve(modQ) --> :Optimal
    @fact getobjectivevalue(modQ) --> roughly(-0.75, 1e-6)
end; end; end

facts("[qcqpmodel] Rotated second-order cones") do
for solver in rsoc_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)

    @variable(mod, x[1:5] >= 0)
    @variable(mod, 0 <= u <= 5)
    @variable(mod, v)

    @objective(mod, Max, v)

    @constraint(mod, norm(x) <= 1)
    @constraint(mod, v^2 <= u * x[1])

    @fact solve(mod) --> :Optimal
    @fact getvalue(x) --> roughly([1,0,0,0,0], 1e-2)
    @fact getvalue(u) --> roughly(5, 1e-4)
    @fact getvalue(v) --> roughly(sqrt(5), 1e-6)
    @fact getvalue(norm(x)) --> roughly(1, 1e-4)
end; end; end
