#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/quadmodel.jl
# Testing quadratic models (objective and constraints)
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using JuMP, Compat.Test, Compat
# If solvers not loaded, load them (i.e running just these tests)
!isdefined(@__MODULE__, :lp_solvers) && include("solvers.jl")

@testset "Quadratics" begin
    @testset "Quad objective (discrete) with $solver" for solver in quad_mip_solvers

        modQ = Model(solver=solver)
        @variable(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
        @objective(modQ, Min, 10*x[1]*x[1] + 3*x[1]*x[2] + 5*x[2]*x[2] + 9*x[3]*x[3])
        @constraint(modQ, x[2] <= 1.7*x[3])
        @constraint(modQ, x[2] >= 0.5*x[1])
        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), 247.0, atol=1e-5)
        @test isapprox(getvalue(x), [2.0, 3.0, 4.0], atol=1e-6)

        # vectorized version
        modV = Model(solver=solver)
        @variable(modV, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
        obj = x'*[10 1.5 0; 1.5 5 0; 0 0 9]*x
        @objective(modV, Min, obj)
        A = [  0  1 -1.7
             0.5 -1    0]
        @constraint(modV, A*x .<= zeros(2))
        @test solve(modV) == :Optimal
        @test isapprox(getobjectivevalue(modV), 247.0, atol=1e-5)
        @test isapprox(getvalue(x), [2.0, 3.0, 4.0], atol=1e-6)

        modQ = Model(solver=solver)
        @variable(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
        @objective(modQ, :Max, -10*x[1]*x[1] - 3*x[1]*x[2] - 5*x[2]*x[2] - 9*x[3]*x[3])
        @constraint(modQ, x[2] <= 1.7*x[3])
        @constraint(modQ, x[2] >= 0.5*x[1])
        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), -247.0, atol=1e-5)
        @test isapprox(getvalue(x), [2.0, 3.0, 4.0], atol=1e-6)

        # sparse vectorized version
        modV = Model(solver=solver)
        @variable(modV, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
        obj = x'*sparse([10 1.5 0; 1.5 5 0; 0 0 9])*x
        @objective(modV, Min, obj)
        A = sparse([  0  1 -1.7
                    0.5 -1    0])
        @constraint(modV, A*x .<= zeros(2))
        @test solve(modV) == :Optimal
        @test isapprox(getobjectivevalue(modV), 247.0, atol=1e-5)
        @test isapprox(getvalue(x), [2.0, 3.0, 4.0], atol=1e-6)

        modQ = Model(solver=solver)
        @variable(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
        @objective(modQ, :Max, -10*x[1]*x[1] - 3*x[1]*x[2] - 5*x[2]*x[2] - 9*x[3]*x[3])
        @constraint(modQ, x[2] <= 1.7*x[3])
        @constraint(modQ, x[2] >= 0.5*x[1])
        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), -247.0, atol=1e-5)
        @test isapprox(getvalue(x), [2.0, 3.0, 4.0], atol=1e-6)

        # vectorized version
        modV = Model(solver=solver)
        @variable(modV, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
        Q = [10 3 0; 0 5 0; 0 0 9]
        obj = x'Q*x
        @objective(modV, Min, obj)
        A = [   0 -1 1.7
             -0.5  1   0]
        @constraint(modV, A*x .>= zeros(2))
        @test solve(modV) == :Optimal
        @test isapprox(getobjectivevalue(modV), 247.0, atol=1e-5)
        @test isapprox(getvalue(x), [2.0, 3.0, 4.0], atol=1e-6)
    end

    @testset "Quad objective (continuous) with $solver" for solver in quad_solvers

        modQ = Model(solver=solver)
        @variable(modQ, 0.5 <= x <= 2 )
        @variable(modQ, 0 <= y <= 30 )
        @objective(modQ, :Min, (x+y)*(x+y) )
        @constraint(modQ, x + y >= 1 )
        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), 1.0, atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), 1.0, atol=1e-6)

        # Vectorized version
        modV = Model(solver=solver)
        @variable(modV, 0.5 <= x <= 2 )
        @variable(modV, 0 <= y <= 30 )
        obj = [x,y]'ones(2,2)*[x,y]
        @objective(modV, Min, obj)
        @constraint(modV, ones(1,2)*[x,y] .>= 1)
        @test solve(modV) == :Optimal
        @test isapprox(getobjectivevalue(modV), 1.0, atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), 1.0, atol=1e-6)

        # Sparse vectorized version
        modV = Model(solver=solver)
        @variable(modV, 0.5 <= x <= 2 )
        @variable(modV, 0 <= y <= 30 )
        obj = [x y]*sparse(ones(2,2))*[x,y]
        @objective(modV, Min, obj[1])
        @constraint(modV, sparse(ones(1,2))*[x,y] .>= 1)
        @test solve(modV) == :Optimal
        @test isapprox(getobjectivevalue(modV), 1.0, atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), 1.0, atol=1e-6)
    end

    @testset "Quad constraints (continuous) with $solver" for solver in quad_solvers

        modQ = Model(solver=solver)
        @variable(modQ, -2 <= x <= 2 )
        @variable(modQ, -2 <= y <= 2 )
        @objective(modQ, Min, x - y )
        @constraint(modQ, x + x*x + x*y + y*y <= 1 )
        @test MathProgBase.numquadconstr(modQ) == 1
        @test MathProgBase.numlinconstr(modQ) == 0
        @test MathProgBase.numconstr(modQ) == 1

        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), -1-4/sqrt(3), atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), -1/3, atol=1e-3)

        # Vectorized version
        modV = Model(solver=solver)
        @variable(modV, -2 <= x <= 2 )
        @variable(modV, -2 <= y <= 2 )
        obj = [1 -1]*[x,y]
        @objective(modV, Min, obj[1])
        A = [1 0.5; 0.5 1]
        @constraint(modV, [x y]*A*[x,y] + [x] .<= [1])
        @test MathProgBase.numquadconstr(modV) == 1
        @test MathProgBase.numlinconstr(modV) == 0
        @test MathProgBase.numconstr(modV) == 1

        @test solve(modV) == :Optimal
        @test isapprox(getobjectivevalue(modV), -1-4/sqrt(3), atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), -1/3, atol=1e-3)

        # Sparse vectorized version
        modV = Model(solver=solver)
        @variable(modV, -2 <= x <= 2 )
        @variable(modV, -2 <= y <= 2 )
        obj = [1 -1]*[x,y]
        @objective(modV, Min, obj[1])
        A = sparse([1 0.5; 0.5 1])
        @constraint(modV, [x y]*A*[x,y] + [x] .<= [1])
        @test MathProgBase.numquadconstr(modV) == 1
        @test MathProgBase.numlinconstr(modV) == 0
        @test MathProgBase.numconstr(modV) == 1

        @test solve(modV) == :Optimal
        @test isapprox(getobjectivevalue(modV), -1-4/sqrt(3), atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), -1/3, atol=1e-3)
    end

    @testset "SOC constraints (continuous) with $solver" for solver in soc_solvers

        # SOC version 1
        modN = Model(solver=solver)
        @variable(modN, x)
        @variable(modN, y)
        @variable(modN, t >= 0)
        @objective(modN, Min, t)
        @constraint(modN, x + y >= 1)
        @constraint(modN, norm([x,y]) <= t)

        # Getter/setters
        @test JuMP.numsocconstr(modN) == 1
        @test MathProgBase.numconstr(modN) == 2

        @test solve(modN) == :Optimal
        @test isapprox(getobjectivevalue(modN), sqrt(1/2), atol=1e-6)
        @test isapprox([getvalue(x), getvalue(y), getvalue(t)], [0.5,0.5,sqrt(1/2)], atol=1e-3)

        # SOC version 2
        modN = Model(solver=solver)
        @variable(modN, x)
        @variable(modN, y)
        @variable(modN, t >= 0)
        @objective(modN, Min, t)
        @constraint(modN, x + y >= 1)
        tmp = [x,y]
        @constraint(modN, norm(tmp[i] for i=1:2) <= t)

        @test solve(modN) == :Optimal
        @test isapprox(getobjectivevalue(modN), sqrt(1/2), atol=1e-6)
        @test isapprox([getvalue(x), getvalue(y), getvalue(t)], [0.5,0.5,sqrt(1/2)], atol=1e-3)
    end

    @testset "SOC constraints (continuous) in quadratic form with $solver" for solver in quad_soc_solvers

        modQ = Model(solver=solver)
        @variable(modQ, x)
        @variable(modQ, y)
        @variable(modQ, t >= 0)
        @objective(modQ, Min, t)
        @constraint(modQ, x+y >= 1)
        @constraint(modQ, x^2 + y^2 <= t^2)

        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), sqrt(1/2), atol=1e-6)
        @test isapprox([getvalue(x), getvalue(y), getvalue(t)], [0.5,0.5,sqrt(1/2)], atol=1e-3)

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

        @test solve(modV) == :Optimal
        @test isapprox(getobjectivevalue(modV), sqrt(1/2), atol=1e-6)
        @test isapprox([getvalue(x), getvalue(y), getvalue(t)], [0.5,0.5,sqrt(1/2)], atol=1e-3)

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

        @test solve(modV) == :Optimal
        @test isapprox(getobjectivevalue(modV), sqrt(1/2), atol=1e-6)
        @test isapprox([getvalue(x), getvalue(y), getvalue(t)], [0.5,0.5,sqrt(1/2)], atol=1e-3)
    end

    @testset "SOC duals with $solver" for solver in soc_solvers
        occursin("MosekSolver", "$(typeof(solver))") && continue # Mosek doesn't support duals with conic-through-quadratic

        modQ = Model(solver=solver)
        @variable(modQ, x >= 0)
        @variable(modQ, y)
        @variable(modQ, z)
        @objective(modQ, Min, -y-z)
        @constraint(modQ, eq_con, x <= 1)
        @constraint(modQ, y^2 + z^2 <= x^2)

        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), -sqrt(2), atol=1e-6)
        @test isapprox(getvalue(y), 1/sqrt(2), atol=1e-6)
        @test isapprox(getvalue(z), 1/sqrt(2), atol=1e-6)
        @test isapprox(getdual(eq_con), -sqrt(2), atol=1e-6)

        @objective(modQ, Max, y+z)
        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), sqrt(2), atol=1e-6)
        @test isapprox(getvalue(y), 1/sqrt(2), atol=1e-6)
        @test isapprox(getvalue(z), 1/sqrt(2), atol=1e-6)
        @test isapprox(getdual(eq_con), sqrt(2), atol=1e-6)

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
        # @test solve(modN) == :Optimal
        # @test modN.objVal == roughly(-sqrt(2), 1e-6)
        # @test getvalue(y) == roughly(1/sqrt(2), 1e-6)
        # @test getvalue(z) == roughly(1/sqrt(2), 1e-6)
        # @test getdual(eq) == roughly(-sqrt(2), 1e-6)
        #
        # @objective(modN, Max, y+z)
        # @test solve(modN) == :Optimal
        # @test modN.objVal == roughly(sqrt(2), 1e-6)
        # @test getvalue(y) == roughly(1/sqrt(2), 1e-6)
        # @test getvalue(z) == roughly(1/sqrt(2), 1e-6)
        # @test getdual(eq) == roughly(sqrt(2), 1e-6)

    end

    @testset "Quad constraints (discrete) with $solver" for solver in quad_mip_solvers
        modQ = Model(solver=solver)
        @variable(modQ, -2 <= x <= 2, Int )
        @variable(modQ, -2 <= y <= 2, Int )
        @objective(modQ, Min, x - y )
        @constraint(modQ, x + x*x + x*y + y*y <= 1 )

        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), -3, atol=1e-6)
        @test isapprox(getvalue(x) + getvalue(y), -1, atol=1e-6)

    end

    @testset "Simple normed problem with $solver" for solver in soc_solvers
        m = Model(solver=solver);
        @variable(m, x[1:3]);
        @constraint(m, 2norm(x[i]-1 for i=1:3) <= 2)
        @objective(m, Max, x[1]+x[2])

        @test solve(m) == :Optimal
        @test isapprox(getobjectivevalue(m), 2+sqrt(2), atol=1e-5)
        @test isapprox(getvalue(x), [1+sqrt(1/2),1+sqrt(1/2),1], atol=1e-6)
        @test isapprox(getvalue(norm(x-1)), 1, atol=1e-5)
        @test isapprox(getvalue(norm(x-1)-2), -1, atol=1e-5)
    end

    @testset "Quad problem modification with $solver" for solver in quad_solvers

        modQ = Model(solver=solver)
        @variable(modQ, x >= 0)
        @constraint(modQ, x*x <= 1)
        @objective(modQ, Max, x)
        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), 1.0, atol=1e-6)

        @constraint(modQ, 2x*x <= 1)
        @test modQ.internalModelLoaded == true
        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), sqrt(0.5), atol=1e-6)

        modQ = Model(solver=solver)
        @variable(modQ,   0 <= x <= 1)
        @variable(modQ, 1/2 <= y <= 1)
        @objective(modQ, Min, x*x - y)
        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), -1.0, atol=1e-6)

        @objective(modQ, Min, y*y - x)
        @test modQ.internalModelLoaded == true
        @test solve(modQ) == :Optimal
        @test isapprox(getobjectivevalue(modQ), -0.75, atol=1e-6)
    end

    @testset "Rotated second-order cones with $solver" for solver in quad_soc_solvers
        mod = Model(solver=solver)

        @variable(mod, x[1:5] >= 0)
        @variable(mod, 0 <= u <= 5)
        @variable(mod, v)
        @variable(mod, t1 == 1)
        @variable(mod, t2 == 1)

        @objective(mod, Max, v)

        @constraint(mod, sum(x.^2) <= t1*t2)
        @constraint(mod, v^2 <= u * x[1])

        @test solve(mod) == :Optimal
        @test isapprox(getvalue(x), [1,0,0,0,0], atol=1e-2)
        @test isapprox(getvalue(u), 5, atol=1e-4)
        @test isapprox(getvalue(v), sqrt(5), atol=1e-6)
        @test isapprox(getvalue(norm(x)), 1, atol=1e-4)
    end
end
