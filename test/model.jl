#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/model.jl
# Testing Model printing, basic solving
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using Compat, Compat.LinearAlgebra, Compat.SparseArrays, Compat.Test
using MathProgBase, JuMP
using OffsetArrays

# If solvers not loaded, load them (i.e running just these tests)
!isdefined(@__MODULE__, :lp_solvers) && include("solvers.jl")

# To ensure the tests work on Windows and Linux/OSX, we need
# to use the correct comparison operators
const leq = JuMP.repl[:leq]
const geq = JuMP.repl[:geq]
const  eq = JuMP.repl[:eq]

const TOL = 1e-4

modPath = joinpath(dirname(@__FILE__), "mod")

mutable struct TemporaryExtensionTestType
    x::Int
end

@testset "Models" begin

    @testset "Check error cases" begin
        @test_throws ErrorException Model(solver=:Foo)
        modErr = Model()
        @test_throws ErrorException setobjectivesense(modErr, :Maximum)
        @variable(modErr, errVar)
        @test isnan(getvalue(errVar))
        @test isnan(getdual(errVar))
        @test_throws ErrorException solve(modErr) # no solver provided

        if length(ip_solvers) > 0
            modErr = Model(solver=ip_solvers[1])
            @variable(modErr, x, Bin)
            @objective(modErr, Max, x)
            con = @constraint(modErr, x <= 0.5)
            solve(modErr)
            @test isnan(getdual(con))
        end

        if length(lp_solvers) > 0
            modErr = Model(solver=lp_solvers[1])
            @variable(modErr, 0 <= x <= 1)
            @objective(modErr, Max, x)
            @constraint(modErr, x <= -1)
            solve(modErr, suppress_warnings=true)
            @test isnan(getvalue(x))
        end
    end

    @testset "Warning on non-symbol variable names" begin
        m = Model()
        x = Dict()
        @variable(m, x[1][1])
    end

    @testset "Performance warnings" begin
        m = Model()
        @variable(m, x[1:2], start=0)
        for i in 1:500
           getvalue(x)
        end
        q = 0
        for i in 1:30000
           q += 3x[1]
        end
    end

    @testset "Test printing a model" begin
        modA = Model()
        x = Variable(modA, 0, Inf, :Cont)
        @variable(modA, y <= 5, Int)
        @variable(modA, 2 <= z <= 4)
        @variable(modA, 0 <= r[i=3:6] <= i)
        @objective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
        @constraintref constraints[1:3]
        constraints[1] = @constraint(modA, 2 <= x+y <= 4)
        constraints[2] = @constraint(modA, sum(r[i] for i=3:5) <= (2 - x)/2.0)
        constraints[3] = @constraint(modA, 6y + y <= z + r[6]/1.9)
        addSOS1(modA, [x, 2y, 3z])
        addSOS2(modA, [30x, 20y, 10z])
        #####################################################################
        # Test LP writer (given names)
        writeLP(modA, modPath * "A.lp", genericnames=false)
        modALP = String[
                 "Maximize",
                 "obj: 0.16666666666666666 col_1 + 0.16666666666666666 y + 1 z + 1 r_3",
                 "Subject To",
                 "c1: 1 col_1 + 1 y >= 2",
                 "c2: 1 col_1 + 1 y <= 4",
                 "c3: 1 r_3 + 1 r_4 + 1 r_5 + 0.5 col_1 <= 1",
                 "c4: 6 y + 1 y - 1 z - 0.5263157894736842 r_6 <= 0",
                 "c5: S1:: col_1:1 y:2 z:3",
                 "c6: S2:: col_1:30 y:20 z:10",
                 "Bounds",
                 "0 <= col_1 <= +inf",
                 "-inf <= y <= 5",
                 "2 <= z <= 4",
                 "0 <= r_3 <= 3",
                 "0 <= r_4 <= 4",
                 "0 <= r_5 <= 5",
                 "0 <= r_6 <= 6",
                 "General",
                 "y",
                 "End"]
        modAfp = open(modPath * "A.lp")
        lineInd = 1
        while !eof(modAfp)
            line = readline(modAfp)
            @test strip(line) == strip(modALP[lineInd])
            lineInd += 1
        end
        close(modAfp)
        #####################################################################
        # Test LP writer (generic names)
        writeLP(modA, modPath * "A.lp", genericnames=true)
        modALP = String[
        "Maximize",
        "obj: 0.16666666666666666 VAR1 + 0.16666666666666666 VAR2 + 1 VAR3 + 1 VAR4",
        "Subject To",
        "c1: 1 VAR1 + 1 VAR2 >= 2",
        "c2: 1 VAR1 + 1 VAR2 <= 4",
        "c3: 1 VAR4 + 1 VAR5 + 1 VAR6 + 0.5 VAR1 <= 1",
        "c4: 6 VAR2 + 1 VAR2 - 1 VAR3 - 0.5263157894736842 VAR7 <= 0",
        "c5: S1:: VAR1:1 VAR2:2 VAR3:3",
        "c6: S2:: VAR1:30 VAR2:20 VAR3:10",
        "Bounds",
        "0 <= VAR1 <= +inf",
        "-inf <= VAR2 <= 5",
        "2 <= VAR3 <= 4",
        "0 <= VAR4 <= 3",
        "0 <= VAR5 <= 4",
        "0 <= VAR6 <= 5",
        "0 <= VAR7 <= 6",
        "General",
        "VAR2",
        "End"]

        modAfp = open(modPath * "A.lp")
        lineInd = 1
        while !eof(modAfp)
            line = readline(modAfp)
            @test strip(line) == strip(modALP[lineInd])
            lineInd += 1
        end
        close(modAfp)
        #####################################################################
        # Test MPS writer
        writeMPS(modA, modPath * "A.mps")
        modAMPS = String[
        "NAME   JuMPModel",
        "ROWS",
        " N  OBJ",
        " E  CON1",
        " L  CON2",
        " L  CON3",
        "COLUMNS",
        "    VAR1  CON1  1",
        "    VAR1  CON2  0.5",
        "    VAR1  OBJ  -0.16666666666666666",
        "    MARKER    'MARKER'                 'INTORG'",
        "    VAR2  CON1  1",
        "    VAR2  CON3  7",
        "    VAR2  OBJ  -0.16666666666666666",
        "    MARKER    'MARKER'                 'INTEND'",
        "    VAR3  CON3  -1",
        "    VAR3  OBJ  -1",
        "    VAR4  CON2  1",
        "    VAR4  OBJ  -1",
        "    VAR5  CON2  1",
        "    VAR5  OBJ  0",
        "    VAR6  CON2  1",
        "    VAR6  OBJ  0",
        "    VAR7  CON3  -0.5263157894736842",
        "    VAR7  OBJ  0",
        "RHS",
        "    rhs    CON1    2",
        "    rhs    CON2    1",
        "    rhs    CON3    0",
        "RANGES",
        "    rhs    CON1    2",
        "BOUNDS",
        "  PL BOUND VAR1",
        "  MI BOUND VAR2",
        "  UP BOUND VAR2 5",
        "  LO BOUND VAR3 2",
        "  UP BOUND VAR3 4",
        "  UP BOUND VAR4 3",
        "  UP BOUND VAR5 4",
        "  UP BOUND VAR6 5",
        "  UP BOUND VAR7 6",
        "ENDATA"]

        modAfp = open(modPath * "A.mps")
        lineInd = 1
        while !eof(modAfp)
            line = readline(modAfp)
            @test chomp(line) == modAMPS[lineInd]
            lineInd += 1
        end
        close(modAfp)

        # Getter/setters
        @test MathProgBase.numvar(modA) == 7
        @test MathProgBase.numlinconstr(modA) == 3
        @test MathProgBase.numquadconstr(modA) == 0
        @test MathProgBase.numconstr(modA) == 5
        @test getobjectivesense(modA) == :Max
        setobjectivesense(modA, :Min)
        @test getobjectivesense(modA) == :Min
    end

    @testset "Quadratic MPS writer" begin
        modQ = Model()
        @variable(modQ, x == 1)
        @variable(modQ, y ≥ 0)
        @objective(modQ, Min, x^2 - 2*x*y + y^2 + x)

        #####################################################################
        # Test MPS writer
        writeMPS(modQ, modPath * "Q.mps")
        modQMPS = String[
        "NAME   JuMPModel",
        "ROWS",
        " N  OBJ",
        "COLUMNS",
        "    VAR1  OBJ  1",
        "    VAR2  OBJ  0",
        "RHS",
        "BOUNDS",
        "  LO BOUND VAR1 1",
        "  UP BOUND VAR1 1",
        "  PL BOUND VAR2",
        "QMATRIX",
        "  VAR1 VAR1  2",
        "  VAR1 VAR2 -2",
        "  VAR2 VAR2  2",
        "ENDATA"]
        modQfp = open(modPath * "Q.mps")
        lineInd = 1
        while !eof(modQfp)
            line = readline(modQfp)
            @test chomp(line) == modQMPS[lineInd]
            lineInd += 1
        end
        close(modQfp)

        # Getter/setters
        @test MathProgBase.numvar(modQ) == 2
        @test MathProgBase.numlinconstr(modQ) == 0
        @test MathProgBase.numquadconstr(modQ) == 0
        @test MathProgBase.numconstr(modQ) == 0
        @test getobjectivesense(modQ) == :Min
    end


    @testset "Solving a MILP with $solver" for solver in ip_solvers
        modA = Model(solver=solver)
        @variable(modA, x >= 0)
        @variable(modA, y <= 5, Int)
        @variable(modA, 2 <= z <= 4)
        @variable(modA, 0 <= r[i=3:6] <= i)
        @objective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
        @constraint(modA, 2 <= x+y)
        @constraint(modA, sum(r[i] for i=3:5) <= (2 - x)/2.0)
        @constraint(modA, 7.0*y <= z + r[6]/1.9)

        @test solve(modA)       == :Optimal
        @test isapprox(modA.objVal, 1.0+4.833334, atol=TOL)
        @test isapprox(getvalue(x), 1.0, atol=TOL)
        @test isapprox(getvalue(y), 1.0, atol=TOL)
        @test isapprox(getvalue(z), 4.0, atol=TOL)
        @test isapprox(getvalue(r)[3], 0.5, atol=TOL)
        @test isapprox(getvalue(r)[4], 0.0, atol=TOL)
        @test isapprox(getvalue(r)[5], 0.0, atol=TOL)
        @test getvalue(r)[6] ≥ 5.7-TOL
        @test getvalue(r)[6] ≤ 6.0+TOL
        @test getobjective(modA).aff == ((x + y)/2.0 + 3.0)/3.0 + z + r[3]
    end

    @testset "Solving an LP (Min) with $solver" for solver in lp_solvers
        modA = Model(solver=solver)
        @variable(modA, x >= 0)
        @variable(modA, y <= 5)
        @variable(modA, 2 <= z <= 4)
        @variable(modA, 0 <= r[i=3:6] <= i)
        @objective(modA, Min, -((x + y)/2.0 + 3.0)/3.0 - z - r[3])
        @constraintref cons[1:3]
        cons[1] = @constraint(modA, x+y >= 2)
        cons[2] = @constraint(modA, sum(r[i] for i=3:5) <= (2 - x)/2.0)
        cons[3] = @constraint(modA, 7.0*y <= z + r[6]/1.9)

        # Solution
        @test solve(modA) == :Optimal
        @test isapprox(getobjectivevalue(modA), -5.8446115, atol=TOL)
        @test isapprox(getvalue(x), 0.9774436, atol=TOL)
        @test isapprox(getvalue(y), 1.0225563, atol=TOL)
        @test isapprox(getvalue(z), 4.0, atol=TOL)
        @test isapprox(getvalue(r)[3], 0.5112781, atol=TOL)
        @test isapprox(getvalue(r)[4], 0.0, atol=TOL)
        @test isapprox(getvalue(r)[5], 0.0, atol=TOL)
        @test isapprox(getvalue(r)[6], 6.0, atol=TOL)

        # Reduced costs
        @test isapprox(getdual(x), 0.0, atol=TOL)
        @test isapprox(getdual(y), 0.0, atol=TOL)
        @test isapprox(getdual(z), -1.0714286, atol=TOL)
        @test isapprox(getdual(r)[3], 0.0, atol=TOL)
        @test isapprox(getdual(r)[4], 1.0, atol=TOL)
        @test isapprox(getdual(r)[5], 1.0, atol=TOL)
        @test isapprox(getdual(r)[6], -0.03759398, atol=TOL)

        # Row duals
        @test isapprox(getdual(cons)[1], 0.333333, atol=TOL)
        @test isapprox(getdual(cons)[2], -1.0, atol=TOL)
        @test isapprox(getdual(cons)[3], -0.0714286, atol=TOL)
    end

    @testset "Test solving an LP (Max) with $solver" for solver in lp_solvers
        modA = Model(solver=solver)
        @variable(modA, x >= 0)
        @variable(modA, y <= 5)
        @variable(modA, 2 <= z <= 4)
        @variable(modA, 0 <= r[i=3:6] <= i)
        @objective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
        @constraintref cons[1:3]
        cons[1] = @constraint(modA, x+y >= 2)
        cons[2] = @constraint(modA, sum(r[i] for i=3:5) <= (2 - x)/2.0)
        cons[3] = @constraint(modA, 7.0*y <= z + r[6]/1.9)

        # Solution
        @test solve(modA) == :Optimal
        @test isapprox(getobjectivevalue(modA), 5.8446115, atol=TOL)
        @test isapprox(getvalue(x), 0.9774436, atol=TOL)
        @test isapprox(getvalue(y), 1.0225563, atol=TOL)
        @test isapprox(getvalue(z), 4.0, atol=TOL)
        @test isapprox(getvalue(r)[3], 0.5112781, atol=TOL)
        @test isapprox(getvalue(r)[4], 0.0, atol=TOL)
        @test isapprox(getvalue(r)[5], 0.0, atol=TOL)
        @test isapprox(getvalue(r)[6], 6.0, atol=TOL)

        # Reduced costs
        @test isapprox(getdual(x), 0.0, atol=TOL)
        @test isapprox(getdual(y), 0.0, atol=TOL)
        @test isapprox(getdual(z), 1.0714286, atol=TOL)
        @test isapprox(getdual(r)[3], 0.0, atol=TOL)
        @test isapprox(getdual(r)[4], -1.0, atol=TOL)
        @test isapprox(getdual(r)[5], -1.0, atol=TOL)
        @test isapprox(getdual(r)[6], 0.03759398, atol=TOL)

        # Row duals
        @test isapprox(getdual(cons)[1], -0.333333, atol=TOL)
        @test isapprox(getdual(cons)[2], 1.0, atol=TOL)
        @test isapprox(getdual(cons)[3], 0.0714286, atol=TOL)
    end



    @testset "Test binary variable handling with $solver" for solver in ip_solvers
        modB = Model(solver=solver)
        @variable(modB, x, Bin)
        @objective(modB, Max, x)
        @constraint(modB, x <= 10)
        status = solve(modB)
        @test status == :Optimal
        @test isapprox(getvalue(x), 1.0, atol=TOL)
    end


    @testset "Test model copying" begin
        source = Model()
        @variable(source, 2 <= x <= 5)
        @variable(source, 0 <= y <= 1, Int)
        @objective(source, Max, 3*x + 1*y)
        @constraint(source, x + 2.0*y <= 6)
        @constraint(source, x*x <= 1)
        addSOS2(source, [x, 2y])
        @SDconstraint(source, x*ones(3,3) >= Matrix(1.0I, 3, 3))
        @SDconstraint(source, ones(3,3) <= 0)
        @variable(source, z[1:3])
        @variable(source, w[2:4]) # JuMPArray
        @variable(source, v[[:red],i=1:3;isodd(i)]) # JuMPDict
        JuMP.setsolvehook(source, m -> :Optimal)

        # uncomment when NLP copying is implemented
        # @NLconstraint(source, c[k=1:3], x^2 + y^3 * sin(x+k*y) >= 1)

        dest = copy(source)

        # uncomment when NLP copying is implemented
        # for name in source.nlpdata
        #     @test source.name == dest.name == true
        # end

        # Obj
        @objective(source, Max, 1x)
        @test length(source.obj.aff.coeffs) == 1
        @test length(dest.obj.aff.coeffs) == 2
        @objective(dest, Max, 1x)
        @test length(source.obj.aff.coeffs) == 1
        @test length(dest.obj.aff.coeffs) == 1
        @objective(dest, Max, 3*x + 1*y)
        @test length(source.obj.aff.coeffs) == 1
        @test length(dest.obj.aff.coeffs) == 2

        # Constraints
        source.linconstr[1].ub = 5.0
        @test dest.linconstr[1].ub == 6.0
        source.quadconstr[1].sense = :(>=)
        @test dest.quadconstr[1].sense == :(<=)
        source.sosconstr[1].terms[1] = Variable(source, 2)
        source.sosconstr[1].terms[2] = Variable(source, 1)
        source.sosconstr[1].weights = [2.0,1.0]
        source.sosconstr[1].sostype = :SOS1
        @test dest.sosconstr[1].terms[1].col == 1
        @test dest.sosconstr[1].terms[2].col == 2
        @test dest.sosconstr[1].weights == [1.0,2.0]
        @test dest.sosconstr[1].sostype == :SOS2
        @test length(dest.sdpconstr) == 2
        xx = copy(x, dest)
        @test all(t -> isequal(t[1],t[2]), zip(dest.sdpconstr[1].terms, xx*ones(3, 3) - Matrix(1.0I, 3, 3)))
        @test all(t -> isequal(t[1],t[2]), zip(dest.sdpconstr[2].terms, convert(Matrix{AffExpr}, -ones(3,3))))

        @test dest.solvehook(dest) == :Optimal

        @test Set(collect(keys(dest.objDict))) == Set([:x,:y,:z,:w,:v])
        @test isequal(dest.objDict[:x], Variable(dest, 1))
        @test isequal(dest.objDict[:y], Variable(dest, 2))
        @test all(t -> isequal(t[1], t[2]), zip(dest.objDict[:z], [Variable(dest, 3), Variable(dest, 4), Variable(dest, 5)]))
        @test all(t -> isequal(t[1], t[2]), zip(dest.objDict[:w].innerArray, [Variable(dest, 6), Variable(dest, 7), Variable(dest, 8)]))
        td = dest.objDict[:v].tupledict
        @test length(td) == 2
        @test isequal(td[:red,1], Variable(dest, 9))
        @test isequal(td[:red,3], Variable(dest, 10))

        # Issue #358
        @test typeof(dest.linconstr)  == Array{JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}},1}
        @test typeof(dest.quadconstr) == Array{JuMP.GenericQuadConstraint{JuMP.GenericQuadExpr{Float64,JuMP.Variable}},1}

        JuMP.setprinthook(source, m -> 2)
        dest2 = copy(source)
        @test dest2.printhook(dest2) == 2

        addlazycallback(source, cb -> 3)
        @test_throws ErrorException copy(source)
    end


    @testset "Test extension copy" begin
        source = Model()
        source.ext[:extensiontype] = TemporaryExtensionTestType(1)
        @test_throws ErrorException copy(source)
        source.ext[:extensiontype] = 1
        dest = copy(source)
        source.ext[:extensiontype] = 2
        @test haskey(dest.ext, :extensiontype)
        @test dest.ext[:extensiontype] == 1
    end


    @testset "Test variable/model 'hygiene'" begin
    @testset "Linear constraint" begin
        modA = Model(); @variable(modA, x)
        modB = Model(); @variable(modB, y)
        @constraint(modA, x+y == 1)
        @test_throws ErrorException solve(modA)
    end
    @testset "Linear objective" begin
        modA = Model(); @variable(modA, 0 <= x <= 1)
        modB = Model(); @variable(modB, 0 <= y <= 1)
        @objective(modA, Min, x + y)
        @test_throws ErrorException solve(modA)
    end
    @testset "Quadratic constraint" begin
        modA = Model(); @variable(modA, x)
        modB = Model(); @variable(modB, y)
        @constraint(modB, x*y >= 1)
        @test_throws ErrorException solve(modB)
    end
    @testset "Affine in quadratic constraint" begin
        modA = Model(); @variable(modA, x)
        modB = Model(); @variable(modB, y)
        @constraint(modB, y*y + x + y <= 1)
        @test_throws ErrorException solve(modB)
    end
    @testset "Quadratic objective" begin
        modA = Model(); @variable(modA, x)
        modB = Model(); @variable(modB, y)
        @objective(modB, Min, x*y)
        @test_throws ErrorException solve(modB)
    end
    end

    if length(lp_solvers) > 0
        @testset "Test NaN checking" begin
            mod = Model(solver=lp_solvers[1])
            @variable(mod, x)
            @objective(mod, Min, NaN*x)
            @test_throws ErrorException solve(mod)
            @objective(mod, Min, NaN*x^2)
            @test_throws ErrorException solve(mod)
            @objective(mod, Min, x)
            @constraint(mod, Min, NaN*x == 0)
            @test_throws ErrorException solve(mod)
        end
    end

    @testset "Test column-wise modeling with $lp_solver" for lp_solver in lp_solvers
        mod = Model(solver=lp_solver)
        @variable(mod, 0 <= x <= 1)
        @variable(mod, 0 <= y <= 1)
        @objective(mod, Max, 5x + 1y)
        @constraint(mod, con[i=1:2], i*x + y <= i+5)
        @variable(mod, 0 <= z1 <= 1, objective=10.0, inconstraints=con, coefficients=[1.0,-2.0])
        @variable(mod, 0 <= z2 <= 1, objective=10.0, inconstraints=Any[con[i] for i in 1:2], coefficients=[1.0,-2.0])
        @test solve(mod) == :Optimal
        @test isapprox(getvalue(z1), 1.0, atol=TOL)
        @test isapprox(getvalue(z2), 1.0, atol=TOL)

        # do a vectorized version as well
        mod = Model(solver=lp_solver)
        @variable(mod, 0 <= x <= 1)
        @variable(mod, 0 <= y <= 1)
        obj = [5,1]'*[x,y]
        @objective(mod, Max, obj)
        A = [1 1
             2 1]
        @constraint(mod, A*[x,y] .<= [6,7])
        @variable(mod, 0 <= z1 <= 1, objective=10.0, inconstraints=con, coefficients=[1.0,-2.0])
        @variable(mod, 0 <= z2 <= 1, objective=10.0, inconstraints=Any[con[i] for i in 1:2], coefficients=[1.0,-2.0])
        @test solve(mod) == :Optimal
        @test isapprox(getvalue(z1), 1.0, atol=TOL)
        @test isapprox(getvalue(z2), 1.0, atol=TOL)

        # vectorized with sparse matrices
        mod = Model(solver=lp_solver)
        @variable(mod, 0 <= x <= 1)
        @variable(mod, 0 <= y <= 1)
        # TODO: get this to work by adding method for Ac_mul_B!
        # obj = sparse([5,1])'*[x,y]
        obj = sparse([5,1]')*[x,y]
        @objective(mod, Max, obj[1])
        A = sparse([1 1
                    2 1])
        @constraint(mod, A*[x,y] .<= [6,7])
        @variable(mod, 0 <= z1 <= 1, objective=10.0, inconstraints=con, coefficients=[1.0,-2.0])
        @variable(mod, 0 <= z2 <= 1, objective=10.0, inconstraints=Any[con[i] for i in 1:2], coefficients=[1.0,-2.0])
        @test solve(mod) == :Optimal
        @test isapprox(getvalue(z1), 1.0, atol=TOL)
        @test isapprox(getvalue(z2), 1.0, atol=TOL)
    end

    @testset "[model] Test all MPS paths" begin
        mod = Model()
        @variable(mod, free_var)
        @variable(mod, int_var, Bin)
        @variable(mod, low_var >= 5)
        @constraint(mod, free_var == int_var)
        @constraint(mod, free_var - int_var >= 0)
        @objective(mod, Max, free_var*int_var + low_var)
        writeMPS(mod,"test.mps")
    end

    @testset "Test all LP paths" begin
        mod = Model()
        @variable(mod, free_var)
        @objective(mod, Max, free_var*free_var)
        @test_throws ErrorException writeLP(mod,"test.lp")
        @objective(mod, Max, free_var)
        @constraint(mod, free_var - 2*free_var == 0)
        @constraint(mod, free_var + 2*free_var >= 1)
        writeLP(mod,"test.lp")
    end

    @testset "Test semi-continuous variables with $solver" for solver in semi_solvers
        mod = Model(solver=solver)
        @variable(mod, x >= 3, SemiCont)
        @variable(mod, y >= 2, SemiCont)
        @constraint(mod, x + y >= 1)
        @objective(mod, Min, x+y)
        solve(mod)
        @test isapprox(getvalue(x), 0.0, atol=TOL)
        @test isapprox(getvalue(y), 2.0, atol=TOL)
    end

    @testset "Test semi-integer variables with $solver" for solver in semi_solvers
        mod = Model(solver=solver)
        @variable(mod, x >= 3, SemiInt)
        @variable(mod, y >= 2, SemiInt)
        @constraint(mod, x + y >= 2.5)
        @objective(mod, Min, x+1.1y)
        solve(mod)
        @test isapprox(getvalue(x), 3.0, atol=TOL)
        @test isapprox(getvalue(y), 0.0)
    end

    @testset "Test fixed variables don't leak through MPB with $solver" for solver in lp_solvers
        mod = Model(solver=solver)
        @variable(mod, 0 <= x[1:3] <= 2)
        @variable(mod, y[k=1:2] == k)
        @objective(mod, Min, x[1] + x[2] + x[3] + y[1] + y[2])
        solve(mod)
        for i in 1:3
            @test isapprox(getvalue(x[i]), 0, atol=TOL)
        end
        for k in 1:2
            @test isapprox(getvalue(y[k]), k, atol=TOL)
        end
    end
    @testset "Fixed variables with $solver" for solver in ip_solvers
        mod = Model(solver=solver)
        @variable(mod, x[1:3], Bin)
        @variable(mod, y[k=1:2] == k)
        JuMP.build(mod)
        @test MathProgBase.getvartype(internalmodel(mod)) == [:Bin,:Bin,:Bin,:Cont,:Cont]
    end


    @testset "Test SOS constraints with $solver" for solver in sos_solvers
        modS = Model(solver=solver)

        @variable(modS, x[1:3], Bin)
        @variable(modS, y[1:5], Bin)
        @variable(modS, z)
        @variable(modS, w)

        @objective(modS, Max, z+w)

        a = [1,2,3]
        b = [5,4,7,2,1]

        @constraint(modS, z == sum(a[i]*x[i] for i=1:3))
        @constraint(modS, w == sum(b[i]*y[i] for i=1:5))

        @test_throws MethodError addSOS1([x[1]+y[1]])
        @test_throws MethodError addSOS1([1z])

        addSOS1(modS, [a[i]x[i] for i in 1:3])
        addSOS2(modS, [b[i]y[i] for i in 1:5])

        @test_throws ErrorException addSOS1(modS, [x[1], x[1]+x[2]])

        @test solve(modS) == :Optimal
        @test isapprox(modS.objVal, 15.0, atol=TOL)
        @test isapprox(getvalue(z),  3.0, atol=TOL)
        @test isapprox(getvalue(w), 12.0, atol=TOL)


        m = Model(solver=solver)
        ub = [1,1,2]
        @variable(m, x[i=1:3] <= ub[i])
        @objective(m, Max, 2x[1]+x[2]+x[3])
        addSOS1(m, [x[1],2x[2]])
        addSOS1(m, [x[1],2x[3]])

        # Getter/setters
        @test JuMP.numsosconstr(m) == 2

        @test solve(m) == :Optimal
        @test isapprox(getvalue(x), [0.0,1.0,2.0], atol=TOL)
        @test isapprox(getobjectivevalue(m), 3.0, atol=TOL)
    end

    @testset "Test vectorized model creation" begin

        A = sprand(50,10,0.15)
        B = sprand(50, 7,0.2)
        modV = Model()
        @variable(modV, x[1:10])
        @variable(modV, y[1:7])
        @constraint(modV, A*x + B*y .<= 1)
        obj = (x'*2A')*(2A*x) + (B*2y)'*(B*(2y))
        @objective(modV, Max, obj)

        modS = Model()
        @variable(modS, x[1:10])
        @variable(modS, y[1:7])
        for i in 1:50
            @constraint(modS, sum(A[i,j]*x[j] for j=1:10) + sum(B[i,k]*y[k] for k=1:7) <= 1)
        end
        AA, BB = 4A'*A, 4B'*A
        @objective(modS, Max, sum(AA[i,j]*x[i]*x[j] for i=1:10,j=1:10) + sum(BB[i,j]*y[i]*y[j] for i=1:7, j=1:7))

        @test JuMP.prepConstrMatrix(modV) == JuMP.prepConstrMatrix(modS)
        @test JuMP.prepAffObjective(modV) == JuMP.prepAffObjective(modS)
        @test JuMP.prepConstrBounds(modV) == JuMP.prepConstrBounds(modS)
    end

    @testset "Test MIQP vectorization" begin
        n = 1000
        p = 4
        function bestsubset(solver,X,y,K,M,integer)
            mod = Model(solver=solver)
            @variable(mod, β[1:p])
            if integer
                @variable(mod, z[1:p], Bin)
            else
                @variable(mod, 0 <= z[1:p] <= 1)
            end
            obj = (y-X*β)'*(y-X*β)
            @objective(mod, Min, obj)
            @constraint(mod, I*β .<=  M*I*z)
            @constraint(mod, I*β .>= -M*I*z)
            @constraint(mod, sum(z) == K)
            solve(mod)
            return getvalue(β)
        end
        include(joinpath("data","miqp_vector.jl")) # loads X and q
        y = X * [100, 50, 10, 1] + 20*q
        @testset "with $solver" for solver in quad_solvers
            @test isapprox(bestsubset(solver,X,y,2,500,false), [101.789,49.414,8.63904,1.72663], atol=10TOL)
            @test isapprox(bestsubset(solver,sparse(X),y,2,500,false), [101.789,49.414,8.63904,1.72663], atol=10TOL)
        end
        @testset "with $solver" for solver in quad_mip_solvers
            y = X * [100, 50, 10, 1] + 20*q
            @test isapprox(bestsubset(solver,X,y,2,500,true), [106.25,53.7799,0.0,0.0], atol=1.0)
            @test isapprox(bestsubset(solver,sparse(X),y,2,500,true), [106.25,53.7799,0.0,0.0], atol=1.0)
        end
    end

    @testset "Test setsolver" begin
        m = Model()
        @variable(m, x[1:5])
        @constraint(m, con[i=1:5], x[6-i] == i)

        for solver in lp_solvers
            setsolver(m, solver)
            @test m.solver == solver
            @test (m.internalModel == nothing) == true
            @test solve(m) == :Optimal
            @test m.solver == solver
            @test (m.internalModel == nothing) == false
        end
    end

    @testset "[model] Setting solve hook" begin
        m = Model()
        @variable(m, x ≥ 0)
        dummy = [1]
        kwarglist = Any[]
        function solvehook(m::Model; kwargs...)
            dummy[1] += 1
            append!(kwarglist, kwargs)
            :DidntDoAnything
        end
        JuMP.setsolvehook(m, solvehook)
        solve(m)
        @test dummy == [2]
        if VERSION >= v"0.7-"
            @test kwarglist == Any[:suppress_warnings => false]
        else
            @test kwarglist == Any[(:suppress_warnings,false)]
        end
    end

    @testset "Setting print hook" begin
        m = Model()
        @variable(m, x ≥ 0)
        dummy = [1]
        function printhook(io::IO, m::Model)
            dummy[1] += 1
        end
        JuMP.setprinthook(m, printhook)
        print(m)
        @test dummy == [2]
    end

    @testset "Test linearindex" begin
        m = Model()
        @variable(m, x[1:5])
        for i in 1:5
            @test isequal(Variable(m,linearindex(x[i])),x[i])
        end
    end

    @testset "Test LinearConstraint from ConstraintRef" begin
        m = Model()
        @variable(m, x)
        @constraint(m, constr, x == 1)
        real_constr = LinearConstraint(constr)
        @test real_constr.terms == @LinearConstraint(x == 1).terms
    end

    @testset "Test getvalue on OneIndexedArrays" begin
        m = Model()
        @variable(m, x[i=1:5], start=i)
        @test typeof(x) == Vector{Variable}
        @test typeof(getvalue(x)) == Vector{Float64}
    end

    @testset "Relaxation keyword argument to solve with $solver" for solver in ip_dual_solvers
        m = Model(solver=solver)
        @variable(m, 1.5 <= y <= 2, Int)
        @variable(m, z, Bin)
        @variable(m, 0.5 <= w <= 1.5, Int)
        @variable(m, 1 <= v <= 2)

        @objective(m, Min, y + z + w + v)

        @test solve(m, relaxation=true) == :Optimal
        @test isapprox(getvalue(y), 1.5, atol=TOL)
        @test isapprox(getvalue(z), 0, atol=TOL)
        @test isapprox(getvalue(w), 0.5, atol=TOL)
        @test isapprox(getvalue(v), 1, atol=TOL)
        @test isapprox(getdual(y), 1, atol=TOL)
        @test isapprox(getdual(z), 1, atol=TOL)
        @test isapprox(getdual(w), 1, atol=TOL)
        @test isapprox(getdual(v), 1, atol=TOL)
        @test isapprox(getobjectivevalue(m), 1.5 + 0 + 0.5 + 1, atol=TOL)

        @test solve(m) == :Optimal
        @test getvalue(y) == 2
        @test getvalue(z) == 0
        @test getvalue(w) == 1
        @test getvalue(v) == 1
        @test getobjectivevalue(m) == 2 + 0 + 1 + 1
    end

    @testset "Relaxation keyword argument to solve (w/ SOS constraints) with $solver" for solver in sos_solvers
        m = Model(solver=solver)
        @variable(m, 1.5 <= y <= 2, Int)
        @variable(m, z, Bin)
        @variable(m, 0.5 <= w <= 1.5, Int)
        @variable(m, 1 <= v <= 2)

        @objective(m, Min, y + z + w + v)
        @variable(m, 1 <= x <= 2, SemiCont)
        @variable(m, -2 <= t <= -1, SemiInt)

        addSOS1(m, [x, 2y, 3z, 4w, 5v, 6t])
        @objective(m, Min, x + y + z + w + v - t)

        @test solve(m, relaxation=true) == :Optimal

        @test getvalue(x) == 0
        @test getvalue(y) == 1.5
        @test getvalue(z) == 0
        @test getvalue(w) == 0.5
        @test getvalue(v) == 1
        @test getvalue(t) == 0
        @test getobjectivevalue(m) == 0 + 1.5 + 0 + 0.5 + 1 + 0
    end

    @testset "Unrecognized keyword argument to solve" begin
        m = Model()
        @test_throws ErrorException solve(m, this_should_throw=true)
    end

    @testset "Solve MIP relaxation with continuous solvers" begin
    @testset "with $solver" for solver in lp_solvers
        m = Model(solver=solver)
        @variable(m, x, Bin)
        @constraint(m, x >= 0.5)
        @objective(m, Min, x)
        @test solve(m, relaxation=true) == :Optimal
        @test isapprox(getvalue(x), 0.5, atol=TOL)
    end
    @testset "with $solver" for solver in nlp_solvers
        m = Model(solver=solver)
        @variable(m, x, Bin)
        @objective(m, Min, x)
        @NLconstraint(m, x >= 0.5)
        @test solve(m, relaxation=true) == :Optimal
        @test isapprox(getvalue(x), 0.5, atol=TOL)
    end
    end

    @testset "Nonliteral exponents in @constraint" begin
        m = Model()
        @variable(m, x)
        foo() = 2
        @constraint(m, x^(foo()) + x^(foo()-1) + x^(foo()-2) == 1)
        @constraint(m, (x-1)^(foo()) + (x-1)^2 + (x-1)^1 + (x-1)^0 == 1)
        @constraint(m, sum(x for i in 1:3)^(foo()) == 1)
        @constraint(m, sum(x for i in 1:3)^(foo()-1) == 1)
        @test m.quadconstr[1].terms == x^2 + x
        @test m.quadconstr[2].terms == x^2 + x^2 - x - x - x - x + x + 1
        @test m.quadconstr[3].terms == x^2 + x^2 + x^2 + x^2 + x^2 + x^2 + x^2 + x^2 + x^2 - 1
        @test m.quadconstr[4].terms == convert(QuadExpr, x + x + x - 1)
    end

    if length(ip_solvers) > 0
        @testset "sets used as indexsets in JuMPArray" begin
            set = BitSet()
            for i in 4:5
                push!(set, i)
            end
            set2 = BitSet()
            for i in 21:23
                push!(set2, i)
            end
            m = Model(solver=ip_solvers[1])
            @variable(m, x[set, set2], Bin)
            @objective(m , Max, sum(sum(x[e,p] for e in set) for p in set2))
            solve(m)
            sol = getvalue(x)
            checked_objval = 0
            for i in keys(sol)
                checked_objval += sol[i...]
            end
            @test checked_objval == 6
        end
    end

    @testset "[model] .^ broadcasting" begin
        m = Model()
        @variable(m, x[1:2])
        @test (x.^2)[1] == x[1]^2
    end

    @testset "Quadratic constraints with zero coefficients with $solver" for solver in quad_solvers
        m = Model(solver=solver)
        @variable(m, 0 <= v <= 2)
        @variable(m, 1 <= x <= 5)
        @constraint(m, v >= 0.0 * x^2 + x)
        @objective(m, Min, v)
        @test solve(m) == :Optimal
        @test isapprox(getvalue(v), 1.0, atol=TOL)
    end

    @testset "Clear duals after MIP -> LP -> MIP with $solver" for solver in ip_dual_solvers
        m = Model(solver=solver)
        @variable(m, x, Bin)
        @objective(m, Max, x)

        @test solve(m) == :Optimal
        @test isnan(getdual(x))

        @test solve(m, relaxation=true) == :Optimal
        @test isapprox(getdual(x), 1.0, atol=TOL)

        @test solve(m) == :Optimal
        @test isnan(getdual(x))
    end

    @testset "Constraints with non-Array AbstractArrays" begin
        m = Model()
        v = @variable(m, [1:3])
        for x in (OffsetArray(v, -length(v)), view(v, :), sparse(v))
            # Version of diagm that works for OffsetArrays:
            A = similar(x, typeof(zero(eltype(x))), (eachindex(x), eachindex(x)))
            for i in eachindex(x), j in eachindex(x)
                A[i, j] = ifelse(i == j, x[i], zero(eltype(x)))
            end

            # No tests, just to make sure that there are no MethodErrors.
            @constraint(m, x + first(x) .== 0)
            @constraint(m, x - first(x) .== 0)
            @constraint(m, (x + 1) + first(x) .== 0)
            @constraint(m, (x + 1) - first(x) .== 0)
            @constraint(m, -x .<= 0)
            @constraint(m, +x .<= 0)
            @SDconstraint(m, A >= 0)

            @test_throws ErrorException @objective(m, Min, x) # vector objective
        end
    end

    @testset "Variable not owned by the conic model" begin
        if !isempty(soc_solvers)
            m = Model()
            M = Model(solver=first(soc_solvers))
            @variable(m, x[1:3])
            @constraint(M, norm(x[1:2]) <= x[3])
            @test_throws JuMP.VariableNotOwnedError solve(M)
        end
    end

   @testset "Variable not owned by the linear model" begin
       if !isempty(lp_solvers)
           m = Model()
           M = Model(solver=first(lp_solvers))
           @variable(m, x)
           @constraint(M, x >= 0)
           @test_throws JuMP.VariableNotOwnedError solve(M)
       end
   end

    @testset "Inconsistency between getobjbound() and getobjectivevalue() #1021" begin
        if !isempty(ip_solvers)
            m = Model(solver=first(ip_solvers))
            @variable(m, 0 <= x <= 5)
            @variable(m, y, Bin)
            @constraint(m, x >= y)
            @objective(m, Min, x - 5)

            @test_throws ErrorException getobjbound(m)

            solve(m)

            @test getobjectivevalue(m) == -5
            @test getobjbound(m) == 0
            @test getobjectivebound(m) == -5

            @objective(m, Min, x - 3)

            @test getobjectivevalue(m) == -5
            @test getobjbound(m) == 0
            @test getobjectivebound(m) == -5
        end
    end

    @testset "getindex for variables and constraints" begin
        m = Model()
        @variable(m, x)
        @test m[:x] == x
        @test_throws KeyError m[:y]
        @constraint(m, c, x <= 1)
        @test m[:c] == c
        @variable(m, c)
        @test_throws Exception m[:c]
    end

    @testset "setindex! for variables and constraints" begin
        m = Model()
        m[:x] = @variable(m)
        @test isa(m[:x], Variable)
        m[:c] = @constraint(m, m[:x] <= 1)
        @test isa(m[:c], JuMP.ConstraintRef)
        m[:c] = @variable(m) # user purposely changes object in m[:c]
        @test isa(m[:c], JuMP.Variable)
    end
end
