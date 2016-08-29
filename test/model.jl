#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/model.jl
# Testing Model printing, basic solving
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using JuMP, FactCheck

# If solvers not loaded, load them (i.e running just these tests)
!isdefined(:lp_solvers) && include("solvers.jl")

# To ensure the tests work on Windows and Linux/OSX, we need
# to use the correct comparison operators
const leq = JuMP.repl[:leq]
const geq = JuMP.repl[:geq]
const  eq = JuMP.repl[:eq]

const TOL = 1e-4

modPath = joinpath(dirname(@__FILE__), "mod")

facts("[model] Check error cases") do
    @fact_throws Model(solver=:Foo)
    modErr = Model()
    @fact_throws setobjectivesense(modErr, :Maximum)
    @variable(modErr, errVar)
    @fact getvalue(errVar) --> isnan
    @fact_throws getdual(errVar)

    modErr = Model()
    @variable(modErr, x, Bin)
    @objective(modErr, Max, x)
    con = @constraint(modErr, x <= 0.5)
    solve(modErr)
    @fact_throws getdual(con)

    modErr = Model()
    @variable(modErr, 0 <= x <= 1)
    @objective(modErr, Max, x)
    @constraint(modErr, x <= -1)
    solve(modErr, suppress_warnings=true)
    @fact getvalue(x) --> isnan
end

facts("[model] Performance warnings") do
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

facts("[model] Test printing a model") do
    modA = Model()
    x = Variable(modA, 0, Inf, :Cont)
    @variable(modA, y <= 5, Int)
    @variable(modA, 2 <= z <= 4)
    @variable(modA, 0 <= r[i=3:6] <= i)
    @objective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @constraintref constraints[1:3]
    constraints[1] = @constraint(modA, 2 <= x+y <= 4)
    constraints[2] = @constraint(modA, sum{r[i],i=3:5} <= (2 - x)/2.0)
    constraints[3] = @constraint(modA, 6y + y <= z + r[6]/1.9)
    #####################################################################
    # Test LP writer
    writeLP(modA, modPath * "A.lp")
    modALP = if VERSION >= v"0.5.0-dev+1866" # leading zero, base Julia PR #14377
        String[
        "Maximize",
        "obj: 0.16666666666666666 VAR1 + 0.16666666666666666 VAR2 + 1 VAR3 + 1 VAR4",
        "Subject To",
        "c1: 1 VAR1 + 1 VAR2 >= 2",
        "c2: 1 VAR1 + 1 VAR2 <= 4",
        "c3: 1 VAR4 + 1 VAR5 + 1 VAR6 + 0.5 VAR1 <= 1",
        "c4: 6 VAR2 + 1 VAR2 - 1 VAR3 - 0.5263157894736842 VAR7 <= 0",
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
    else
        String[
        "Maximize",
        "obj: .16666666666666666 VAR1 + .16666666666666666 VAR2 + 1 VAR3 + 1 VAR4",
        "Subject To",
        "c1: 1 VAR1 + 1 VAR2 >= 2",
        "c2: 1 VAR1 + 1 VAR2 <= 4",
        "c3: 1 VAR4 + 1 VAR5 + 1 VAR6 + .5 VAR1 <= 1",
        "c4: 6 VAR2 + 1 VAR2 - 1 VAR3 - .5263157894736842 VAR7 <= 0",
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
    end
    modAfp = open(modPath * "A.lp")
    lineInd = 1
    while !eof(modAfp)
        line = readline(modAfp)
        @fact strip(line) --> strip(modALP[lineInd])
        lineInd += 1
    end
    close(modAfp)
    #####################################################################
    # Test MPS writer
    writeMPS(modA, modPath * "A.mps")
    if VERSION >= v"0.5.0-dev+1866" # leading zero, base Julia PR #14377
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
    else
        modAMPS = String[
        "NAME   JuMPModel",
        "ROWS",
        " N  OBJ",
        " E  CON1",
        " L  CON2",
        " L  CON3",
        "COLUMNS",
        "    VAR1  CON1  1",
        "    VAR1  CON2  .5",
        "    VAR1  OBJ  -.16666666666666666",
        "    MARKER    'MARKER'                 'INTORG'",
        "    VAR2  CON1  1",
        "    VAR2  CON3  7",
        "    VAR2  OBJ  -.16666666666666666",
        "    MARKER    'MARKER'                 'INTEND'",
        "    VAR3  CON3  -1",
        "    VAR3  OBJ  -1",
        "    VAR4  CON2  1",
        "    VAR4  OBJ  -1",
        "    VAR5  CON2  1",
        "    VAR5  OBJ  0",
        "    VAR6  CON2  1",
        "    VAR6  OBJ  0",
        "    VAR7  CON3  -.5263157894736842",
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
    end
    modAfp = open(modPath * "A.mps")
    lineInd = 1
    while !eof(modAfp)
        line = readline(modAfp)
        @fact chomp(line) --> modAMPS[lineInd]
        lineInd += 1
    end
    close(modAfp)

    # Getter/setters
    @fact MathProgBase.numvar(modA) --> 7
    @fact MathProgBase.numlinconstr(modA) --> 3
    @fact MathProgBase.numquadconstr(modA) --> 0
    @fact MathProgBase.numconstr(modA) --> 3
    @fact getobjectivesense(modA) --> :Max
    setobjectivesense(modA, :Min)
    @fact getobjectivesense(modA) --> :Min
end

facts("[model] Quadratic MPS writer") do
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
        @fact chomp(line) --> modQMPS[lineInd]
        lineInd += 1
    end
    close(modQfp)

    # Getter/setters
    @fact MathProgBase.numvar(modQ) --> 2
    @fact MathProgBase.numlinconstr(modQ) --> 0
    @fact MathProgBase.numquadconstr(modQ) --> 0
    @fact MathProgBase.numconstr(modQ) --> 0
    @fact getobjectivesense(modQ) --> :Min
end


facts("[model] Test solving a MILP") do
for solver in ip_solvers
context("With solver $(typeof(solver))") do
    modA = Model(solver=solver)
    @variable(modA, x >= 0)
    @variable(modA, y <= 5, Int)
    @variable(modA, 2 <= z <= 4)
    @variable(modA, 0 <= r[i=3:6] <= i)
    @objective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @constraint(modA, 2 <= x+y)
    @constraint(modA, sum{r[i],i=3:5} <= (2 - x)/2.0)
    @constraint(modA, 7.0*y <= z + r[6]/1.9)

    @fact solve(modA)       --> :Optimal
    @fact modA.objVal       --> roughly(1.0+4.833334, TOL)
    @fact getvalue(x)       --> roughly(1.0, TOL)
    @fact getvalue(y)       --> roughly(1.0, TOL)
    @fact getvalue(z)       --> roughly(4.0, TOL)
    @fact getvalue(r)[3]    --> roughly(0.5, TOL)
    @fact getvalue(r)[4]    --> roughly(0.0, TOL)
    @fact getvalue(r)[5]    --> roughly(0.0, TOL)
    @fact getvalue(r)[6]    --> greater_than_or_equal(5.7-TOL)
    @fact getvalue(r)[6]    --> less_than_or_equal(6.0+TOL)
    @fact getobjective(modA).aff --> ((x + y)/2.0 + 3.0)/3.0 + z + r[3]
end # solver context
end # loop over solvers
end # facts block

facts("[model] Test solving an LP (Min)") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    modA = Model(solver=solver)
    @variable(modA, x >= 0)
    @variable(modA, y <= 5)
    @variable(modA, 2 <= z <= 4)
    @variable(modA, 0 <= r[i=3:6] <= i)
    @objective(modA, Min, -((x + y)/2.0 + 3.0)/3.0 - z - r[3])
    @constraintref cons[1:3]
    cons[1] = @constraint(modA, x+y >= 2)
    cons[2] = @constraint(modA, sum{r[i],i=3:5} <= (2 - x)/2.0)
    cons[3] = @constraint(modA, 7.0*y <= z + r[6]/1.9)

    # Solution
    @fact solve(modA) --> :Optimal
    @fact getobjectivevalue(modA) --> roughly(-5.8446115, TOL)
    @fact getvalue(x)       --> roughly(0.9774436, TOL)
    @fact getvalue(y)       --> roughly(1.0225563, TOL)
    @fact getvalue(z)       --> roughly(4.0, TOL)
    @fact getvalue(r)[3]    --> roughly(0.5112781, TOL)
    @fact getvalue(r)[4]    --> roughly(0.0, TOL)
    @fact getvalue(r)[5]    --> roughly(0.0, TOL)
    @fact getvalue(r)[6]    --> roughly(6.0, TOL)

    # Reduced costs
    @fact getdual(x)    --> roughly( 0.0, TOL)
    @fact getdual(y)    --> roughly( 0.0, TOL)
    @fact getdual(z)    --> roughly(-1.0714286, TOL)
    @fact getdual(r)[3] --> roughly( 0.0, TOL)
    @fact getdual(r)[4] --> roughly(1.0, TOL)
    @fact getdual(r)[5] --> roughly(1.0, TOL)
    @fact getdual(r)[6] --> roughly(-0.03759398, TOL)

    # Row duals
    @fact getdual(cons)[1] --> roughly( 0.333333, TOL)
    @fact getdual(cons)[2] --> roughly(-1.0, TOL)
    @fact getdual(cons)[3] --> roughly(-0.0714286, TOL)
end # solver context
end # loop over solvers
end # facts block

facts("[model] Test solving an LP (Max)") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    modA = Model(solver=solver)
    @variable(modA, x >= 0)
    @variable(modA, y <= 5)
    @variable(modA, 2 <= z <= 4)
    @variable(modA, 0 <= r[i=3:6] <= i)
    @objective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @constraintref cons[1:3]
    cons[1] = @constraint(modA, x+y >= 2)
    cons[2] = @constraint(modA, sum{r[i],i=3:5} <= (2 - x)/2.0)
    cons[3] = @constraint(modA, 7.0*y <= z + r[6]/1.9)

    # Solution
    @fact solve(modA) --> :Optimal
    @fact getobjectivevalue(modA) --> roughly(5.8446115, TOL)
    @fact getvalue(x)       --> roughly(0.9774436, TOL)
    @fact getvalue(y)       --> roughly(1.0225563, TOL)
    @fact getvalue(z)       --> roughly(4.0, TOL)
    @fact getvalue(r)[3]    --> roughly(0.5112781, TOL)
    @fact getvalue(r)[4]    --> roughly(0.0, TOL)
    @fact getvalue(r)[5]    --> roughly(0.0, TOL)
    @fact getvalue(r)[6]    --> roughly(6.0, TOL)

    # Reduced costs
    @fact getdual(x)    --> roughly( 0.0, TOL)
    @fact getdual(y)    --> roughly( 0.0, TOL)
    @fact getdual(z)    --> roughly( 1.0714286, TOL)
    @fact getdual(r)[3] --> roughly( 0.0, TOL)
    @fact getdual(r)[4] --> roughly(-1.0, TOL)
    @fact getdual(r)[5] --> roughly(-1.0, TOL)
    @fact getdual(r)[6] --> roughly( 0.03759398, TOL)

    # Row duals
    @fact getdual(cons)[1] --> roughly(-0.333333, TOL)
    @fact getdual(cons)[2] --> roughly( 1.0, TOL)
    @fact getdual(cons)[3] --> roughly( 0.0714286, TOL)
end # solver context
end # loop over solvers
end # facts block



facts("[model] Test binary variable handling") do
for solver in ip_solvers
context("With solver $(typeof(solver))") do
    modB = Model(solver=solver)
    @variable(modB, x, Bin)
    @objective(modB, Max, x)
    @constraint(modB, x <= 10)
    status = solve(modB)
    @fact status --> :Optimal
    @fact getvalue(x) --> roughly(1.0, TOL)
end
end
end


facts("[model] Test model copying") do
    source = Model()
    @variable(source, 2 <= x <= 5)
    @variable(source, 0 <= y <= 1, Int)
    @objective(source, Max, 3*x + 1*y)
    @constraint(source, x + 2.0*y <= 6)
    @constraint(source, x*x <= 1)
    addSOS2(source, [x, 2y])
    @SDconstraint(source, x*ones(3,3) >= eye(3,3))
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
    #     @fact source.name == dest.name --> true
    # end

    # Obj
    @objective(source, Max, 1x)
    @fact length(source.obj.aff.coeffs) --> 1
    @fact length(dest.obj.aff.coeffs) --> 2
    @objective(dest, Max, 1x)
    @fact length(source.obj.aff.coeffs) --> 1
    @fact length(dest.obj.aff.coeffs) --> 1
    @objective(dest, Max, 3*x + 1*y)
    @fact length(source.obj.aff.coeffs) --> 1
    @fact length(dest.obj.aff.coeffs) --> 2

    # Constraints
    source.linconstr[1].ub = 5.0
    @fact dest.linconstr[1].ub --> 6.0
    source.quadconstr[1].sense = :(>=)
    @fact dest.quadconstr[1].sense --> :(<=)
    source.sosconstr[1].terms[1] = Variable(source, 2)
    source.sosconstr[1].terms[2] = Variable(source, 1)
    source.sosconstr[1].weights = [2.0,1.0]
    source.sosconstr[1].sostype = :SOS1
    @fact dest.sosconstr[1].terms[1].col --> 1
    @fact dest.sosconstr[1].terms[2].col --> 2
    @fact dest.sosconstr[1].weights --> [1.0,2.0]
    @fact dest.sosconstr[1].sostype --> :SOS2
    @fact length(dest.sdpconstr) --> 2
    xx = copy(x, dest)
    @fact all(t -> isequal(t[1],t[2]), zip(dest.sdpconstr[1].terms, xx*ones(3,3) - eye(3,3))) --> true
    @fact all(t -> isequal(t[1],t[2]), zip(dest.sdpconstr[2].terms, convert(Matrix{AffExpr}, -ones(3,3)))) --> true

    @fact dest.solvehook(dest) --> :Optimal

    @fact Set(collect(keys(dest.varDict))) --> Set([:x,:y,:z,:w,:v])
    @fact isequal(dest.varDict[:x], Variable(dest, 1)) --> true
    @fact isequal(dest.varDict[:y], Variable(dest, 2)) --> true
    @fact all(t -> isequal(t[1], t[2]), zip(dest.varDict[:z], [Variable(dest, 3), Variable(dest, 4), Variable(dest, 5)])) --> true
    @fact all(t -> isequal(t[1], t[2]), zip(dest.varDict[:w].innerArray, [Variable(dest, 6), Variable(dest, 7), Variable(dest, 8)])) --> true
    td = dest.varDict[:v].tupledict
    @fact length(td) --> 2
    @fact isequal(td[:red,1], Variable(dest, 9))  --> true
    @fact isequal(td[:red,3], Variable(dest, 10)) --> true

    # Issue #358
    @fact typeof(dest.linconstr)  --> Array{JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}},1}
    @fact typeof(dest.quadconstr) --> Array{JuMP.GenericQuadConstraint{JuMP.GenericQuadExpr{Float64,JuMP.Variable}},1}

    JuMP.setprinthook(source, m -> 2)
    dest2 = copy(source)
    @fact dest2.printhook(dest2) --> 2

    # Test copying model with multiple variables of same name
    @variable(source, x)
    dest3 = copy(source)
    @fact_throws getvariable(dest3, :x)
    @fact dest3.varDict[:x] --> nothing

    addlazycallback(source, cb -> 3)
    @fact_throws copy(source)
end

type TemporaryExtensionTestType
    x::Int
end
facts("[model] Test extension copy") do
    source = Model()
    source.ext[:extensiontype] = TemporaryExtensionTestType(1)
    @fact_throws copy(source)
    source.ext[:extensiontype] = 1
    dest = copy(source)
    source.ext[:extensiontype] = 2
    @fact haskey(dest.ext, :extensiontype) --> true
    @fact dest.ext[:extensiontype] --> 1
end


facts("[model] Test variable/model 'hygiene'") do
context("Linear constraint") do
    modA = Model(); @variable(modA, x)
    modB = Model(); @variable(modB, y)
    @constraint(modA, x+y == 1)
    @fact_throws solve(modA)
end
context("Linear objective") do
    modA = Model(); @variable(modA, 0 <= x <= 1)
    modB = Model(); @variable(modB, 0 <= y <= 1)
    @objective(modA, Min, x + y)
    @fact_throws solve(modA)
end
context("Quadratic constraint") do
    modA = Model(); @variable(modA, x)
    modB = Model(); @variable(modB, y)
    @constraint(modB, x*y >= 1)
    @fact_throws solve(modB)
end
context("Affine in quadratic constraint") do
    modA = Model(); @variable(modA, x)
    modB = Model(); @variable(modB, y)
    @constraint(modB, y*y + x + y <= 1)
    @fact_throws solve(modB)
end
context("Quadratic objective") do
    modA = Model(); @variable(modA, x)
    modB = Model(); @variable(modB, y)
    @objective(modB, Min, x*y)
    @fact_throws solve(modB)
end
end

facts("[model] Test NaN checking") do
    mod = Model()
    @variable(mod, x)
    @objective(mod, Min, NaN*x)
    @fact_throws solve(mod)
    @objective(mod, Min, NaN*x^2)
    @fact_throws solve(mod)
    @objective(mod, Min, x)
    @constraint(mod, Min, NaN*x == 0)
    @fact_throws solve(mod)
end

facts("[model] Test column-wise modeling") do
    mod = Model()
    @variable(mod, 0 <= x <= 1)
    @variable(mod, 0 <= y <= 1)
    @objective(mod, Max, 5x + 1y)
    @constraint(mod, con[i=1:2], i*x + y <= i+5)
    @variable(mod, 0 <= z1 <= 1, objective=10.0, inconstraints=con, coefficients=[1.0,-2.0])
    @variable(mod, 0 <= z2 <= 1, objective=10.0, inconstraints=Any[con[i] for i in 1:2], coefficients=[1.0,-2.0])
    @fact solve(mod) --> :Optimal
    @fact getvalue(z1) --> roughly(1.0, TOL)
    @fact getvalue(z2) --> roughly(1.0, TOL)

    # do a vectorized version as well
    mod = Model()
    @variable(mod, 0 <= x <= 1)
    @variable(mod, 0 <= y <= 1)
    obj = [5,1]'*[x,y]
    @objective(mod, Max, obj[1])
    A = [1 1
         2 1]
    @constraint(mod, A*[x,y] .<= [6,7])
    @variable(mod, 0 <= z1 <= 1, objective=10.0, inconstraints=con, coefficients=[1.0,-2.0])
    @variable(mod, 0 <= z2 <= 1, objective=10.0, inconstraints=Any[con[i] for i in 1:2], coefficients=[1.0,-2.0])
    @fact solve(mod) --> :Optimal
    @fact getvalue(z1) --> roughly(1.0, TOL)
    @fact getvalue(z2) --> roughly(1.0, TOL)

    # vectorized with sparse matrices
    mod = Model()
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
    @fact solve(mod) --> :Optimal
    @fact getvalue(z1) --> roughly(1.0, TOL)
    @fact getvalue(z2) --> roughly(1.0, TOL)
end

facts("[model] Test all MPS paths") do
    mod = Model()
    @variable(mod, free_var)
    @variable(mod, int_var, Bin)
    @variable(mod, low_var >= 5)
    @constraint(mod, free_var == int_var)
    @constraint(mod, free_var - int_var >= 0)
    @objective(mod, Max, free_var*int_var + low_var)
    writeMPS(mod,"test.mps")
end

facts("[model] Test all LP paths") do
    mod = Model()
    @variable(mod, free_var)
    @objective(mod, Max, free_var*free_var)
    @fact_throws writeLP(mod,"test.lp")
    @objective(mod, Max, free_var)
    @constraint(mod, free_var - 2*free_var == 0)
    @constraint(mod, free_var + 2*free_var >= 1)
    writeLP(mod,"test.lp")
end

facts("[model] Test semi-continuous variables") do
for solver in semi_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @variable(mod, x >= 3, SemiCont)
    @variable(mod, y >= 2, SemiCont)
    @constraint(mod, x + y >= 1)
    @objective(mod, Min, x+y)
    solve(mod)
    @fact getvalue(x) --> roughly(0.0, TOL)
    @fact getvalue(y) --> roughly(2.0, TOL)
end; end; end

facts("[model] Test semi-integer variables") do
for solver in semi_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @variable(mod, x >= 3, SemiInt)
    @variable(mod, y >= 2, SemiInt)
    @constraint(mod, x + y >= 2.5)
    @objective(mod, Min, x+1.1y)
    solve(mod)
    @fact getvalue(x) --> roughly(3.0, TOL)
    @fact getvalue(y) --> 0.0
end; end; end

facts("[model] Test fixed variables don't leak through MPB") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @variable(mod, 0 <= x[1:3] <= 2)
    @variable(mod, y[k=1:2] == k)
    @objective(mod, Min, x[1] + x[2] + x[3] + y[1] + y[2])
    solve(mod)
    for i in 1:3
        @fact getvalue(x[i]) --> roughly(0, TOL)
    end
    for k in 1:2
        @fact getvalue(y[k]) --> roughly(k, TOL)
    end
end; end
for solver in ip_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @variable(mod, x[1:3], Bin)
    @variable(mod, y[k=1:2] == k)
    JuMP.build(mod)
    @fact MathProgBase.getvartype(internalmodel(mod)) --> [:Bin,:Bin,:Bin,:Cont,:Cont]
end; end; end


facts("[model] Test SOS constraints") do
for solver in sos_solvers
context("With solver $(typeof(solver))") do
    modS = Model(solver=solver)

    @variable(modS, x[1:3], Bin)
    @variable(modS, y[1:5], Bin)
    @variable(modS, z)
    @variable(modS, w)

    @objective(modS, Max, z+w)

    a = [1,2,3]
    b = [5,4,7,2,1]

    @constraint(modS, z == sum{a[i]*x[i], i=1:3})
    @constraint(modS, w == sum{b[i]*y[i], i=1:5})

    @fact_throws constructSOS([x[1]+y[1]])
    @fact_throws constructSOS([1z])

    addSOS1(modS, [a[i]x[i] for i in 1:3])
    addSOS2(modS, [b[i]y[i] for i in 1:5])

    @fact_throws addSOS1(modS, [x[1], x[1]+x[2]])

    @fact solve(modS) --> :Optimal
    @fact modS.objVal --> roughly(15.0, TOL)
    @fact getvalue(z) --> roughly( 3.0, TOL)
    @fact getvalue(w) --> roughly(12.0, TOL)


    m = Model(solver=solver)
    ub = [1,1,2]
    @variable(m, x[i=1:3] <= ub[i])
    @objective(m, Max, 2x[1]+x[2]+x[3])
    addSOS1(m, [x[1],2x[2]])
    addSOS1(m, [x[1],2x[3]])

    @fact solve(m) --> :Optimal
    @fact getvalue(x)[:] --> roughly([0.0,1.0,2.0], TOL)
    @fact getobjectivevalue(m) --> roughly(3.0, TOL)
end; end; end

facts("[model] Test vectorized model creation") do

    A = sprand(50,10,0.15)
    B = sprand(50, 7,0.2)
    modV = Model()
    @variable(modV, x[1:10])
    @variable(modV, y[1:7])
    @constraint(modV, A*x + B*y .<= 1)
    obj = (x'*2A')*(2A*x) + (B*2y)'*(B*(2y))
    @objective(modV, Max, obj[1])

    modS = Model()
    @variable(modS, x[1:10])
    @variable(modS, y[1:7])
    for i in 1:50
        @constraint(modS, sum{A[i,j]*x[j], j=1:10} + sum{B[i,k]*y[k], k=1:7} <= 1)
    end
    AA, BB = 4A'*A, 4B'*A
    @objective(modS, Max, sum{AA[i,j]*x[i]*x[j], i=1:10,j=1:10} + sum{BB[i,j]*y[i]*y[j], i=1:7, j=1:7})

    @fact JuMP.prepConstrMatrix(modV) --> JuMP.prepConstrMatrix(modS)
    @fact JuMP.prepProblemBounds(modV) --> JuMP.prepProblemBounds(modS)
end

facts("[model] Test MIQP vectorization") do
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
        @objective(mod, Min, obj[1])
        @constraint(mod, eye(p)*β .<=  M*eye(p)*z)
        @constraint(mod, eye(p)*β .>= -M*eye(p)*z)
        @constraint(mod, sum(z) == K)
        solve(mod)
        return getvalue(β)[:]
    end
    include(joinpath("data","miqp_vector.jl")) # loads X and q
    y = X * [100, 50, 10, 1] + 20*q
    for solver in quad_solvers
        @fact bestsubset(solver,X,y,2,500,false) --> roughly([101.789,49.414,8.63904,1.72663], 10TOL)
        @fact bestsubset(solver,sparse(X),y,2,500,false) --> roughly([101.789,49.414,8.63904,1.72663], 10TOL)
    end
    for solver in quad_mip_solvers
        y = X * [100, 50, 10, 1] + 20*q
        @fact bestsubset(solver,X,y,2,500,true) --> roughly([106.25,53.7799,0.0,0.0], 1.0)
        @fact bestsubset(solver,sparse(X),y,2,500,true) --> roughly([106.25,53.7799,0.0,0.0], 1.0)
    end
end

facts("[model] Test setsolver") do
    m = Model()
    @variable(m, x[1:5])
    @constraint(m, con[i=1:5], x[6-i] == i)
    @fact solve(m) --> :Optimal

    for solver in lp_solvers
        setsolver(m, solver)
        @fact m.solver --> solver
        @fact (m.internalModel == nothing) --> true
        @fact solve(m) --> :Optimal
        @fact m.solver --> solver
        @fact (m.internalModel == nothing) --> false
    end
end

facts("[model] Setting solve hook") do
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
    @fact dummy --> [2]
    @fact kwarglist --> Any[(:suppress_warnings,false)]
end

facts("[model] Setting print hook") do
    m = Model()
    @variable(m, x ≥ 0)
    dummy = [1]
    function printhook(io::IO, m::Model)
        dummy[1] += 1
    end
    JuMP.setprinthook(m, printhook)
    print(m)
    @fact dummy --> [2]
end

facts("[model] Test linearindex") do
    m = Model()
    @variable(m, x[1:5])
    for i in 1:5
        @fact isequal(Variable(m,linearindex(x[i])),x[i]) --> true
    end
end

facts("[model] Test LinearConstraint from ConstraintRef") do
    m = Model()
    @variable(m, x)
    @constraint(m, constr, x == 1)
    real_constr = LinearConstraint(constr)
    @fact real_constr.terms --> @LinearConstraint(x == 1).terms
end

facts("[model] Test getvalue on OneIndexedArrays") do
    m = Model()
    @variable(m, x[i=1:5], start=i)
    @fact typeof(x) --> Vector{Variable}
    @fact typeof(getvalue(x)) --> Vector{Float64}
end

facts("[model] Relaxation keyword argument to solve") do
    m = Model()
    @variable(m, 1.5 <= y <= 2, Int)
    @variable(m, z, Bin)
    @variable(m, 0.5 <= w <= 1.5, Int)
    @variable(m, 1 <= v <= 2)

    @objective(m, Min, y + z + w + v)

    # Force LP solver since not all MIP solvers
    # return duals (i.e. Cbc)

    for solver in lp_solvers
    context("With solver $(typeof(solver))") do
        setsolver(m, solver)
        @fact solve(m, relaxation=true) --> :Optimal
        @fact getvalue(y) --> roughly(1.5, TOL)
        @fact getvalue(z) --> roughly(0, TOL)
        @fact getvalue(w) --> roughly(0.5, TOL)
        @fact getvalue(v) --> roughly(1, TOL)
        @fact getdual(y) --> roughly(1, TOL)
        @fact getdual(z) --> roughly(1, TOL)
        @fact getdual(w) --> roughly(1, TOL)
        @fact getdual(v) --> roughly(1, TOL)
        @fact getobjectivevalue(m) --> roughly(1.5 + 0 + 0.5 + 1, TOL)
    end; end

    # Let JuMP choose solver again
    setsolver(m, JuMP.UnsetSolver())
    @fact solve(m) --> :Optimal
    @fact getvalue(y) --> 2
    @fact getvalue(z) --> 0
    @fact getvalue(w) --> 1
    @fact getvalue(v) --> 1
    @fact getobjectivevalue(m) --> 2 + 0 + 1 + 1

    @variable(m, 1 <= x <= 2, SemiCont)
    @variable(m, -2 <= t <= -1, SemiInt)

    addSOS1(m, [x, 2y, 3z, 4w, 5v, 6t])
    @objective(m, Min, x + y + z + w + v - t)

    @fact solve(m, relaxation=true) --> :Optimal

    @fact getvalue(x) --> 0
    @fact getvalue(y) --> 1.5
    @fact getvalue(z) --> 0
    @fact getvalue(w) --> 0.5
    @fact getvalue(v) --> 1
    @fact getvalue(t) --> 0
    @fact getobjectivevalue(m) --> 0 + 1.5 + 0 + 0.5 + 1 + 0
end

facts("[model] Unrecognized keyword argument to solve") do
    m = Model()
    @fact_throws solve(m, this_should_throw=true)
end

facts("[model] Solve MIP relaxation with continuous solvers") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, x, Bin)
    @constraint(m, x >= 0.5)
    @objective(m, Min, x)
    @fact solve(m, relaxation=true) --> :Optimal
    @fact getvalue(x) --> roughly(0.5, TOL)
end; end
for solver in nlp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @variable(m, x, Bin)
    @objective(m, Min, x)
    @NLconstraint(m, x >= 0.5)
    @fact solve(m, relaxation=true) --> :Optimal
    @fact getvalue(x) --> roughly(0.5, TOL)
end; end
end

facts("[model] Nonliteral exponents in @constraint") do
    m = Model()
    @variable(m, x)
    foo() = 2
    @constraint(m, x^(foo()) + x^(foo()-1) + x^(foo()-2) == 1)
    @constraint(m, (x-1)^(foo()) + (x-1)^2 + (x-1)^1 + (x-1)^0 == 1)
    @constraint(m, sum{x, i in 1:3}^(foo()) == 1)
    @constraint(m, sum{x, i in 1:3}^(foo()-1) == 1)
    @fact m.quadconstr[1].terms --> x^2 + x
    @fact m.quadconstr[2].terms --> x^2 + x^2 - x - x - x - x + x + 1
    @fact m.quadconstr[3].terms --> x^2 + x^2 + x^2 + x^2 + x^2 + x^2 + x^2 + x^2 + x^2 - 1
    @fact m.quadconstr[4].terms --> QuadExpr(x + x + x - 1)
end
