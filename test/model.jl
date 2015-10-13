#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
using Compat

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
    @fact_throws setObjectiveSense(modErr, :Maximum)
    @defVar(modErr, errVar)
    @fact getValue(errVar) --> isnan
    @fact_throws getDual(errVar)

    modErr = Model()
    @defVar(modErr, x, Bin)
    @setObjective(modErr, Max, x)
    con = @addConstraint(modErr, x <= 0.5)
    solve(modErr)
    @fact_throws getDual(con)

    modErr = Model()
    @defVar(modErr, 0 <= x <= 1)
    @setObjective(modErr, Max, x)
    @addConstraint(modErr, x <= -1)
    solve(modErr, suppress_warnings=true)
    @fact getValue(x) --> isnan
end

facts("[model] Performance warnings") do
    m = Model()
    @defVar(m, x[1:2], start=0)
    for i in 1:500
       getValue(x)
    end
    q = 0
    for i in 1:30000
       q += 3x[1]
    end
end

facts("[model] Test printing a model") do
    modA = Model()
    @defVar(modA, x >= 0)
    @defVar(modA, y <= 5, Int)
    @defVar(modA, 2 <= z <= 4)
    @defVar(modA, 0 <= r[i=3:6] <= i)
    @setObjective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @defConstrRef constraints[1:3]
    constraints[1] = @addConstraint(modA, 2 <= x+y <= 4)
    constraints[2] = @addConstraint(modA, sum{r[i],i=3:5} <= (2 - x)/2.0)
    constraints[3] = @addConstraint(modA, 6y + y <= z + r[6]/1.9)
    #####################################################################
    # Test LP writer
    writeLP(modA, modPath * "A.lp")
    modALP = ASCIIString[
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
    modAMPS = ASCIIString[
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
        @fact chomp(line) --> modAMPS[lineInd]
        lineInd += 1
    end
    close(modAfp)

    # Getter/setters
    @fact MathProgBase.numvar(modA) --> 7
    @fact MathProgBase.numlinconstr(modA) --> 3
    @fact MathProgBase.numquadconstr(modA) --> 0
    @fact MathProgBase.numconstr(modA) --> 3
    @fact getObjectiveSense(modA) --> :Max
    setObjectiveSense(modA, :Min)
    @fact getObjectiveSense(modA) --> :Min
end

facts("[model] Quadratic MPS writer") do
    modQ = Model()
    @defVar(modQ, x == 1)
    @defVar(modQ, y ≥ 0)
    @setObjective(modQ, Min, x^2 - 2*x*y + y^2 + x)

    #####################################################################
    # Test MPS writer
    writeMPS(modQ, modPath * "Q.mps")
    modQMPS = ASCIIString[
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
    @fact getObjectiveSense(modQ) --> :Min
end


facts("[model] Test solving a MILP") do
for solver in ip_solvers
context("With solver $(typeof(solver))") do
    modA = Model(solver=solver)
    @defVar(modA, x >= 0)
    @defVar(modA, y <= 5, Int)
    @defVar(modA, 2 <= z <= 4)
    @defVar(modA, 0 <= r[i=3:6] <= i)
    @setObjective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @addConstraint(modA, 2 <= x+y)
    @addConstraint(modA, sum{r[i],i=3:5} <= (2 - x)/2.0)
    @addConstraint(modA, 7.0*y <= z + r[6]/1.9)

    @fact solve(modA)       --> :Optimal
    @fact modA.objVal       --> roughly(1.0+4.833334, TOL)
    @fact getValue(x)       --> roughly(1.0, TOL)
    @fact getValue(y)       --> roughly(1.0, TOL)
    @fact getValue(z)       --> roughly(4.0, TOL)
    @fact getValue(r)[3]    --> roughly(0.5, TOL)
    @fact getValue(r)[4]    --> roughly(0.0, TOL)
    @fact getValue(r)[5]    --> roughly(0.0, TOL)
    @fact getValue(r)[6]    --> roughly(6.0, TOL)
end # solver context
end # loop over solvers
end # facts block

facts("[model] Test solving an LP (Min)") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    modA = Model(solver=solver)
    @defVar(modA, x >= 0)
    @defVar(modA, y <= 5)
    @defVar(modA, 2 <= z <= 4)
    @defVar(modA, 0 <= r[i=3:6] <= i)
    @setObjective(modA, Min, -((x + y)/2.0 + 3.0)/3.0 - z - r[3])
    @defConstrRef cons[1:3]
    cons[1] = @addConstraint(modA, x+y >= 2)
    cons[2] = @addConstraint(modA, sum{r[i],i=3:5} <= (2 - x)/2.0)
    cons[3] = @addConstraint(modA, 7.0*y <= z + r[6]/1.9)

    # Solution
    @fact solve(modA) --> :Optimal
    @fact getObjectiveValue(modA) --> roughly(-5.8446115, TOL)
    @fact getValue(x)       --> roughly(0.9774436, TOL)
    @fact getValue(y)       --> roughly(1.0225563, TOL)
    @fact getValue(z)       --> roughly(4.0, TOL)
    @fact getValue(r)[3]    --> roughly(0.5112781, TOL)
    @fact getValue(r)[4]    --> roughly(0.0, TOL)
    @fact getValue(r)[5]    --> roughly(0.0, TOL)
    @fact getValue(r)[6]    --> roughly(6.0, TOL)

    # Reduced costs
    @fact getDual(x)    --> roughly( 0.0, TOL)
    @fact getDual(y)    --> roughly( 0.0, TOL)
    @fact getDual(z)    --> roughly(-1.0714286, TOL)
    @fact getDual(r)[3] --> roughly( 0.0, TOL)
    @fact getDual(r)[4] --> roughly(1.0, TOL)
    @fact getDual(r)[5] --> roughly(1.0, TOL)
    @fact getDual(r)[6] --> roughly(-0.03759398, TOL)

    # Row duals
    @fact getDual(cons)[1] --> roughly( 0.333333, TOL)
    @fact getDual(cons)[2] --> roughly(-1.0, TOL)
    @fact getDual(cons)[3] --> roughly(-0.0714286, TOL)
end # solver context
end # loop over solvers
end # facts block

facts("[model] Test solving an LP (Max)") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    modA = Model(solver=solver)
    @defVar(modA, x >= 0)
    @defVar(modA, y <= 5)
    @defVar(modA, 2 <= z <= 4)
    @defVar(modA, 0 <= r[i=3:6] <= i)
    @setObjective(modA, Max, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
    @defConstrRef cons[1:3]
    cons[1] = @addConstraint(modA, x+y >= 2)
    cons[2] = @addConstraint(modA, sum{r[i],i=3:5} <= (2 - x)/2.0)
    cons[3] = @addConstraint(modA, 7.0*y <= z + r[6]/1.9)

    # Solution
    @fact solve(modA) --> :Optimal
    @fact getObjectiveValue(modA) --> roughly(5.8446115, TOL)
    @fact getValue(x)       --> roughly(0.9774436, TOL)
    @fact getValue(y)       --> roughly(1.0225563, TOL)
    @fact getValue(z)       --> roughly(4.0, TOL)
    @fact getValue(r)[3]    --> roughly(0.5112781, TOL)
    @fact getValue(r)[4]    --> roughly(0.0, TOL)
    @fact getValue(r)[5]    --> roughly(0.0, TOL)
    @fact getValue(r)[6]    --> roughly(6.0, TOL)

    # Reduced costs
    @fact getDual(x)    --> roughly( 0.0, TOL)
    @fact getDual(y)    --> roughly( 0.0, TOL)
    @fact getDual(z)    --> roughly( 1.0714286, TOL)
    @fact getDual(r)[3] --> roughly( 0.0, TOL)
    @fact getDual(r)[4] --> roughly(-1.0, TOL)
    @fact getDual(r)[5] --> roughly(-1.0, TOL)
    @fact getDual(r)[6] --> roughly( 0.03759398, TOL)

    # Row duals
    @fact getDual(cons)[1] --> roughly(-0.333333, TOL)
    @fact getDual(cons)[2] --> roughly( 1.0, TOL)
    @fact getDual(cons)[3] --> roughly( 0.0714286, TOL)
end # solver context
end # loop over solvers
end # facts block



facts("[model] Test binary variable handling") do
for solver in ip_solvers
context("With solver $(typeof(solver))") do
    modB = Model(solver=solver)
    @defVar(modB, x, Bin)
    @setObjective(modB, Max, x)
    @addConstraint(modB, x <= 10)
    status = solve(modB)
    @fact status --> :Optimal
    @fact getValue(x) --> roughly(1.0, TOL)
end
end
end


facts("[model] Test model copying") do
    source = Model()
    @defVar(source, 2 <= x <= 5)
    @defVar(source, 0 <= y <= 1, Int)
    @setObjective(source, Max, 3*x + 1*y)
    @addConstraint(source, x + 2.0*y <= 6)
    @addConstraint(source, x*x <= 1)
    addSOS2(source, [x, 2y])
    @addSDPConstraint(source, x*ones(3,3) >= eye(3,3))
    @addSDPConstraint(source, ones(3,3) <= 0)
    @defVar(source, z[1:3])
    @defVar(source, w[2:4]) # JuMPArray
    @defVar(source, v[[:red],1:3]) # JuMPDict
    setSolveHook(source, m -> 1)

    # uncomment when NLP copying is implemented
    # @addNLConstraint(source, c[k=1:3], x^2 + y^3 * sin(x+k*y) >= 1)

    dest = copy(source)

    # uncomment when NLP copying is implemented
    # for name in source.nlpdata
    #     @fact source.name == dest.name --> true
    # end

    # Obj
    @setObjective(source, Max, 1x)
    @fact length(source.obj.aff.coeffs) --> 1
    @fact length(dest.obj.aff.coeffs) --> 2
    @setObjective(dest, Max, 1x)
    @fact length(source.obj.aff.coeffs) --> 1
    @fact length(dest.obj.aff.coeffs) --> 1
    @setObjective(dest, Max, 3*x + 1*y)
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

    @fact dest.solvehook(dest) --> 1

    @fact Set(collect(keys(dest.varDict))) --> Set([:x,:y,:z,:w,:v])
    @fact isequal(dest.varDict[:x], Variable(dest, 1)) --> true
    @fact isequal(dest.varDict[:y], Variable(dest, 2)) --> true
    @fact all(t -> isequal(t[1], t[2]), zip(dest.varDict[:z], [Variable(dest, 3), Variable(dest, 4), Variable(dest, 5)])) --> true
    @fact all(t -> isequal(t[1], t[2]), zip(dest.varDict[:w].innerArray, [Variable(dest, 6), Variable(dest, 7), Variable(dest, 8)])) --> true
    td = dest.varDict[:v].tupledict
    @fact length(td) --> 3
    @fact isequal(td[:red,1], Variable(dest, 9))  --> true
    @fact isequal(td[:red,2], Variable(dest, 10)) --> true
    @fact isequal(td[:red,3], Variable(dest, 11)) --> true

    # Issue #358
    @fact typeof(dest.linconstr)  --> Array{JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}},1}
    @fact typeof(dest.quadconstr) --> Array{JuMP.GenericQuadConstraint{JuMP.GenericQuadExpr{Float64,JuMP.Variable}},1}

    setPrintHook(source, m -> 2)
    dest2 = copy(source)
    @fact dest2.printhook(dest2) --> 2

    addLazyCallback(source, cb -> 3)
    @fact_throws copy(source)
end


facts("[model] Test variable/model 'hygiene'") do
context("Linear constraint") do
    modA = Model(); @defVar(modA, x)
    modB = Model(); @defVar(modB, y)
    @addConstraint(modA, x+y == 1)
    @fact_throws solve(modA)
end
context("Linear objective") do
    modA = Model(); @defVar(modA, 0 <= x <= 1)
    modB = Model(); @defVar(modB, 0 <= y <= 1)
    @setObjective(modA, Min, x + y)
    @fact_throws solve(modA)
end
context("Quadratic constraint") do
    modA = Model(); @defVar(modA, x)
    modB = Model(); @defVar(modB, y)
    @addConstraint(modB, x*y >= 1)
    @fact_throws solve(modB)
end
context("Affine in quadratic constraint") do
    modA = Model(); @defVar(modA, x)
    modB = Model(); @defVar(modB, y)
    @addConstraint(modB, y*y + x + y <= 1)
    @fact_throws solve(modB)
end
context("Quadratic objective") do
    modA = Model(); @defVar(modA, x)
    modB = Model(); @defVar(modB, y)
    @setObjective(modB, Min, x*y)
    @fact_throws solve(modB)
end
end

facts("[model] Test NaN checking") do
    mod = Model()
    @defVar(mod, x)
    @setObjective(mod, Min, NaN*x)
    @fact_throws solve(mod)
    setObjective(mod, :Min, NaN*x^2)
    @fact_throws solve(mod)
    @setObjective(mod, Min, x)
    @addConstraint(mod, Min, NaN*x == 0)
    @fact_throws solve(mod)
end

facts("[model] Test column-wise modeling") do
    mod = Model()
    @defVar(mod, 0 <= x <= 1)
    @defVar(mod, 0 <= y <= 1)
    @setObjective(mod, Max, 5x + 1y)
    @addConstraint(mod, con[i=1:2], i*x + y <= i+5)
    @defVar(mod, 0 <= z1 <= 1, objective=10.0, inconstraints=con, coefficients=[1.0,-2.0])
    @defVar(mod, 0 <= z2 <= 1, objective=10.0, inconstraints=Any[con[i] for i in 1:2], coefficients=[1.0,-2.0])
    @fact solve(mod) --> :Optimal
    @fact getValue(z1) --> roughly(1.0, TOL)
    @fact getValue(z2) --> roughly(1.0, TOL)

    # do a vectorized version as well
    mod = Model()
    @defVar(mod, 0 <= x <= 1)
    @defVar(mod, 0 <= y <= 1)
    obj = [5,1]'*[x,y]
    @setObjective(mod, Max, obj[1])
    A = [1 1
         2 1]
    @addConstraint(mod, A*[x,y] .<= [6,7])
    @defVar(mod, 0 <= z1 <= 1, objective=10.0, inconstraints=con, coefficients=[1.0,-2.0])
    @defVar(mod, 0 <= z2 <= 1, objective=10.0, inconstraints=Any[con[i] for i in 1:2], coefficients=[1.0,-2.0])
    @fact solve(mod) --> :Optimal
    @fact getValue(z1) --> roughly(1.0, TOL)
    @fact getValue(z2) --> roughly(1.0, TOL)

    # vectorized with sparse matrices
    mod = Model()
    @defVar(mod, 0 <= x <= 1)
    @defVar(mod, 0 <= y <= 1)
    # TODO: get this to work by adding method for Ac_mul_B!
    # obj = sparse([5,1])'*[x,y]
    obj = sparse([5,1]')*[x,y]
    @setObjective(mod, Max, obj[1])
    A = sparse([1 1
                2 1])
    @addConstraint(mod, A*[x,y] .<= [6,7])
    @defVar(mod, 0 <= z1 <= 1, objective=10.0, inconstraints=con, coefficients=[1.0,-2.0])
    @defVar(mod, 0 <= z2 <= 1, objective=10.0, inconstraints=Any[con[i] for i in 1:2], coefficients=[1.0,-2.0])
    @fact solve(mod) --> :Optimal
    @fact getValue(z1) --> roughly(1.0, TOL)
    @fact getValue(z2) --> roughly(1.0, TOL)
end

facts("[model] Test all MPS paths") do
    mod = Model()
    @defVar(mod, free_var)
    @defVar(mod, int_var, Bin)
    @defVar(mod, low_var >= 5)
    @addConstraint(mod, free_var == int_var)
    @addConstraint(mod, free_var - int_var >= 0)
    setObjective(mod, :Max, free_var*int_var + low_var)
    writeMPS(mod,"test.mps")
end

facts("[model] Test all LP paths") do
    mod = Model()
    @defVar(mod, free_var)
    setObjective(mod, :Max, free_var*free_var)
    @fact_throws writeLP(mod,"test.lp")
    @setObjective(mod, Max, free_var)
    @addConstraint(mod, free_var - 2*free_var == 0)
    @addConstraint(mod, free_var + 2*free_var >= 1)
    writeLP(mod,"test.lp")
end

facts("[model] Test semi-continuous variables") do
for solver in semi_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @defVar(mod, x >= 3, SemiCont)
    @defVar(mod, y >= 2, SemiCont)
    @addConstraint(mod, x + y >= 1)
    @setObjective(mod, Min, x+y)
    solve(mod)
    @fact getValue(x) --> roughly(0.0, TOL)
    @fact getValue(y) --> roughly(2.0, TOL)
end; end; end

facts("[model] Test semi-integer variables") do
for solver in semi_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @defVar(mod, x >= 3, SemiInt)
    @defVar(mod, y >= 2, SemiInt)
    @addConstraint(mod, x + y >= 2.5)
    @setObjective(mod, Min, x+1.1y)
    solve(mod)
    @fact getValue(x) --> roughly(3.0, TOL)
    @fact getValue(y) --> 0.0
end; end; end

facts("[model] Test fixed variables don't leak through MPB") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @defVar(mod, 0 <= x[1:3] <= 2)
    @defVar(mod, y[k=1:2] == k)
    @setObjective(mod, Min, x[1] + x[2] + x[3] + y[1] + y[2])
    solve(mod)
    for i in 1:3
        @fact getValue(x[i]) --> roughly(0, TOL)
    end
    for k in 1:2
        @fact getValue(y[k]) --> roughly(k, TOL)
    end
end; end
for solver in ip_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @defVar(mod, x[1:3], Bin)
    @defVar(mod, y[k=1:2] == k)
    buildInternalModel(mod)
    @fact MathProgBase.getvartype(getInternalModel(mod)) --> [:Bin,:Bin,:Bin,:Cont,:Cont]
end; end; end


facts("[model] Test SOS constraints") do
for solver in sos_solvers
context("With solver $(typeof(solver))") do
    modS = Model(solver=solver)

    @defVar(modS, x[1:3], Bin)
    @defVar(modS, y[1:5], Bin)
    @defVar(modS, z)
    @defVar(modS, w)

    setObjective(modS, :Max, z+w)

    a = [1,2,3]
    b = [5,4,7,2,1]

    @addConstraint(modS, z == sum{a[i]*x[i], i=1:3})
    @addConstraint(modS, w == sum{b[i]*y[i], i=1:5})

    @fact_throws constructSOS([x[1]+y[1]])
    @fact_throws constructSOS([1z])

    addSOS1(modS, [a[i]x[i] for i in 1:3])
    addSOS2(modS, [b[i]y[i] for i in 1:5])

    @fact_throws addSOS1(modS, [x[1], x[1]+x[2]])

    @fact solve(modS) --> :Optimal
    @fact modS.objVal --> roughly(15.0, TOL)
    @fact getValue(z) --> roughly( 3.0, TOL)
    @fact getValue(w) --> roughly(12.0, TOL)


    m = Model(solver=solver)
    ub = [1,1,2]
    @defVar(m, x[i=1:3] <= ub[i])
    @setObjective(m, Max, 2x[1]+x[2]+x[3])
    addSOS1(m, [x[1],2x[2]])
    addSOS1(m, [x[1],2x[3]])

    @fact solve(m) --> :Optimal
    @fact getValue(x)[:] --> roughly([0.0,1.0,2.0], TOL)
    @fact getObjectiveValue(m) --> roughly(3.0, TOL)
end; end; end

facts("[model] Test vectorized model creation") do

    A = sprand(50,10,0.15)
    B = sprand(50, 7,0.2)
    modV = Model()
    @defVar(modV, x[1:10])
    @defVar(modV, y[1:7])
    @addConstraint(modV, A*x + B*y .<= 1)
    obj = (x'*2A')*(2A*x) + (B*2y)'*(B*(2y))
    @setObjective(modV, Max, obj[1])

    modS = Model()
    @defVar(modS, x[1:10])
    @defVar(modS, y[1:7])
    for i in 1:50
        @addConstraint(modS, sum{A[i,j]*x[j], j=1:10} + sum{B[i,k]*y[k], k=1:7} <= 1)
    end
    AA, BB = 4A'*A, 4B'*A
    @setObjective(modS, Max, sum{AA[i,j]*x[i]*x[j], i=1:10,j=1:10} + sum{BB[i,j]*y[i]*y[j], i=1:7, j=1:7})

    @fact JuMP.prepConstrMatrix(modV) --> JuMP.prepConstrMatrix(modS)
    @fact JuMP.prepProblemBounds(modV) --> JuMP.prepProblemBounds(modS)
end

facts("[model] Test MIQP vectorization") do
    n = 1000
    p = 4
    function bestsubset(solver,X,y,K,M,integer)
        mod = Model(solver=solver)
        @defVar(mod, β[1:p])
        if integer
            @defVar(mod, z[1:p], Bin)
        else
            @defVar(mod, 0 <= z[1:p] <= 1)
        end
        obj = (y-X*β)'*(y-X*β)
        @setObjective(mod, Min, obj[1])
        @addConstraint(mod, eye(p)*β .<=  M*eye(p)*z)
        @addConstraint(mod, eye(p)*β .>= -M*eye(p)*z)
        @addConstraint(mod, sum(z) == K)
        solve(mod)
        return getValue(β)[:]
    end
    include(joinpath("data","miqp_vector.jl")) # loads X and q
    y = X * [100, 50, 10, 1] + 20*q
    for solver in quad_solvers
        @fact bestsubset(solver,X,y,2,500,false) --> roughly([101.789,49.414,8.63904,1.72663], 10TOL)
        @fact bestsubset(solver,sparse(X),y,2,500,false) --> roughly([101.789,49.414,8.63904,1.72663], 10TOL)
    end
    for solver in quad_mip_solvers
        y = X * [100, 50, 10, 1] + 20*q
        @fact bestsubset(solver,X,y,2,500,true) --> roughly([106.25,53.7799,0.0,0.0], 10TOL)
        @fact bestsubset(solver,sparse(X),y,2,500,true) --> roughly([106.25,53.7799,0.0,0.0], 10TOL)
    end
end

facts("[model] Test setSolver") do
    m = Model()
    @defVar(m, x[1:5])
    @addConstraint(m, con[i=1:5], x[6-i] == i)
    @fact solve(m) --> :Optimal

    for solver in lp_solvers
        setSolver(m, solver)
        @fact m.solver --> solver
        @fact (m.internalModel == nothing) --> true
        @fact solve(m) --> :Optimal
        @fact m.solver --> solver
        @fact (m.internalModel == nothing) --> false
    end
end

facts("[model] Setting solve hook") do
    m = Model()
    @defVar(m, x ≥ 0)
    dummy = [1]
    kwarglist = Any[]
    function solvehook(m::Model; kwargs...)
        dummy[1] += 1
        append!(kwarglist, kwargs)
        nothing
    end
    setSolveHook(m, solvehook)
    solve(m)
    @fact dummy --> [2]
    @fact kwarglist --> Any[(:suppress_warnings,false)]
end

facts("[model] Setting print hook") do
    m = Model()
    @defVar(m, x ≥ 0)
    dummy = [1]
    function printhook(io::IO, m::Model)
        dummy[1] += 1
    end
    setPrintHook(m, printhook)
    print(m)
    @fact dummy --> [2]
end

facts("[model] Test getLinearIndex") do
    m = Model()
    @defVar(m, x[1:5])
    for i in 1:5
        @fact isequal(Variable(m,getLinearIndex(x[i])),x[i]) --> true
    end
end

facts("[model] Test LinearConstraint from ConstraintRef") do
    m = Model()
    @defVar(m, x)
    @addConstraint(m, constr, x == 1)
    real_constr = LinearConstraint(constr)
    @fact real_constr.terms --> @LinearConstraint(x == 1).terms
end

facts("[model] Test getValue on OneIndexedArrays") do
    m = Model()
    @defVar(m, x[i=1:5], start=i)
    @fact typeof(x) --> Vector{Variable}
    @fact typeof(getValue(x)) --> Vector{Float64}
end

facts("[model] Relaxation keyword argument to solve") do
    m = Model()
    @defVar(m, 1.5 <= y <= 2, Int)
    @defVar(m, z, Bin)
    @defVar(m, 0.5 <= w <= 1.5, Int)
    @defVar(m, 1 <= v <= 2)

    @setObjective(m, Min, y + z + w + v)

    @fact solve(m, relaxation=true) --> :Optimal
    @fact getValue(y) --> 1.5
    @fact getValue(z) --> 0
    @fact getValue(w) --> 0.5
    @fact getValue(v) --> 1
    @fact getObjectiveValue(m) --> 1.5 + 0 + 0.5 + 1

    @fact solve(m) --> :Optimal
    @fact getValue(y) --> 2
    @fact getValue(z) --> 0
    @fact getValue(w) --> 1
    @fact getValue(v) --> 1
    @fact getObjectiveValue(m) --> 2 + 0 + 1 + 1

    @defVar(m, 1 <= x <= 2, SemiCont)
    @defVar(m, -2 <= t <= -1, SemiInt)

    addSOS1(m, [x, 2y, 3z, 4w, 5v, 6t])
    @setObjective(m, Min, x + y + z + w + v - t)

    @fact solve(m, relaxation=true) --> :Optimal

    @fact getValue(x) --> 0
    @fact getValue(y) --> 1.5
    @fact getValue(z) --> 0
    @fact getValue(w) --> 0.5
    @fact getValue(v) --> 1
    @fact getValue(t) --> 0
    @fact getObjectiveValue(m) --> 0 + 1.5 + 0 + 0.5 + 1 + 0
end
