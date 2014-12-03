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

# To ensure the tests work on Windows and Linux/OSX, we need
# to use the correct comparison operators
const leq = JuMP.repl_leq
const geq = JuMP.repl_geq
const  eq = JuMP.repl_eq

modPath = joinpath(Pkg.dir("JuMP"),"test","mod")

facts("[model] Check error cases") do
    @fact_throws Model(solver=:Foo)
    modErr = Model()
    @fact_throws setObjectiveSense(modErr, :Maximum)
    @defVar(modErr, errVar)
    @fact getValue(errVar) => isnan
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
    @fact getValue(x) => isnan
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
    constraints[3] = @addConstraint(modA, 7.0*y <= z + r[6]/1.9)
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
    "c4: 7 VAR2 + -1 VAR3 + -.5263157894736842 VAR7 <= 0",
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
        @fact strip(line) => strip(modALP[lineInd])
        lineInd += 1
    end
    close(modAfp)
    #####################################################################
    # Test MPS writer
    writeMPS(modA, modPath * "A.mps")
    modAMPS = ASCIIString[
    "NAME   MathProgModel",
    "ROWS",
    " N  CON4",
    " E  CON1",
    " L  CON2",
    " L  CON3",
    "COLUMNS",
    "    VAR1  CON1  1",
    "    VAR1  CON2  .5",
    "    VAR1  CON4  -.16666666666666666",
    "    MARKER    'MARKER'                 'INTORG'",
    "    VAR2  CON1  1",
    "    VAR2  CON3  7",
    "    VAR2  CON4  -.16666666666666666",
    "    MARKER    'MARKER'                 'INTEND'",
    "    VAR3  CON3  -1",
    "    VAR3  CON4  -1",
    "    VAR4  CON2  1",
    "    VAR4  CON4  -1",
    "    VAR5  CON2  1",
    "    VAR6  CON2  1",
    "    VAR7  CON3  -.5263157894736842",
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
        @fact chomp(line) => modAMPS[lineInd]
        lineInd += 1
    end
    close(modAfp)

    # Getter/setters
    @fact getNumVars(modA) => 7
    @fact getNumConstraints(modA) => 3
    @fact getObjectiveSense(modA) => :Max
    setObjectiveSense(modA, :Min)
    @fact getObjectiveSense(modA) => :Min
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
 
    @fact solve(modA)       => :Optimal
    @fact modA.objVal       => roughly(1.0+4.833334, 1e-6)
    @fact getValue(x)       => roughly(1.0, 1e-6)
    @fact getValue(y)       => roughly(1.0, 1e-6)
    @fact getValue(z)       => roughly(4.0, 1e-6)
    @fact getValue(r)[3]    => roughly(0.5, 1e-6)
    @fact getValue(r)[4]    => roughly(0.0, 1e-6)
    @fact getValue(r)[5]    => roughly(0.0, 1e-6)
    @fact getValue(r)[6]    => roughly(6.0, 1e-6)
end # solver context
end # loop over solvers
end # facts block



facts("[model] Test solving an LP") do
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
    @fact solve(modA) => :Optimal  
    @fact getObjectiveValue(modA) => roughly(5.8446115, 1e-6)
    @fact getValue(x)       => roughly(0.9774436, 1e-6)
    @fact getValue(y)       => roughly(1.0225563, 1e-6)
    @fact getValue(z)       => roughly(4.0, 1e-6)
    @fact getValue(r)[3]    => roughly(0.5112781, 1e-6)
    @fact getValue(r)[4]    => roughly(0.0, 1e-6)
    @fact getValue(r)[5]    => roughly(0.0, 1e-6)
    @fact getValue(r)[6]    => roughly(6.0, 1e-6)

    # Reduced costs
    @fact getDual(x)    => roughly( 0.0, 1e-6)
    @fact getDual(y)    => roughly( 0.0, 1e-6)
    @fact getDual(z)    => roughly( 1.0714286, 1e-6)
    @fact getDual(r)[3] => roughly( 0.0, 1e-6)
    @fact getDual(r)[4] => roughly(-1.0, 1e-6)
    @fact getDual(r)[5] => roughly(-1.0, 1e-6)
    @fact getDual(r)[6] => roughly( 0.03759398, 1e-6)

    # Row duals
    @fact getDual(cons)[1] => roughly(-0.333333, 1e-6)
    @fact getDual(cons)[2] => roughly( 1.0, 1e-6)
    @fact getDual(cons)[3] => roughly( 0.0714286, 1e-6)
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
    @fact status => :Optimal
    @fact getValue(x) => roughly(1.0, 1e-6)
end
end
end


facts("[model] Test model copying") do
    source = Model()
    @defVar(source, 2 <= x <= 5)
    @defVar(source, 0 <= y <= 1, Int)
    @setObjective(source, Max, 3*x + 1*y)
    @addConstraint(source, x + 2.0*y <= 6)
    addConstraint(source, x*x <= 1)

    dest = copy(source)

    # Obj
    @setObjective(source, Max, 1x)
    @fact length(source.obj.aff.coeffs) => 1
    @fact length(dest.obj.aff.coeffs) => 2
    @setObjective(dest, Max, 1x)
    @fact length(source.obj.aff.coeffs) => 1
    @fact length(dest.obj.aff.coeffs) => 1
    @setObjective(dest, Max, 3*x + 1*y)
    @fact length(source.obj.aff.coeffs) => 1
    @fact length(dest.obj.aff.coeffs) => 2

    # Constraints
    source.linconstr[1].ub = 5.0
    @fact dest.linconstr[1].ub => 6.0
    source.quadconstr[1].sense = :(>=)
    @fact dest.quadconstr[1].sense => :(<=)
end


facts("[model] Test variable/model 'hygiene'") do
    modA = Model()
    modB = Model()
    @defVar(modA, x)
    @defVar(modB, y)
    @addConstraint(modA, x+y == 1)
    @fact_throws solve(modA)

    addConstraint(modB, x*y >= 1)
    @fact_throws solve(modB)
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
    @defVar(mod, 0 <= z1 <= 0, 0.0, con, [1.0,-2.0]) # coverage for deprecated syntax
    @defVar(mod, 0 <= z1 <= 1, objective=10.0, inconstraints=con, coefficients=[1.0,-2.0])
    @defVar(mod, 0 <= z2 <= 1, objective=10.0, inconstraints=Any[con[i] for i in [1:2]], coefficients=[1.0,-2.0])
    @fact solve(mod) => :Optimal
    @fact getValue(z1) => roughly(1.0, 1e-6)
    @fact getValue(z2) => roughly(1.0, 1e-6)
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
    @fact getValue(x) => roughly(0.0, 1e-6)
    @fact getValue(y) => roughly(2.0, 1e-6)
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
    @fact getValue(x) => roughly(3.0, 1e-6)
    @fact getValue(y) => 0.0
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

    @fact solve(modS) => :Optimal
    @fact modS.objVal => roughly(15.0, 1e-6)
    @fact getValue(z) => roughly( 3.0, 1e-6)
    @fact getValue(w) => roughly(12.0, 1e-6)


    m = Model(solver=solver)
    ub = [1,1,2]
    @defVar(m, x[i=1:3] <= ub[i])
    @setObjective(m, Max, 2x[1]+x[2]+x[3])
    addSOS1(m, [x[1],2x[2]])
    addSOS1(m, [x[1],2x[3]])

    @fact solve(m) => :Optimal
    @fact getValue(x)[:] => roughly([0.0,1.0,2.0], 1e-6)
    @fact getObjectiveValue(m) => roughly(3.0, 1e-6)
end; end; end

facts("[model] Test setSolver") do
    m = Model()
    @defVar(m, x[1:5])
    @addConstraint(m, con[i=1:5], x[6-i] == i)
    @fact solve(m) => :Optimal

    for solver in lp_solvers
        setSolver(m, solver)
        @fact m.solver => solver
        @fact m.internalModel == nothing => true
        @fact solve(m) => :Optimal
        @fact m.solver => solver
        @fact m.internalModel == nothing => false
    end
end
