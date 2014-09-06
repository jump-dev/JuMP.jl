# model.jl
# Test coverage for Model - writing it to files, and solving

modPath = joinpath(Pkg.dir("JuMP"),"test","mod")

# Check error cases
let
    @test_throws ErrorException Model(solver=:Foo)
    modErr = Model()
    @test_throws ErrorException setObjectiveSense(modErr, :Maximum)
    @defVar(modErr, errVar)
    @test isnan(getValue(errVar))
    @test_throws ErrorException getDual(errVar)

    modErr = Model()
    @defVar(modErr, x, Bin)
    @setObjective(modErr, Max, x)
    con = @addConstraint(modErr, x <= 0.5)
    solve(modErr)
    @test_throws ErrorException getDual(con)

    modErr = Model()
    @defVar(modErr, 0 <= x <= 1)
    @setObjective(modErr, Max, x)
    @addConstraint(modErr, x <= -1)
    solve(modErr)
    @test isnan(getValue(x))
end

###############################################################################
# Test Model A
#####################################################################
# Build it
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
    @test strip(line) == strip(modALP[lineInd])
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
    @test chomp(line) == modAMPS[lineInd]
    lineInd += 1
end
close(modAfp)
#####################################################################
# Test printing
s = IOBuffer()
show(s, modA)
@test takebuf_string(s) == "Maximization problem with:\n * 3 linear constraints\n * 7 variables: 1 integer\nSolver set to Default"
#####################################################################
# Test solution (MIP)
status = solve(modA)
@test status == :Optimal
@test_approx_eq_eps modA.objVal (1.0+4.833334) 1e-6
@test_approx_eq_eps getValue(x) 1.0 1e-6
@test_approx_eq_eps getValue(y) 1.0 1e-6
@test_approx_eq_eps getValue(z) 4.0 1e-6
@test_approx_eq_eps getValue(r)[3] 0.5 1e-6
@test_approx_eq_eps getValue(r)[4] 0.0 1e-6
@test_approx_eq_eps getValue(r)[5] 0.0 1e-6
@test_approx_eq_eps getValue(r)[6] 6.0 1e-6
#####################################################################
# Test solution (LP)
modA.colCat[2] = :Cont
status = solve(modA)
@test status == :Optimal  
@test_approx_eq_eps modA.objVal 5.844611528822055 1e-6
@test_approx_eq_eps getValue(x) 0.9774436090225564 1e-6
@test_approx_eq_eps getValue(y) 1.0225563909774436 1e-6
@test_approx_eq_eps getValue(z) 4.0 1e-6
@test_approx_eq_eps getValue(r)[3] 0.5112781954887218 1e-6
@test_approx_eq_eps getValue(r)[4] 0.0 1e-6
@test_approx_eq_eps getValue(r)[5] 0.0 1e-6
@test_approx_eq_eps getValue(r)[6] 6.0 1e-6

# reduced costs
@test_approx_eq_eps getDual(x) 0.0 1e-6
@test_approx_eq_eps getDual(y) 0.0 1e-6
@test_approx_eq_eps getDual(z) 1.0714285714285714 1e-6
@test_approx_eq_eps getDual(r)[3] 0.0 1e-6
@test_approx_eq_eps getDual(r)[4] -1.0 1e-6
@test_approx_eq_eps getDual(r)[5] -1.0 1e-6
@test_approx_eq_eps getDual(r)[6] 0.03759398496240601 1e-6

# row duals
@test_approx_eq_eps getDual(constraints)[1] -0.333333 1e-6
@test_approx_eq_eps getDual(constraints)[2] 1.0 1e-6
@test_approx_eq_eps getDual(constraints)[3] 0.07142857142857142 1e-6


println("  !!TODO: test external solvers for reading LP and MPS files")


# Getter/setters
@test getNumVars(modA) == 7
@test getNumConstraints(modA) == 3
@test_approx_eq_eps getObjectiveValue(modA) 5.844611528822055 1e-6
@test getObjectiveSense(modA) == :Max
setObjectiveSense(modA, :Min)
@test getObjectiveSense(modA) == :Min

#####################################################################
# Test binary variable handling
let
    modB = Model()
    @defVar(modB, x, Bin)
    @setObjective(modB, Max, x)
    @addConstraint(modB, x <= 10)
    status = solve(modB)
    @test status == :Optimal
    @test_approx_eq_eps getValue(x) 1.0 1e-6
end

#####################################################################
# Test model copying
let
    source = Model()
    @defVar(source, 2 <= x <= 5)
    @defVar(source, 0 <= y <= 1, Int)
    @setObjective(source, Max, 3*x + 1*y)
    @addConstraint(source, x + 2.0*y <= 6)
    addConstraint(source, x*x <= 1)

    dest = copy(source)

    # Obj
    @setObjective(source, Max, 1x)
    @test length(source.obj.aff.coeffs) == 1
    @test length(dest.obj.aff.coeffs) == 2
    @setObjective(dest, Max, 1x)
    @test length(source.obj.aff.coeffs) == 1
    @test length(dest.obj.aff.coeffs) == 1
    @setObjective(dest, Max, 3*x + 1*y)
    @test length(source.obj.aff.coeffs) == 1
    @test length(dest.obj.aff.coeffs) == 2

    # Constraints
    source.linconstr[1].ub = 5.0
    @test dest.linconstr[1].ub == 6.0
    source.quadconstr[1].sense = :(>=)
    @test dest.quadconstr[1].sense == :(<=)
end  

#####################################################################
# Test variable/model "hygiene"
let 
    modA = Model()
    modB = Model()
    @defVar(modA, x)
    @defVar(modB, y)
    @addConstraint(modA, x+y == 1)
    @test_throws ErrorException solve(modA)

    addConstraint(modB, x*y >= 1)
    @test_throws ErrorException solve(modB)
end

######################################################################
# Test nonlinear printing
let
    mod = Model()
    @defVar(mod, x[1:5])
    @addNLConstraint(mod, x[1]*x[2] == 1)
    @addNLConstraint(mod, x[3]*x[4] == 1)
    @addNLConstraint(mod, x[5]*x[1] == 1)
    @setNLObjective(mod, Min, x[1]*x[3])
    str = string(mod)
    @test str == "Min (nonlinear expression)\nSubject to \n3 nonlinear constraints\nx[i], for all i in {1..5} free\n"
end


######################################################################
# Test NaN checking
let
    mod = Model()
    @defVar(mod, x)
    @setObjective(mod, Min, NaN*x)
    @test_throws ErrorException solve(mod)
    setObjective(mod, :Min, NaN*x^2)
    @test_throws ErrorException solve(mod)
    @setObjective(mod, Min, x)
    @addConstraint(mod, Min, NaN*x == 0)
    @test_throws ErrorException solve(mod)
end

######################################################################
# Test all MPS paths
let
    mod = Model()
    @defVar(mod, free_var)
    @defVar(mod, int_var, Bin)
    @defVar(mod, low_var >= 5)
    @addConstraint(mod, free_var == int_var)
    @addConstraint(mod, free_var - int_var >= 0)
    setObjective(mod, :Max, free_var*int_var + low_var)
    writeMPS(mod,"test.mps")
end

######################################################################
# Test all LP paths
let
    mod = Model()
    @defVar(mod, free_var)
    setObjective(mod, :Max, free_var*free_var)
    @test_throws MethodError writeLP("test.lp")
    @setObjective(mod, Max, free_var)
    @addConstraint(mod, free_var - 2*free_var == 0)
    @addConstraint(mod, free_var + 2*free_var >= 1)
    writeLP(mod,"test.lp")
end

######################################################################
# Test semi-continuous variables
if Pkg.installed("Gurobi") != nothing
    let 
        using Gurobi
        mod = Model(solver=GurobiSolver())
        @defVar(mod, x >= 3, SemiCont)
        @defVar(mod, y >= 2, SemiCont)
        @addConstraint(mod, x + y >= 1)
        @setObjective(mod, Min, x+y)
        solve(mod)
        @test getValue(x) == 0.0
        @test getValue(y) == 2.0
    end
end
if Pkg.installed("CPLEX") != nothing
    let 
        using CPLEX
        mod = Model(solver=CplexSolver())
        @defVar(mod, x >= 3, SemiCont)
        @defVar(mod, y >= 2, SemiCont)
        @addConstraint(mod, x + y >= 1)
        @setObjective(mod, Min, x+y)
        solve(mod)
        @test getValue(x) == 0.0
        @test getValue(y) == 2.0
    end
end

######################################################################
# Test semi-integer variables
if Pkg.installed("Gurobi") != nothing
    let 
        using Gurobi
        mod = Model(solver=GurobiSolver())
        @defVar(mod, x >= 3, SemiInt)
        @defVar(mod, y >= 2, SemiInt)
        @addConstraint(mod, x + y >= 2.5)
        @setObjective(mod, Min, x+1.1y)
        solve(mod)
        @test getValue(x) == 3.0
        @test getValue(y) == 0.0
    end
end
if Pkg.installed("CPLEX") != nothing
    let 
        using CPLEX
        mod = Model(solver=CplexSolver())
        @defVar(mod, x >= 3, SemiInt)
        @defVar(mod, y >= 2, SemiInt)
        @addConstraint(mod, x + y >= 2.5)
        @setObjective(mod, Min, x+1.1y)
        solve(mod)
        @test getValue(x) == 3.0
        @test getValue(y) == 0.0
    end
end
