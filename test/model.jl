# model.jl
# Test coverage for Model - writing it to files, and solving

modPath = joinpath(Pkg.dir("JuMP"),"test","mod")

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
"obj: 0.166667 VAR1 + 0.166667 VAR2 + 1.000000 VAR3 + 1.000000 VAR4",
"Subject To",
"c1: 1.000000 VAR1 + 1.000000 VAR2 >= 2.000000",
"c2: 1.000000 VAR1 + 1.000000 VAR2 <= 4.000000",
"c3: 1.000000 VAR4 + 1.000000 VAR5 + 1.000000 VAR6 + 0.500000 VAR1 <= 1.000000",
"c4: 7.000000 VAR2 + -1.000000 VAR3 + -0.526316 VAR7 <= 0.000000",
"Bounds",
"0.000000 <= VAR1 <= +inf",
"-inf <= VAR2 <= 5.000000",
"2.000000 <= VAR3 <= 4.000000",
"0.000000 <= VAR4 <= 3.000000",
"0.000000 <= VAR5 <= 4.000000",
"0.000000 <= VAR6 <= 5.000000",
"0.000000 <= VAR7 <= 6.000000",
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
"    VAR1  CON1  1.000000",
"    VAR1  CON2  0.500000",
"    VAR1  CON4  -0.166667",
"    MARKER    'MARKER'                 'INTORG'",
"    VAR2  CON1  1.000000",
"    VAR2  CON3  7.000000",
"    VAR2  CON4  -0.166667",
"    MARKER    'MARKER'                 'INTEND'",
"    VAR3  CON3  -1.000000",
"    VAR3  CON4  -1.000000",
"    VAR4  CON2  1.000000",
"    VAR4  CON4  -1.000000",
"    VAR5  CON2  1.000000",
"    VAR6  CON2  1.000000",
"    VAR7  CON3  -0.526316",
"RHS",
"    rhs    CON1    2.000000",
"    rhs    CON2    1.000000",
"    rhs    CON3    0.000000",
"RANGES",
"    rhs    CON1    2.000000",
"BOUNDS",
"  MI BOUND VAR2",
"  UP BOUND VAR2 5.000000",
"  LO BOUND VAR3 2.000000",
"  UP BOUND VAR3 4.000000",
"  UP BOUND VAR4 3.000000",
"  UP BOUND VAR5 4.000000",
"  UP BOUND VAR6 5.000000",
"  UP BOUND VAR7 6.000000",
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
# Test solution (MIP)
status = solve(modA)
@test status == :Optimal
@test_approx_eq_eps modA.objVal (1.0+4.833334) 1e-6
@test_approx_eq getValue(x) 1.0
@test_approx_eq getValue(y) 1.0
@test_approx_eq getValue(z) 4.0
@test_approx_eq getValue(r)[3] 0.5
@test_approx_eq getValue(r)[4] 0.0
@test_approx_eq getValue(r)[5] 0.0
@test_approx_eq getValue(r)[6] 6.0
#####################################################################
# Test solution (LP)
modA.colCat[2] = 0
status = solve(modA)
@test status == :Optimal  
@test_approx_eq_eps modA.objVal 5.844611528822055 1e-6
@test_approx_eq getValue(x) 0.9774436090225564
@test_approx_eq getValue(y) 1.0225563909774436
@test_approx_eq getValue(z) 4.0
@test_approx_eq getValue(r)[3] 0.5112781954887218
@test_approx_eq getValue(r)[4] 0.0
@test_approx_eq getValue(r)[5] 0.0
@test_approx_eq getValue(r)[6] 6.0

# reduced costs
@test_approx_eq getDual(x) 0.0
@test_approx_eq getDual(y) 0.0
@test_approx_eq getDual(z) 1.0714285714285714
@test_approx_eq getDual(r)[3] 0.0
@test_approx_eq getDual(r)[4] -1.0
@test_approx_eq getDual(r)[5] -1.0
@test_approx_eq getDual(r)[6] 0.03759398496240601

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
    @test_approx_eq getValue(x) 1.0
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
    source.quadconstr[1].sense == :(>=)
    @test dest.quadconstr[1].sense == :(<=)
end  

#####################################################################
# Test model printing
let
    modC = Model()
    @defVar(modC, a>=1)
    @defVar(modC, b<=1)
    @defVar(modC, -1<=c<=1)
    @defVar(modC, a1>=1,Int)
    @defVar(modC, b1<=1,Int)
    @defVar(modC, -1<=c1<=1,Int)
    @defVar(modC, x, Bin)
    @defVar(modC, y)
    @defVar(modC, z, Int)
    @setObjective(modC, Max, a - b + 2a1 - 10x)
    @addConstraint(modC, a + b - 10c - 2x + c1 <= 1)

    str = string(modC)

    @test str == "Max 1.0 a - 1.0 b + 2.0 a1 - 10.0 x\nSubject to: \n1.0 a + 1.0 b - 10.0 c - 2.0 x + 1.0 c1 <= 1.0\na ≥ 1.0\nb ≤ 1.0\n-1.0 ≤ c ≤ 1.0\na1 ≥ 1.0, a1 integer\nb1 ≤ 1.0, b1 integer\n-1.0 ≤ c1 ≤ 1.0, c1 integer\nx binary\ny free\nz free integer\n"
end
