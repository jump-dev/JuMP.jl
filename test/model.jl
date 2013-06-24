# model.jl
# Test coverage for Model - writing it to files, and solving

modPath = joinpath(Pkg.dir("MathProg"),"test","mod")

###############################################################################
# Test Model A
#####################################################################
# Build it
modA = Model(:Max)
@defVar(modA, x >= 0)
@defVar(modA, y <= 5.5, Int)
@defVar(modA, 2 <= z <= 4)
@defVar(modA, 0 <= r[i=1:10] <= i)
@setObjective(modA, ((x + y)/2.0 + 3.0)/3.0 + z + r[3])
@addConstraint(modA, -1 <= x+y <= 4)
@addConstraint(modA, sum{r[i],i=3:5} <= (2 - x)/2.0)
#@addConstraint(modA, y*(2+5.0) <= z + r[9])
@addConstraint(modA, 7.0*y <= z + r[9])
#####################################################################
# Test LP writer
writeLP(modA, modPath * "A.lp")
modALP = ASCIIString[
"Maximize",
"obj: 0.166667 VAR1 + 0.166667 VAR2 + 1.000000 VAR3 + 1.000000 VAR6",
"Subject To",
"c1: 1.000000 VAR1 + 1.000000 VAR2 >= -1.000000",
"c2: 1.000000 VAR1 + 1.000000 VAR2 <= 4.000000",
"c3: 1.000000 VAR6 + 1.000000 VAR7 + 1.000000 VAR8 + 0.500000 VAR1 <= 1.000000",
"c4: 7.000000 VAR2 + -1.000000 VAR3 + -1.000000 VAR12 <= 0.000000",
"Bounds",
"0.000000 <= VAR1 <= +inf",
"-inf <= VAR2 <= 5.500000",
"2.000000 <= VAR3 <= 4.000000",
"0.000000 <= VAR4 <= 1.000000",
"0.000000 <= VAR5 <= 2.000000",
"0.000000 <= VAR6 <= 3.000000",
"0.000000 <= VAR7 <= 4.000000",
"0.000000 <= VAR8 <= 5.000000",
"0.000000 <= VAR9 <= 6.000000",
"0.000000 <= VAR10 <= 7.000000",
"0.000000 <= VAR11 <= 8.000000",
"0.000000 <= VAR12 <= 9.000000",
"0.000000 <= VAR13 <= 10.000000",
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
"    VAR1  CON4  0.166667",
"    MARKER    'MARKER'                 'INTORG'",
"    VAR2  CON1  1.000000",
"    VAR2  CON3  7.000000",
"    VAR2  CON4  0.166667",
"    MARKER    'MARKER'                 'INTEND'",
"    VAR3  CON3  -1.000000",
"    VAR3  CON4  1.000000",
"    VAR6  CON2  1.000000",
"    VAR6  CON4  1.000000",
"    VAR7  CON2  1.000000",
"    VAR8  CON2  1.000000",
"    VAR12  CON3  -1.000000",
"RHS",
"    rhs    CON1    -1.000000",
"    rhs    CON2    1.000000",
"    rhs    CON3    0.000000",
"RANGES",
"    rhs    CON1    5.000000",
"BOUNDS",
"  MI BOUND VAR2",
"  UP BOUND VAR2 5.500000",
"  LO BOUND VAR3 2.000000",
"  UP BOUND VAR3 4.000000",
"  UP BOUND VAR4 1.000000",
"  UP BOUND VAR5 2.000000",
"  UP BOUND VAR6 3.000000",
"  UP BOUND VAR7 4.000000",
"  UP BOUND VAR8 5.000000",
"  UP BOUND VAR9 6.000000",
"  UP BOUND VAR10 7.000000",
"  UP BOUND VAR11 8.000000",
"  UP BOUND VAR12 9.000000",
"  UP BOUND VAR13 10.000000",
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
@test_approx_eq modA.objVal (6 + 1.0/6)
@test_approx_eq getValue(x) 0.0
@test_approx_eq getValue(y) 1.0
@test_approx_eq getValue(z) 4.0
@test_approx_eq getValue(r)[1] 0.0
@test_approx_eq getValue(r)[3] 1.0
@test_approx_eq getValue(r)[9] 9.0
@test_approx_eq getValue(r)[8] 0.0
#####################################################################
# Test solution (LP)
modA.colCat[2] = 0
status = solve(modA)
@test status == :Optimal  
@test_approx_eq modA.objVal (1.0 + 5.3095238095)
@test_approx_eq getValue(x) 0.0
@test_approx_eq getValue(y) 1.85714285714286
@test_approx_eq getValue(z) 4.0
@test_approx_eq getValue(r)[1] 0.0
@test_approx_eq getValue(r)[3] 1.0
@test_approx_eq getValue(r)[9] 9.0
@test_approx_eq getValue(r)[8] 0.0
println("  !!TODO: test external solvers for reading LP and MPS files")
