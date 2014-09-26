# variable.jl
# Test coverage for Variable

# Constructors
mcon = Model()
@defVar(mcon, nobounds)
@defVar(mcon, lbonly >= 0)
@defVar(mcon, ubonly <= 1)
@defVar(mcon, 0 <= bothb <= 1)
@defVar(mcon, 0 <= onerange[-5:5] <= 10)
@defVar(mcon, onerangeub[-7:1] <= 10, Int)
@defVar(mcon, manyrangelb[0:1,10:20,1:1] >= 2)
@test getLower(manyrangelb[0,15,1]) == 2
s = ["Green","Blue"]
@defVar(mcon, x[-10:10,s] <= 5.5, Int)
@test getUpper(x[-4,"Green"]) == 5.5

# Bounds
m = Model()
@defVar(m, 0 <= x <= 2)
@test getLower(x) == 0
@test getUpper(x) == 2
setLower(x, 1)
@test getLower(x) == 1
setUpper(x, 3)
@test getUpper(x) == 3
@defVar(m, y, Bin)
@test getLower(y) == 0
@test getUpper(y) == 1

# Test printing of getValue(JuMPDict{Float64})
valmod = Model()
@defVar(valmod, i*j <= foobar[i=9:10, [:Apple,5,:Banana], j=-1:+1] <= i*j)
solve(valmod)
buf = IOBuffer()
println(buf, getValue(foobar))
result = takebuf_string(buf)
if hash(5) < hash(:Apple) # hashing is different on some platforms, affects output order
    @test result == "3-dimensional JuMPDict with 18 entries:\n foobar[9,5,-1]       = -9.0\n foobar[9,5,0]        = 0.0\n foobar[9,5,1]        = 9.0\n foobar[9,Apple,-1]   = -9.0\n foobar[9,Apple,0]    = 0.0\n foobar[9,Apple,1]    = 9.0\n foobar[9,Banana,-1]  = -9.0\n foobar[9,Banana,0]   = 0.0\n foobar[9,Banana,1]   = 9.0\n foobar[10,5,-1]      = -10.0\n foobar[10,5,0]       = 0.0\n foobar[10,5,1]       = 10.0\n foobar[10,Apple,-1]  = -10.0\n foobar[10,Apple,0]   = 0.0\n foobar[10,Apple,1]   = 10.0\n foobar[10,Banana,-1] = -10.0\n foobar[10,Banana,0]  = 0.0\n foobar[10,Banana,1]  = 10.0\n"
else
    @test result == "3-dimensional JuMPDict with 18 entries:\n foobar[9,Apple,-1]   = -9.0\n foobar[9,Apple,0]    = 0.0\n foobar[9,Apple,1]    = 9.0\n foobar[9,Banana,-1]  = -9.0\n foobar[9,Banana,0]   = 0.0\n foobar[9,Banana,1]   = 9.0\n foobar[9,5,-1]       = -9.0\n foobar[9,5,0]        = 0.0\n foobar[9,5,1]        = 9.0\n foobar[10,Apple,-1]  = -10.0\n foobar[10,Apple,0]   = 0.0\n foobar[10,Apple,1]   = 10.0\n foobar[10,Banana,-1] = -10.0\n foobar[10,Banana,0]  = 0.0\n foobar[10,Banana,1]  = 10.0\n foobar[10,5,-1]      = -10.0\n foobar[10,5,0]       = 0.0\n foobar[10,5,1]       = 10.0\n"
 end

# Repeated elements in index set (issue #199)
repeatmod = Model()
s = [:x,:x,:y]
@defVar(repeatmod, x[s])
@test getNumVars(repeatmod) == 3

# Test conditionals in variable definition
# condmod = Model()
# @defVar(condmod, x[i=1:10]; iseven(i))
# @defVar(condmod, y[j=1:10,k=3:2:9]; isodd(j+k) && k <= 8)
# @test JuMP.dictstring(x, :REPL)   == "x[i], for all i in {1..10} s.t. iseven(i) free"
# @test JuMP.dictstring(x, :IJulia) == "x_{i} \\quad \\forall i \\in \\{ 1..10 \\} s.t. iseven(i) free"
# @test JuMP.dictstring(y, :REPL)   == "y[i,j], for all i in {1..10}, j in {3,5..7,9} s.t. isodd(j + k) and k <= 8 free"
# @test JuMP.dictstring(y, :IJulia) == "y_{i,j} \\quad \\forall i \\in \\{ 1..10 \\}, j \\in \\{ 3,5..7,9 \\} s.t. isodd(j + k) and k <= 8 free"
# @test string(condmod) == "Min 0\nSubject to \nx[i], for all i in {1..10} s.t. iseven(i) free\ny[i,j], for all i in {1..10}, j in {3,5..7,9} s.t. isodd(j + k) and k <= 8 free\n" 
