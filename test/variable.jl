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

# Test setters/getters
# Name
m = Model()
@defVar(m, 0 <= x <= 2)
@defVar(m, y, Bin)
@test getName(x) == "x"
setName(x, "x2")
@test getName(x) == "x2"
setName(x, "")
@test getName(x) == "_col1"
@defVar(m, z[1:2,3:5])
@defVar(m, w[3:9,["red","blue","green"]])
rng = 2:5
@defVar(m, v[rng,rng,rng,rng,rng,rng,rng])
@test getName(z[1,3]) == "z[1,3]"
@test getName(z[2,4]) == "z[2,4]"
@test getName(z[2,5]) == "z[2,5]"
@test getName(w[7,"green"]) == "w[7,green]"
@test getName(v[4,5,2,3,2,2,4]) == "v[4,5,2,3,2,2,4]"

# Bounds
@test getLower(x) == 0
@test getUpper(x) == 2
setLower(x, 1)
@test getLower(x) == 1
setUpper(x, 3)
@test getUpper(x) == 3
@test getLower(y) == 0
@test getUpper(y) == 1

# Test long-form bounds printing code in print.jl
mprint = Model()
@defVar(mprint, 0 <= a[3:6] <= 5)
@test JuMP.dictstring(a, :REPL)   == "0 ≤ a[i] ≤ 5, for all i in {3..6}"
@test JuMP.dictstring(a, :IJulia) == "0 \\leq a_{i} \\leq 5 \\quad \\forall i \\in \\{ 3..6 \\}"
@defVar(mprint, b[-2:3] >= 2)
@test JuMP.dictstring(b, :REPL)   == "b[i] ≥ 2, for all i in {-2..3}"
@test JuMP.dictstring(b, :IJulia) == "b_{i} \\geq 2 \\quad \\forall i \\in \\{ -2..3 \\}"
@defVar(mprint, c[-2:3] <= 7)
@test JuMP.dictstring(c, :REPL)   == "c[i] ≤ 7, for all i in {-2..3}"
@test JuMP.dictstring(c, :IJulia) == "c_{i} \\leq 7 \\quad \\forall i \\in \\{ -2..3 \\}"
@defVar(mprint, 0 <= d[-5:2] <= 5, Int)
@test JuMP.dictstring(d, :REPL)   == "0 ≤ d[i] ≤ 5, for all i in {-5..2}, integer"
@test JuMP.dictstring(d, :IJulia) == "0 \\leq d_{i} \\leq 5 \\quad \\forall i \\in \\{ -5..2 \\}, integer"
@defVar(mprint, z[-5:2], Bin)
@test JuMP.dictstring(z, :REPL)   == "z[i], for all i in {-5..2}, binary"
@test JuMP.dictstring(z, :IJulia) == "z_{i} \\quad \\forall i \\in \\{ -5..2 \\}, binary"
@defVar(mprint, f[6:9] >= 3, Int)
@test JuMP.dictstring(f, :REPL)   == "f[i] ≥ 3, for all i in {6..9}, integer"
@test JuMP.dictstring(f, :IJulia) == "f_{i} \\geq 3 \\quad \\forall i \\in \\{ 6..9 \\}, integer"
@defVar(mprint, g[1:5:3] >= 0)
@test JuMP.dictstring(g, :REPL)   == "g[i] ≥ 0, for all i in {1}"
@test JuMP.dictstring(g, :IJulia) == "g_{i} \\geq 0 \\quad \\forall i \\in \\{ 1 \\}"
@defVar(mprint, h[0:2:2] >= 0)
@test JuMP.dictstring(h, :REPL)   == "h[i] ≥ 0, for all i in {0,2}"
@test JuMP.dictstring(h, :IJulia) == "h_{i} \\geq 0 \\quad \\forall i \\in \\{ 0,2 \\}"
@defVar(mprint, i[0:2:4] >= 0)
@test JuMP.dictstring(i, :REPL)   == "i[i] ≥ 0, for all i in {0,2,4}"
@test JuMP.dictstring(i, :IJulia) == "i_{i} \\geq 0 \\quad \\forall i \\in \\{ 0,2,4 \\}"
@defVar(mprint, j[1:2:20] >= 0)
@test JuMP.dictstring(j, :REPL)   == "j[i] ≥ 0, for all i in {1,3..17,19}"
@test JuMP.dictstring(j, :IJulia) == "j_{i} \\geq 0 \\quad \\forall i \\in \\{ 1,3..17,19 \\}"
@defVar(mprint, k[[:a,:b,:c]] >= 0)
@test JuMP.dictstring(k, :REPL)   == "k[i] ≥ 0, for all i in {a,b,c}"
@test JuMP.dictstring(k, :IJulia) == "k_{i} \\geq 0 \\quad \\forall i \\in \\{ a,b,c \\}"
@defVar(mprint, l[[:apple,:banana,:carrot,:diamonds]] >= 0)
@test JuMP.dictstring(l, :REPL)   == "l[i] ≥ 0, for all i in {apple,banana..}"
@test JuMP.dictstring(l, :IJulia) == "l_{i} \\geq 0 \\quad \\forall i \\in \\{ apple,banana.. \\}"
@defVar(mprint, m[1:5,2:2:6,[:cats,:dogs,:dinosaurs],["gah",5,2,:symbol]] >= 0)
@test JuMP.dictstring(m, :REPL)   == "m[i,j,k,l] ≥ 0, for all i in {1..5}, j in {2,4,6}, k in {cats,dogs..}, l in {gah,5,2,symbol}"
@test JuMP.dictstring(m, :IJulia) == "m_{i,j,k,l} \\geq 0 \\quad \\forall i \\in \\{ 1..5 \\}, j \\in \\{ 2,4,6 \\}, k \\in \\{ cats,dogs.. \\}, l \\in \\{ gah,5,2,symbol \\}"
idx1 = [[1:10]]
idx2 = [[2:2:20]]
idx3 = [[2:2:20],21]
@defVar(mprint, p[idx1,idx2])
@test JuMP.dictstring(p, :REPL)   == "p[i,j], for all i in {1..10}, j in {2,4..18,20} free"
@test JuMP.dictstring(p, :IJulia) == "p_{i,j} \\quad \\forall i \\in \\{ 1..10 \\}, j \\in \\{ 2,4..18,20 \\} free"
@defVar(mprint, q[idx2,idx3])
@test JuMP.dictstring(q, :REPL)   == "q[i,j], for all i in {2,4..18,20}, j in {2,4,6,8,10,12..} free"
@test JuMP.dictstring(q, :IJulia) == "q_{i,j} \\quad \\forall i \\in \\{ 2,4..18,20 \\}, j \\in \\{ 2,4,6,8,10,12.. \\} free"
@defVar(mprint, r[idx1,idx2] >= 2, SemiCont)
@test JuMP.dictstring(r, :REPL)   == "r[i,j] ≥ 2, for all i in {1..10}, j in {2,4..18,20}, semicontinuous"
@test JuMP.dictstring(r, :IJulia) == "r_{i,j} \\geq 2 \\quad \\forall i \\in \\{ 1..10 \\}, j \\in \\{ 2,4..18,20 \\}, semicontinuous"
@defVar(mprint, s[idx2,idx3] <= -3, SemiInt)
@test JuMP.dictstring(s, :REPL)   == "s[i,j] ≤ -3, for all i in {2,4..18,20}, j in {2,4,6,8,10,12..}, semi-integer"
@test JuMP.dictstring(s, :IJulia) == "s_{i,j} \\leq -3 \\quad \\forall i \\in \\{ 2,4..18,20 \\}, j \\in \\{ 2,4,6,8,10,12.. \\}, semi-integer"
@defVar(mprint, t[1:5,1:10])
@test JuMP.dictstring(t, :REPL)   == "t[i,j], for all i in {1..5}, j in {1..10} free"
@test JuMP.dictstring(t, :IJulia) == "t_{i,j} \\quad \\forall i \\in \\{ 1..5 \\}, j \\in \\{ 1..10 \\} free"

# test empty JuMPDict printing (issue #124)
@defVar(mprint, xx[1:0])
@test JuMP.dictstring(xx, :REPL) == JuMP.dictstring(xx, :IJulia) == ""

# Test printing of getValue(JuMPDict{Float64})
valmod = Model()
@defVar(valmod, i*j <= foobar[i=9:10, [:Apple,5,:Banana], j=-1:+1] <= i*j)
solve(valmod)
buf = IOBuffer()
println(buf, getValue(foobar))
result = takebuf_string(buf)
@test result == "foobar\n[ 9,:,:]\n  [ 9, Apple,:]\n    [ 9, Apple,-1] = -9.0\n    [ 9, Apple, 0] = 0.0\n    [ 9, Apple, 1] = 9.0\n  [ 9,     5,:]\n    [ 9,     5,-1] = -9.0\n    [ 9,     5, 0] = 0.0\n    [ 9,     5, 1] = 9.0\n  [ 9,Banana,:]\n    [ 9,Banana,-1] = -9.0\n    [ 9,Banana, 0] = 0.0\n    [ 9,Banana, 1] = 9.0\n[10,:,:]\n  [10, Apple,:]\n    [10, Apple,-1] = -10.0\n    [10, Apple, 0] = 0.0\n    [10, Apple, 1] = 10.0\n  [10,     5,:]\n    [10,     5,-1] = -10.0\n    [10,     5, 0] = 0.0\n    [10,     5, 1] = 10.0\n  [10,Banana,:]\n    [10,Banana,-1] = -10.0\n    [10,Banana, 0] = 0.0\n    [10,Banana, 1] = 10.0\n\n"

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
