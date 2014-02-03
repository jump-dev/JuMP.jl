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
