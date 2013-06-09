# variable.jl
# Test coverage for Variable

# Test setters/getters
# Name
m = Model("max")
@defVar(m, 0 <= x <= 2)
@assert getName(x) == "x"
setName(x, "x2")
@assert getName(x) == "x2"
setName(x, "")
@assert getName(x) == "_col1"


# Bounds
@assert getLower(x) == 0
@assert getUpper(x) == 2
setLower(x, 1)
@assert getLower(x) == 1
setUpper(x, 3)
@assert getUpper(x) == 3
