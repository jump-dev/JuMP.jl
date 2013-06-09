# variable.jl
# Test coverage for Variable

# Test setters/getters
m = Model("max")
@defVar(m, 0 <= x <= 2)
@assert getName(x) == "x"
setName(x, "x2")
@assert getName(x) == "x2"
# is this functionality a good thing?
setName(x, "")
@assert getName(x) == "_col1"
