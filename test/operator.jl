# operator.jl
# Test coverage for all operator overloading

m = Model("max")

# Different objects that must all interact:
# 1. Number
# 2. Variable
# 3. AffExpr
# 4. QuadExpr
# 5. Constraint (for comparison ops)
