# operator.jl
# Test coverage for all operator overloading

m = Model("max")
@defVar(m, w)
@defVar(m, x)
@defVar(m, y)
@defVar(m, z)
aff = 7.1 * x + 2.5
@test affToStr(aff) == "7.1 x + 2.5"
q = 2.5 * y * z + aff
@test quadToStr(q) == "2.5 z*y + 7.1 x + 2.5"

# Different objects that must all interact:
# 1. Number
# 2. Variable
# 3. AffExpr
# 4. QuadExpr
# 5. Constraint (for comparison ops)

# 1. Number tests
# 1-1 Number--Number - nope!
# 1-2 Number--Variable
@test affToStr(4.13 + w) == "1.0 w + 4.13"
@test affToStr(3.16 - w) == "-1.0 w + 3.16"
@test affToStr(5.23 * w) == "5.23 w"
@test_fails 2.94 / w
# 1-3 Number--AffExpr
@test affToStr(1.5 + aff) == "7.1 x + 4.0"
@test affToStr(1.5 - aff) == "-7.1 x + -1.0"
@test affToStr(2.0 * aff) == "14.2 x + 5.0"
@test_fails 2.0 / aff
# 1-4 Number--QuadExpr
@test quadToStr(1.5 + q) == "2.5 z*y + 7.1 x + 4.0"
@test quadToStr(1.5 - q) == "-2.5 z*y + -7.1 x + -1.0"
@test quadToStr(2.0 * q) == "5.0 z*y + 14.2 x + 5.0"
@test_fails 2.0 / q

# 2. Variable tests
# 2-1 Variable--Number
@test affToStr(w + 4.13) == "1.0 w + 4.13"
@test affToStr(w - 4.13) == "1.0 w + -4.13"
@test affToStr(w * 4.13) == "4.13 w"
@test affToStr(w / 2.00) == "0.5 w"
# 2-2 Variable--Variable
@test affToStr(w + x) == "1.0 w + 1.0 x"
@test affToStr(w - x) == "1.0 w + -1.0 x"
println("  TODO: Trailing + on quadToStr? Maybe ok cos of aff")
println("        e.g. quadToStr(w * x) = 1.0 w*x + 0.0")
@test quadToStr(w * x) == "1.0 w*x + 0.0"
@test_fails w / x
# 2-3 Variable--AffExpr
@test affToStr(z + aff) == "7.1 x + 1.0 z + 2.5"
@test affToStr(z - aff) == "-7.1 x + 1.0 z + -2.5"
@test quadToStr(z * aff) == "7.1 z*x + 2.5 z"
@test_fails z / aff
# 2-4 Variable--QuadExpr
@test quadToStr(w + q) == "2.5 z*y + 7.1 x + 1.0 w + 2.5"
@test quadToStr(w - q) == "-2.5 z*y + -7.1 x + 1.0 w + -2.5"
@test_fails w*q
@test_fails w/q
