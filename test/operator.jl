# operator.jl
# Test coverage for all operator overloading
using JuMP
using Base.Test

m = Model()
@defVar(m, w)
@defVar(m, x)
@defVar(m, y)
@defVar(m, z)
aff = 7.1 * x + 2.5
@test affToStr(aff) == "7.1 x + 2.5"
aff2 = 1.2 * y + 1.2
@test affToStr(aff2) == "1.2 y + 1.2"
q = 2.5 * y * z + aff
@test quadToStr(q) == "2.5 z*y + 7.1 x + 2.5"
q2 = 8.0 * x * z + aff2
@test quadToStr(q2) == "8.0 z*x + 1.2 y + 1.2"
q3 = 2.0 * x * x + 1.0 * y * y + z + 3.0
@test quadToStr(q3) == "2.0 x² + 1.0 y² + 1.0 z + 3.0"

# Different objects that must all interact:
# 1. Number
# 2. Variable
# 3. AffExpr
# 4. QuadExpr

# 1. Number tests
# 1-1 Number--Number - nope!
# 1-2 Number--Variable
@test affToStr(4.13 + w) == "1.0 w + 4.13"
@test affToStr(3.16 - w) == "-1.0 w + 3.16"
@test affToStr(5.23 * w) == "5.23 w"
@test_throws 2.94 / w
@test conToStr(2.1 <= w) == "1.0 w >= 2.1"
@test conToStr(2.1 == w) == "1.0 w == 2.1"
@test conToStr(2.1 >= w) == "1.0 w <= 2.1"
# 1-3 Number--AffExpr
@test affToStr(1.5 + aff) == "7.1 x + 4.0"
@test affToStr(1.5 - aff) == "-7.1 x - 1.0"
@test affToStr(2.0 * aff) == "14.2 x + 5.0"
@test_throws 2.0 / aff
@test conToStr(1.0 <= aff) == "7.1 x >= -1.5"
@test conToStr(1.0 == aff) == "7.1 x == -1.5"
@test conToStr(1.0 >= aff) == "7.1 x <= -1.5"
# 1-4 Number--QuadExpr
@test quadToStr(1.5 + q) == "2.5 z*y + 7.1 x + 4.0"
@test quadToStr(1.5 - q) == "-2.5 z*y - 7.1 x - 1.0"
@test quadToStr(2.0 * q) == "5.0 z*y + 14.2 x + 5.0"
@test_throws 2.0 / q
@test conToStr(1.0 <= q) == "2.5 z*y + 7.1 x + 1.5 >= 0"
@test conToStr(1.0 == q) == "2.5 z*y + 7.1 x + 1.5 == 0"
@test conToStr(1.0 >= q) == "2.5 z*y + 7.1 x + 1.5 <= 0"

# 2. Variable tests
# 2-1 Variable--Number
@test affToStr(w + 4.13) == "1.0 w + 4.13"
@test affToStr(w - 4.13) == "1.0 w - 4.13"
@test affToStr(w * 4.13) == "4.13 w"
@test affToStr(w / 2.00) == "0.5 w"
@test conToStr(w <= 1.0) == "1.0 w <= 1.0"
@test conToStr(w == 1.0) == "1.0 w == 1.0"
@test conToStr(w >= 1.0) == "1.0 w >= 1.0"
@test conToStr(x*y <= 1.0) == "1.0 x*y - 1.0 <= 0"
@test conToStr(x*y == 1.0) == "1.0 x*y - 1.0 == 0"
@test conToStr(x*y >= 1.0) == "1.0 x*y - 1.0 >= 0"
# 2-2 Variable--Variable
@test affToStr(w + x) == "1.0 w + 1.0 x"
@test affToStr(w - x) == "1.0 w - 1.0 x"
@test quadToStr(w * x) == "1.0 w*x"
@test affToStr(x - x) == "0.0"
@test_throws w / x
@test conToStr(w <= x) == "1.0 w - 1.0 x <= 0.0"
@test conToStr(w == x) == "1.0 w - 1.0 x == 0.0"
@test conToStr(w >= x) == "1.0 w - 1.0 x >= 0.0"
@test conToStr(y*z <= x) == "1.0 y*z - 1.0 x <= 0"
@test conToStr(y*z == x) == "1.0 y*z - 1.0 x == 0"
@test conToStr(y*z >= x) == "1.0 y*z - 1.0 x >= 0"
@test conToStr(x <= x) == "0.0 <= 0.0"
@test conToStr(x == x) == "0.0 == 0.0"
@test conToStr(x >= x) == "0.0 >= 0.0"
# 2-3 Variable--AffExpr
@test affToStr(z + aff) == "7.1 x + 1.0 z + 2.5"
@test affToStr(z - aff) == "-7.1 x + 1.0 z - 2.5"
@test quadToStr(z * aff) == "7.1 z*x + 2.5 z"
@test_throws z / aff
@test conToStr(z <= aff) == "-7.1 x + 1.0 z <= 2.5"
@test conToStr(z == aff) == "-7.1 x + 1.0 z == 2.5"
@test conToStr(z >= aff) == "-7.1 x + 1.0 z >= 2.5"
@test conToStr(7.1 * x - aff <= 0) == "0.0 <= 2.5"
@test conToStr(7.1 * x - aff == 0) == "0.0 == 2.5"
@test conToStr(7.1 * x - aff >= 0) == "0.0 >= 2.5"
# 2-4 Variable--QuadExpr
@test quadToStr(w + q) == "2.5 z*y + 7.1 x + 1.0 w + 2.5"
@test quadToStr(w - q) == "-2.5 z*y - 7.1 x + 1.0 w - 2.5"
@test_throws w*q
@test_throws w/q
@test conToStr(w <= q) == "-2.5 z*y - 7.1 x + 1.0 w - 2.5 <= 0"
@test conToStr(w == q) == "-2.5 z*y - 7.1 x + 1.0 w - 2.5 == 0"
@test conToStr(w >= q) == "-2.5 z*y - 7.1 x + 1.0 w - 2.5 >= 0"

# 3. AffExpr tests
# 3-1 AffExpr--Number
@test affToStr(aff + 1.5) == "7.1 x + 4.0"
@test affToStr(aff - 1.5) == "7.1 x + 1.0"
@test affToStr(aff * 2.0) == "14.2 x + 5.0"
@test affToStr(aff / 2.0) == "3.55 x + 1.25"
@test conToStr(aff <= 1.0) == "7.1 x <= -1.5"
@test conToStr(aff == 1.0) == "7.1 x == -1.5"
@test conToStr(aff >= 1.0) == "7.1 x >= -1.5"

# 3-2 AffExpr--Variable
@test affToStr(aff + z) == "7.1 x + 1.0 z + 2.5"
@test affToStr(aff - z) == "7.1 x - 1.0 z + 2.5"
@test quadToStr(aff * z) == "7.1 z*x + 2.5 z"
@test_throws aff/z
@test conToStr(aff <= z) == "7.1 x - 1.0 z <= -2.5"
@test conToStr(aff == z) == "7.1 x - 1.0 z == -2.5"
@test conToStr(aff >= z) == "7.1 x - 1.0 z >= -2.5"
@test conToStr(aff - 7.1 * x <= 0) == "0.0 <= -2.5"
@test conToStr(aff - 7.1 * x == 0) == "0.0 == -2.5"
@test conToStr(aff - 7.1 * x >= 0) == "0.0 >= -2.5"


# 3-3 AffExpr--AffExpr
@test affToStr(aff + aff2) == "7.1 x + 1.2 y + 3.7"
@test affToStr(aff - aff2) == "7.1 x - 1.2 y + 1.3"
@test quadToStr(aff * aff2) == "8.52 x*y + 3.0 y + 8.52 x + 3.0"
@test_throws aff/aff2
@test conToStr(aff <= aff2) == "7.1 x - 1.2 y <= -1.3"
@test conToStr(aff == aff2) == "7.1 x - 1.2 y == -1.3"
@test conToStr(aff >= aff2) == "7.1 x - 1.2 y >= -1.3"
@test conToStr(aff-aff <= 0) == "0.0 <= 0.0"
@test conToStr(aff-aff == 0) == "0.0 == 0.0"
@test conToStr(aff-aff >= 0) == "0.0 >= 0.0"
# 3-4 AffExpr--QuadExpr
@test quadToStr(aff2 + q) == "2.5 z*y + 1.2 y + 7.1 x + 3.7"
@test quadToStr(aff2 - q) == "-2.5 z*y + 1.2 y - 7.1 x - 1.3"
@test_throws aff2 * q
@test_throws aff2 / q
@test conToStr(aff2 <= q) == "-2.5 z*y + 1.2 y - 7.1 x - 1.3 <= 0"
@test conToStr(aff2 == q) == "-2.5 z*y + 1.2 y - 7.1 x - 1.3 == 0"
@test conToStr(aff2 >= q) == "-2.5 z*y + 1.2 y - 7.1 x - 1.3 >= 0"

# 4. QuadExpr
# 4-1 QuadExpr--Number
@test quadToStr(q + 1.5) == "2.5 z*y + 7.1 x + 4.0"
@test quadToStr(q - 1.5) == "2.5 z*y + 7.1 x + 1.0"
@test quadToStr(q * 2.0) == "5.0 z*y + 14.2 x + 5.0"
@test quadToStr(q / 2.0) == "1.25 z*y + 3.55 x + 1.25"
@test conToStr(q >= 1.0) == "2.5 z*y + 7.1 x + 1.5 >= 0"
@test conToStr(q == 1.0) == "2.5 z*y + 7.1 x + 1.5 == 0"
@test conToStr(q <= 1.0) == "2.5 z*y + 7.1 x + 1.5 <= 0"
# 4-2 QuadExpr--Variable
@test quadToStr(q + w) == "2.5 z*y + 7.1 x + 1.0 w + 2.5"
@test quadToStr(q - w) == "2.5 z*y + 7.1 x - 1.0 w + 2.5"
@test_throws w*q
@test_throws w/q
@test conToStr(q <= w) == "2.5 z*y + 7.1 x - 1.0 w + 2.5 <= 0"
@test conToStr(q == w) == "2.5 z*y + 7.1 x - 1.0 w + 2.5 == 0"
@test conToStr(q >= w) == "2.5 z*y + 7.1 x - 1.0 w + 2.5 >= 0"
# 4-3 QuadExpr--AffExpr
@test quadToStr(q + aff2) == "2.5 z*y + 7.1 x + 1.2 y + 3.7"
@test quadToStr(q - aff2) == "2.5 z*y + 7.1 x - 1.2 y + 1.3"
@test_throws q * aff2
@test_throws q / aff2
@test conToStr(q <= aff2) == "2.5 z*y + 7.1 x - 1.2 y + 1.3 <= 0"
@test conToStr(q == aff2) == "2.5 z*y + 7.1 x - 1.2 y + 1.3 == 0"
@test conToStr(q >= aff2) == "2.5 z*y + 7.1 x - 1.2 y + 1.3 >= 0"
# 4-4 QuadExpr--QuadExpr
@test quadToStr(q + q2) == "2.5 z*y + 8.0 z*x + 7.1 x + 1.2 y + 3.7"
@test quadToStr(q - q2) == "2.5 z*y - 8.0 z*x + 7.1 x - 1.2 y + 1.3"
@test_throws q * q2
@test_throws q / q2
@test conToStr(q <= q2) == "2.5 z*y - 8.0 z*x + 7.1 x - 1.2 y + 1.3 <= 0"
@test conToStr(q == q2) == "2.5 z*y - 8.0 z*x + 7.1 x - 1.2 y + 1.3 == 0"
@test conToStr(q >= q2) == "2.5 z*y - 8.0 z*x + 7.1 x - 1.2 y + 1.3 >= 0"
