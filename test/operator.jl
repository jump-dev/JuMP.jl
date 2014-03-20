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
@test quadToStr(q) == "2.5 y*z + 7.1 x + 2.5"
q2 = 8 * x * z + aff2
@test quadToStr(q2) == "8 x*z + 1.2 y + 1.2"
q3 = 2 * x * x + 1 * y * y + z + 3
@test quadToStr(q3) == "2 x² + y² + z + 3"

# Different objects that must all interact:
# 1. Number
# 2. Variable
# 3. AffExpr
# 4. QuadExpr

# 1. Number tests
# 1-1 Number--Number - nope!
# 1-2 Number--Variable
@test affToStr(4.13 + w) == "w + 4.13"
@test affToStr(3.16 - w) == "-w + 3.16"
@test affToStr(5.23 * w) == "5.23 w"
@test_throws 2.94 / w
@test conToStr(2.1 <= w) == "w >= 2.1"
@test conToStr(2.1 == w) == "w == 2.1"
@test conToStr(2.1 >= w) == "w <= 2.1"
# 1-3 Number--AffExpr
@test affToStr(1.5 + aff) == "7.1 x + 4"
@test affToStr(1.5 - aff) == "-7.1 x - 1"
@test affToStr(2 * aff) == "14.2 x + 5"
@test_throws 2 / aff
@test conToStr(1 <= aff) == "7.1 x >= -1.5"
@test conToStr(1 == aff) == "7.1 x == -1.5"
@test conToStr(1 >= aff) == "7.1 x <= -1.5"
# 1-4 Number--QuadExpr
@test quadToStr(1.5 + q) == "2.5 y*z + 7.1 x + 4"
@test quadToStr(1.5 - q) == "-2.5 y*z - 7.1 x - 1"
@test quadToStr(2 * q) == "5 y*z + 14.2 x + 5"
@test_throws 2 / q
@test conToStr(1 <= q) == "2.5 y*z + 7.1 x + 1.5 >= 0"
@test conToStr(1 == q) == "2.5 y*z + 7.1 x + 1.5 == 0"
@test conToStr(1 >= q) == "2.5 y*z + 7.1 x + 1.5 <= 0"

# 2. Variable tests
# 2-1 Variable--Number
@test affToStr(w + 4.13) == "w + 4.13"
@test affToStr(w - 4.13) == "w - 4.13"
@test affToStr(w * 4.13) == "4.13 w"
@test affToStr(w / 2.00) == "0.5 w"
@test conToStr(w <= 1) == "w <= 1"
@test conToStr(w == 1) == "w == 1"
@test conToStr(w >= 1) == "w >= 1"
@test conToStr(x*y <= 1) == "x*y - 1 <= 0"
@test conToStr(x*y == 1) == "x*y - 1 == 0"
@test conToStr(x*y >= 1) == "x*y - 1 >= 0"
# 2-2 Variable--Variable
@test affToStr(w + x) == "w + x"
@test affToStr(w - x) == "w - x"
@test quadToStr(w * x) == "w*x"
@test affToStr(x - x) == "0"
@test_throws w / x
@test conToStr(w <= x) == "w - x <= 0"
@test conToStr(w == x) == "w - x == 0"
@test conToStr(w >= x) == "w - x >= 0"
@test conToStr(y*z <= x) == "y*z - x <= 0"
@test conToStr(y*z == x) == "y*z - x == 0"
@test conToStr(y*z >= x) == "y*z - x >= 0"
@test conToStr(x <= x) == "0 <= 0"
@test conToStr(x == x) == "0 == 0"
@test conToStr(x >= x) == "0 >= 0"
# 2-3 Variable--AffExpr
@test affToStr(z + aff) == "7.1 x + z + 2.5"
@test affToStr(z - aff) == "-7.1 x + z - 2.5"
@test quadToStr(z * aff) == "7.1 x*z + 2.5 z"
@test_throws z / aff
@test conToStr(z <= aff) == "-7.1 x + z <= 2.5"
@test conToStr(z == aff) == "-7.1 x + z == 2.5"
@test conToStr(z >= aff) == "-7.1 x + z >= 2.5"
@test conToStr(7.1 * x - aff <= 0) == "0 <= 2.5"
@test conToStr(7.1 * x - aff == 0) == "0 == 2.5"
@test conToStr(7.1 * x - aff >= 0) == "0 >= 2.5"
# 2-4 Variable--QuadExpr
@test quadToStr(w + q) == "2.5 y*z + 7.1 x + w + 2.5"
@test quadToStr(w - q) == "-2.5 y*z - 7.1 x + w - 2.5"
@test_throws w*q
@test_throws w/q
@test conToStr(w <= q) == "-2.5 y*z - 7.1 x + w - 2.5 <= 0"
@test conToStr(w == q) == "-2.5 y*z - 7.1 x + w - 2.5 == 0"
@test conToStr(w >= q) == "-2.5 y*z - 7.1 x + w - 2.5 >= 0"

# 3. AffExpr tests
# 3-1 AffExpr--Number
@test affToStr(aff + 1.5) == "7.1 x + 4"
@test affToStr(aff - 1.5) == "7.1 x + 1"
@test affToStr(aff * 2) == "14.2 x + 5"
@test affToStr(aff / 2) == "3.55 x + 1.25"
@test conToStr(aff <= 1) == "7.1 x <= -1.5"
@test conToStr(aff == 1) == "7.1 x == -1.5"
@test conToStr(aff >= 1) == "7.1 x >= -1.5"

# 3-2 AffExpr--Variable
@test affToStr(aff + z) == "7.1 x + z + 2.5"
@test affToStr(aff - z) == "7.1 x - z + 2.5"
@test quadToStr(aff * z) == "7.1 x*z + 2.5 z"
@test_throws aff/z
@test conToStr(aff <= z) == "7.1 x - z <= -2.5"
@test conToStr(aff == z) == "7.1 x - z == -2.5"
@test conToStr(aff >= z) == "7.1 x - z >= -2.5"
@test conToStr(aff - 7.1 * x <= 0) == "0 <= -2.5"
@test conToStr(aff - 7.1 * x == 0) == "0 == -2.5"
@test conToStr(aff - 7.1 * x >= 0) == "0 >= -2.5"


# 3-3 AffExpr--AffExpr
@test affToStr(aff + aff2) == "7.1 x + 1.2 y + 3.7"
@test affToStr(aff - aff2) == "7.1 x - 1.2 y + 1.3"
@test quadToStr(aff * aff2) == "8.52 x*y + 3 y + 8.52 x + 3"
@test_throws aff/aff2
@test conToStr(aff <= aff2) == "7.1 x - 1.2 y <= -1.3"
@test conToStr(aff == aff2) == "7.1 x - 1.2 y == -1.3"
@test conToStr(aff >= aff2) == "7.1 x - 1.2 y >= -1.3"
@test conToStr(aff-aff <= 0) == "0 <= 0"
@test conToStr(aff-aff == 0) == "0 == 0"
@test conToStr(aff-aff >= 0) == "0 >= 0"
# 3-4 AffExpr--QuadExpr
@test quadToStr(aff2 + q) == "2.5 y*z + 1.2 y + 7.1 x + 3.7"
@test quadToStr(aff2 - q) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
@test_throws aff2 * q
@test_throws aff2 / q
@test conToStr(aff2 <= q) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3 <= 0"
@test conToStr(aff2 == q) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3 == 0"
@test conToStr(aff2 >= q) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3 >= 0"

# 4. QuadExpr
# 4-1 QuadExpr--Number
@test quadToStr(q + 1.5) == "2.5 y*z + 7.1 x + 4"
@test quadToStr(q - 1.5) == "2.5 y*z + 7.1 x + 1"
@test quadToStr(q * 2) == "5 y*z + 14.2 x + 5"
@test quadToStr(q / 2) == "1.25 y*z + 3.55 x + 1.25"
@test conToStr(q >= 1) == "2.5 y*z + 7.1 x + 1.5 >= 0"
@test conToStr(q == 1) == "2.5 y*z + 7.1 x + 1.5 == 0"
@test conToStr(q <= 1) == "2.5 y*z + 7.1 x + 1.5 <= 0"
# 4-2 QuadExpr--Variable
@test quadToStr(q + w) == "2.5 y*z + 7.1 x + w + 2.5"
@test quadToStr(q - w) == "2.5 y*z + 7.1 x - w + 2.5"
@test_throws w*q
@test_throws w/q
@test conToStr(q <= w) == "2.5 y*z + 7.1 x - w + 2.5 <= 0"
@test conToStr(q == w) == "2.5 y*z + 7.1 x - w + 2.5 == 0"
@test conToStr(q >= w) == "2.5 y*z + 7.1 x - w + 2.5 >= 0"
# 4-3 QuadExpr--AffExpr
@test quadToStr(q + aff2) == "2.5 y*z + 7.1 x + 1.2 y + 3.7"
@test quadToStr(q - aff2) == "2.5 y*z + 7.1 x - 1.2 y + 1.3"
@test_throws q * aff2
@test_throws q / aff2
@test conToStr(q <= aff2) == "2.5 y*z + 7.1 x - 1.2 y + 1.3 <= 0"
@test conToStr(q == aff2) == "2.5 y*z + 7.1 x - 1.2 y + 1.3 == 0"
@test conToStr(q >= aff2) == "2.5 y*z + 7.1 x - 1.2 y + 1.3 >= 0"
# 4-4 QuadExpr--QuadExpr
@test quadToStr(q + q2) == "8 x*z + 2.5 y*z + 7.1 x + 1.2 y + 3.7"
@test quadToStr(q - q2) == "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3"
@test_throws q * q2
@test_throws q / q2
@test conToStr(q <= q2) == "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 <= 0"
@test conToStr(q == q2) == "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 == 0"
@test conToStr(q >= q2) == "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 >= 0"


# Higher-level operators
# sum
let
    sum_m = Model()
    @defVar(sum_m, 0 <= matrix[1:3,1:3] <= 1)
    # sum(j::JuMPDict{Variable}) 
    @test affToStr(sum(matrix)) == "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"
    # sum(j::JuMPDict{Variable}) in a macro
    @setObjective(sum_m, Max, sum(matrix))
    @test quadToStr(sum_m.obj) == "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"
    solve(sum_m)
    # sum{T<:Real}(j::JuMPDict{T})
    @test_approx_eq_eps sum(getValue(matrix)) 9 1e-6
    # sum(j::Array{Variable})
    @test affToStr(sum(matrix[1:3,1:3])) == affToStr(sum(matrix))
    # sum(affs::Array{AffExpr})
    @test affToStr(sum([2*matrix[i,j] for i in 1:3, j in 1:3])) == "2 matrix[1,1] + 2 matrix[2,1] + 2 matrix[3,1] + 2 matrix[1,2] + 2 matrix[2,2] + 2 matrix[3,2] + 2 matrix[1,3] + 2 matrix[2,3] + 2 matrix[3,3]"
end

# dot
let
    dot_m = Model()
    @defVar(dot_m, 0 <= x[1:3] <= 1)
    c = [1:3]
    @test affToStr(dot(c,x)) == "x[1] + 2 x[2] + 3 x[3]"
    @test affToStr(dot(x,c)) == "x[1] + 2 x[2] + 3 x[3]"

    A = [1 3 ; 2 4]
    @defVar(dot_m, 1 <= y[1:2,1:2] <= 1)
    @test affToStr(dot(A,y)) == "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"
    @test affToStr(dot(y,A)) == "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"

    B = ones(2,2,2)
    @defVar(dot_m, 0 <= z[1:2,1:2,1:2] <= 1)
    @test affToStr(dot(B,z)) == "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"
    @test affToStr(dot(z,B)) == "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"

    @setObjective(dot_m, Max, dot(x, ones(3)) - dot(y, ones(2,2)) )
    solve(dot_m)
    @test_approx_eq_eps dot(c, getValue(x))   6  1e-6
    @test_approx_eq_eps dot(A, getValue(y))  10  1e-6
end

let
    slice_m = Model()
    C = [:cat,:dog]
    @defVar(slice_m, x[-1:1,C])
    catcoef = [1,2,3]
    dogcoef = [3,4,5]

    @test affToStr(dot(catcoef, x[:,:cat])) == "x[-1,cat] + 2 x[0,cat] + 3 x[1,cat]"
    @test affToStr(dot(dogcoef, x[:,:dog])) == "3 x[-1,dog] + 4 x[0,dog] + 5 x[1,dog]"
end

