#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/operator.jl
# Testing operator overloading is correct
#############################################################################
using JuMP, FactCheck

# To ensure the tests work on both Windows and Linux/OSX, we need
# to use the correct comparison operators in the strings
const leq = JuMP.repl_leq
const geq = JuMP.repl_geq
const  eq = JuMP.repl_eq

facts("[operator] Testing basic operator overloads") do
    m = Model()
    @defVar(m, w)
    @defVar(m, x)
    @defVar(m, y)
    @defVar(m, z)
    aff = 7.1 * x + 2.5
    @fact affToStr(aff) => "7.1 x + 2.5"
    aff2 = 1.2 * y + 1.2
    @fact affToStr(aff2) => "1.2 y + 1.2"
    q = 2.5 * y * z + aff
    @fact quadToStr(q) => "2.5 y*z + 7.1 x + 2.5"
    q2 = 8 * x * z + aff2
    @fact quadToStr(q2) => "8 x*z + 1.2 y + 1.2"
    q3 = 2 * x * x + 1 * y * y + z + 3
    @fact quadToStr(q3) => "2 x² + y² + z + 3"

    # Different objects that must all interact:
    # 1. Number
    # 2. Variable
    # 3. AffExpr
    # 4. QuadExpr

    # 1. Number tests
    context("Number--???") do
    # 1-1 Number--Number - nope!
    # 1-2 Number--Variable
    @fact affToStr(4.13 + w) => "w + 4.13"
    @fact affToStr(3.16 - w) => "-w + 3.16"
    @fact affToStr(5.23 * w) => "5.23 w"
    @fact_throws  2.94 / w
    @fact conToStr(2.1 ≤ w) => "w $geq 2.1"
    @fact conToStr(2.1 == w) => "w $eq 2.1"
    @fact conToStr(2.1 ≥ w) => "w $leq 2.1"
    # 1-3 Number--AffExpr
    @fact affToStr(1.5 + aff) => "7.1 x + 4"
    @fact affToStr(1.5 - aff) => "-7.1 x - 1"
    @fact affToStr(2 * aff) => "14.2 x + 5"
    @fact_throws  2 / aff
    @fact conToStr(1 ≤ aff) => "7.1 x $geq -1.5"
    @fact conToStr(1 == aff) => "7.1 x $eq -1.5"
    @fact conToStr(1 ≥ aff) => "7.1 x $leq -1.5"
    # 1-4 Number--QuadExpr
    @fact quadToStr(1.5 + q) => "2.5 y*z + 7.1 x + 4"
    @fact quadToStr(1.5 - q) => "-2.5 y*z - 7.1 x - 1"
    @fact quadToStr(2 * q) => "5 y*z + 14.2 x + 5"
    @fact_throws  2 / q
    @fact conToStr(1 ≤ q) => "2.5 y*z + 7.1 x + 1.5 $geq 0"
    @fact conToStr(1 == q) => "2.5 y*z + 7.1 x + 1.5 $eq 0"
    @fact conToStr(1 ≥ q) => "2.5 y*z + 7.1 x + 1.5 $leq 0"
    end 

    # 2. Variable tests
    context("Variable--???") do
    # 2-0 Variable unary
    @fact (+x) => exactly(x)
    @fact affToStr(-x) => "-x"
    # 2-1 Variable--Number
    @fact affToStr(w + 4.13) => "w + 4.13"
    @fact affToStr(w - 4.13) => "w - 4.13"
    @fact affToStr(w * 4.13) => "4.13 w"
    @fact affToStr(w / 2.00) => "0.5 w"
    @fact conToStr(w ≤ 1) => "w $leq 1"
    @fact conToStr(w == 1) => "w $eq 1"
    @fact conToStr(w ≥ 1) => "w $geq 1"
    @fact conToStr(x*y ≤ 1) => "x*y - 1 $leq 0"
    @fact conToStr(x*y == 1) => "x*y - 1 $eq 0"
    @fact conToStr(x*y ≥ 1) => "x*y - 1 $geq 0"
    # 2-2 Variable--Variable
    @fact affToStr(w + x) => "w + x"
    @fact affToStr(w - x) => "w - x"
    @fact quadToStr(w * x) => "w*x"
    @fact affToStr(x - x) => "0"
    @fact_throws  w / x
    @fact conToStr(w ≤ x) => "w - x $leq 0"
    @fact conToStr(w == x) => "w - x $eq 0"
    @fact conToStr(w ≥ x) => "w - x $geq 0"
    @fact conToStr(y*z ≤ x) => "y*z - x $leq 0"
    @fact conToStr(y*z == x) => "y*z - x $eq 0"
    @fact conToStr(y*z ≥ x) => "y*z - x $geq 0"
    @fact conToStr(x ≤ x) => "0 $leq 0"
    @fact conToStr(x == x) => "0 $eq 0"
    @fact conToStr(x ≥ x) => "0 $geq 0"
    # 2-3 Variable--AffExpr
    @fact affToStr(z + aff) => "7.1 x + z + 2.5"
    @fact affToStr(z - aff) => "-7.1 x + z - 2.5"
    @fact quadToStr(z * aff) => "7.1 x*z + 2.5 z"
    @fact_throws  z / aff
    @fact conToStr(z ≤ aff) => "-7.1 x + z $leq 2.5"
    @fact conToStr(z == aff) => "-7.1 x + z $eq 2.5"
    @fact conToStr(z ≥ aff) => "-7.1 x + z $geq 2.5"
    @fact conToStr(7.1 * x - aff ≤ 0) => "0 $leq 2.5"
    @fact conToStr(7.1 * x - aff == 0) => "0 $eq 2.5"
    @fact conToStr(7.1 * x - aff ≥ 0) => "0 $geq 2.5"
    # 2-4 Variable--QuadExpr
    @fact quadToStr(w + q) => "2.5 y*z + 7.1 x + w + 2.5"
    @fact quadToStr(w - q) => "-2.5 y*z - 7.1 x + w - 2.5"
    @fact_throws  w*q
    @fact_throws  w/q
    @fact conToStr(w ≤ q) => "-2.5 y*z - 7.1 x + w - 2.5 $leq 0"
    @fact conToStr(w == q) => "-2.5 y*z - 7.1 x + w - 2.5 $eq 0"
    @fact conToStr(w ≥ q) => "-2.5 y*z - 7.1 x + w - 2.5 $geq 0"
    end

    # 3. AffExpr tests
    context("AffExpr--???") do
    # 3-0 AffExpr unary
    @fact affToStr(+aff) => "7.1 x + 2.5"
    @fact affToStr(-aff) => "-7.1 x - 2.5"
    # 3-1 AffExpr--Number
    @fact affToStr(aff + 1.5) => "7.1 x + 4"
    @fact affToStr(aff - 1.5) => "7.1 x + 1"
    @fact affToStr(aff * 2) => "14.2 x + 5"
    @fact affToStr(aff / 2) => "3.55 x + 1.25"
    @fact conToStr(aff ≤ 1) => "7.1 x $leq -1.5"
    @fact conToStr(aff == 1) => "7.1 x $eq -1.5"
    @fact conToStr(aff ≥ 1) => "7.1 x $geq -1.5"
    # 3-2 AffExpr--Variable
    @fact affToStr(aff + z) => "7.1 x + z + 2.5"
    @fact affToStr(aff - z) => "7.1 x - z + 2.5"
    @fact quadToStr(aff * z) => "7.1 x*z + 2.5 z"
    @fact_throws  aff/z
    @fact conToStr(aff ≤ z) => "7.1 x - z $leq -2.5"
    @fact conToStr(aff == z) => "7.1 x - z $eq -2.5"
    @fact conToStr(aff ≥ z) => "7.1 x - z $geq -2.5"
    @fact conToStr(aff - 7.1 * x ≤ 0) => "0 $leq -2.5"
    @fact conToStr(aff - 7.1 * x == 0) => "0 $eq -2.5"
    @fact conToStr(aff - 7.1 * x ≥ 0) => "0 $geq -2.5"
    # 3-3 AffExpr--AffExpr
    @fact affToStr(aff + aff2) => "7.1 x + 1.2 y + 3.7"
    @fact affToStr(aff - aff2) => "7.1 x - 1.2 y + 1.3"
    @fact quadToStr(aff * aff2) => "8.52 x*y + 3 y + 8.52 x + 3"
    @fact quadToStr((x+x)*(x+3)) => quadToStr((x+3)*(x+x))  # Issue #288
    @fact_throws  aff/aff2
    @fact conToStr(aff ≤ aff2) => "7.1 x - 1.2 y $leq -1.3"
    @fact conToStr(aff == aff2) => "7.1 x - 1.2 y $eq -1.3"
    @fact conToStr(aff ≥ aff2) => "7.1 x - 1.2 y $geq -1.3"
    @fact conToStr(aff-aff ≤ 0) => "0 $leq 0"
    @fact conToStr(aff-aff == 0) => "0 $eq 0"
    @fact conToStr(aff-aff ≥ 0) => "0 $geq 0"
    # 3-4 AffExpr--QuadExpr
    @fact quadToStr(aff2 + q) => "2.5 y*z + 1.2 y + 7.1 x + 3.7"
    @fact quadToStr(aff2 - q) => "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
    @fact_throws  aff2 * q
    @fact_throws  aff2 / q
    @fact conToStr(aff2 ≤ q) => "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $leq 0"
    @fact conToStr(aff2 == q) => "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $eq 0"
    @fact conToStr(aff2 ≥ q) => "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $geq 0"
    end

    # 4. QuadExpr
    context("QuadExpr--???") do
    # 4-0 QuadExpr unary
    @fact quadToStr(+q) => "2.5 y*z + 7.1 x + 2.5"
    @fact quadToStr(-q) => "-2.5 y*z - 7.1 x - 2.5"
    # 4-1 QuadExpr--Number
    @fact quadToStr(q + 1.5) => "2.5 y*z + 7.1 x + 4"
    @fact quadToStr(q - 1.5) => "2.5 y*z + 7.1 x + 1"
    @fact quadToStr(q * 2) => "5 y*z + 14.2 x + 5"
    @fact quadToStr(q / 2) => "1.25 y*z + 3.55 x + 1.25"
    @fact conToStr(q ≥ 1) => "2.5 y*z + 7.1 x + 1.5 $geq 0"
    @fact conToStr(q == 1) => "2.5 y*z + 7.1 x + 1.5 $eq 0"
    @fact conToStr(q ≤ 1) => "2.5 y*z + 7.1 x + 1.5 $leq 0"
    # 4-2 QuadExpr--Variable
    @fact quadToStr(q + w) => "2.5 y*z + 7.1 x + w + 2.5"
    @fact quadToStr(q - w) => "2.5 y*z + 7.1 x - w + 2.5"
    @fact_throws q*w
    @fact_throws q/w
    @fact conToStr(q ≤ w) => "2.5 y*z + 7.1 x - w + 2.5 $leq 0"
    @fact conToStr(q == w) => "2.5 y*z + 7.1 x - w + 2.5 $eq 0"
    @fact conToStr(q ≥ w) => "2.5 y*z + 7.1 x - w + 2.5 $geq 0"
    # 4-3 QuadExpr--AffExpr
    @fact quadToStr(q + aff2) => "2.5 y*z + 7.1 x + 1.2 y + 3.7"
    @fact quadToStr(q - aff2) => "2.5 y*z + 7.1 x - 1.2 y + 1.3"
    @fact_throws  q * aff2
    @fact_throws  q / aff2
    @fact conToStr(q ≤ aff2) => "2.5 y*z + 7.1 x - 1.2 y + 1.3 $leq 0"
    @fact conToStr(q == aff2) => "2.5 y*z + 7.1 x - 1.2 y + 1.3 $eq 0"
    @fact conToStr(q ≥ aff2) => "2.5 y*z + 7.1 x - 1.2 y + 1.3 $geq 0"
    # 4-4 QuadExpr--QuadExpr
    @fact quadToStr(q + q2) => "8 x*z + 2.5 y*z + 7.1 x + 1.2 y + 3.7"
    @fact quadToStr(q - q2) => "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3"
    @fact_throws  q * q2
    @fact_throws  q / q2
    @fact conToStr(q ≤ q2) => "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 $leq 0"
    @fact conToStr(q == q2) => "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 $eq 0"
    @fact conToStr(q ≥ q2) => "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 $geq 0"
    end
end

facts("[operator] Higher-level operators") do
context("sum") do
    sum_m = Model()
    @defVar(sum_m, 0 ≤ matrix[1:3,1:3] ≤ 1, start = 1)
    # sum(j::JuMPArray{Variable})
    @fact affToStr(sum(matrix)) => "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"
    # sum(j::JuMPArray{Variable}) in a macro
    @setObjective(sum_m, Max, sum(matrix))
    @fact quadToStr(sum_m.obj) => "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"

    # sum{T<:Real}(j::JuMPArray{T})
    @fact sum(getValue(matrix)) => roughly(9, 1e-6)
    # sum(j::Array{Variable})
    @fact affToStr(sum(matrix[1:3,1:3])) => affToStr(sum(matrix))
    # sum(affs::Array{AffExpr})
    @fact affToStr(sum([2*matrix[i,j] for i in 1:3, j in 1:3])) => "2 matrix[1,1] + 2 matrix[2,1] + 2 matrix[3,1] + 2 matrix[1,2] + 2 matrix[2,2] + 2 matrix[3,2] + 2 matrix[1,3] + 2 matrix[2,3] + 2 matrix[3,3]"

    S = [1,3]
    @defVar(sum_m, x[S], start=1)
    # sum(j::JuMPDict{Variable})
    @fact length(affToStr(sum(x))) => 11 # order depends on hashing
    @fact contains(affToStr(sum(x)),"x[1]") => true
    @fact contains(affToStr(sum(x)),"x[3]") => true
    # sum{T<:Real}(j::JuMPDict{T})
    @fact sum(getValue(x)) => 2
end

context("dot") do
    dot_m = Model()
    @defVar(dot_m, 0 ≤ x[1:3] ≤ 1)
    c = [1:3]
    @fact affToStr(dot(c,x)) => "x[1] + 2 x[2] + 3 x[3]"
    @fact affToStr(dot(x,c)) => "x[1] + 2 x[2] + 3 x[3]"

    A = [1 3 ; 2 4]
    @defVar(dot_m, 1 ≤ y[1:2,1:2] ≤ 1)
    @fact affToStr(dot(A,y)) => "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"
    @fact affToStr(dot(y,A)) => "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"

    B = ones(2,2,2)
    @defVar(dot_m, 0 ≤ z[1:2,1:2,1:2] ≤ 1)
    @fact affToStr(dot(B,z)) => "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"
    @fact affToStr(dot(z,B)) => "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"

    @setObjective(dot_m, Max, dot(x, ones(3)) - dot(y, ones(2,2)) )
    #solve(dot_m)
    for i in 1:3
        setValue(x[i], 1)
    end
    for i in 1:2, j in 1:2
        setValue(y[i,j], 1)
    end
    @fact dot(c, getValue(x)) => roughly( 6, 1e-6)
    @fact dot(A, getValue(y)) => roughly(10, 1e-6)
end
end


# The behavior in this test is no longer well-defined
#let
#    dot_m = Model()
#    @defVar(dot_m, 0 <= v[1:3,1:4] <= 1)
#    @defVar(dot_m, 0 <= w[3:2:7,7:-2:1] <= 1)
#    @fact quadToStr(dot(v,w)) => "v[1,1]*w[3,7] + v[1,2]*w[3,5] + v[1,3]*w[3,3] + v[1,4]*w[3,1] + v[2,1]*w[5,7] + v[2,2]*w[5,5] + v[2,3]*w[5,3] + v[2,4]*w[5,1] + v[3,1]*w[7,7] + v[3,2]*w[7,5] + v[3,3]*w[7,3] + v[3,4]*w[7,1]"
#
#    @addConstraint(dot_m, sum(v) + sum(w) => 2)
#    @setObjective(dot_m, :Max, v[2,3] + w[5,3])
#    solve(dot_m)
#    @fact dot(getValue(v),getValue(w)) => 1.0
#end

# Ditto
#let
#    slice_m = Model()
#    C = [:cat,:dog]
#    @defVar(slice_m, x[-1:1,C])
#    catcoef = [1,2,3]
#    dogcoef = [3,4,5]
#
#    @fact affToStr(dot(catcoef, x[:,:cat])) => "x[-1,cat] + 2 x[0,cat] + 3 x[1,cat]"
#    @fact affToStr(dot(dogcoef, x[:,:dog])) => "3 x[-1,dog] + 4 x[0,dog] + 5 x[1,dog]"
#end
