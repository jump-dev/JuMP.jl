# macros.jl
# Tests for macros
using JuMP
using Base.Test

# To ensure the tests work on Windows and Linux/OSX, we need
# to use the correct comparison operators
const leq = JuMP.repl_leq
const geq = JuMP.repl_geq
const  eq = JuMP.repl_eq

# Check for changes in Julia's expression parsing
sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j = 1:M; i != j})
#@test string(sumexpr) == "sum{\$(Expr(:parameters, :(i != j))),x[i,j] * y[i,j],i = 1:N,j = 1:M}"
@test sumexpr.head == :curly
@test length(sumexpr.args) == 5
@test sumexpr.args[1] == :sum
@test sumexpr.args[2].head == :parameters
@test sumexpr.args[3] == :(x[i,j] * y[i,j])
@test sumexpr.args[4].head == :(=)
@test sumexpr.args[5].head == :(=)


sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j = 1:M})
#@test string(sumexpr) == "sum{x[i,j] * y[i,j],i = 1:N,j = 1:M}"
@test sumexpr.head == :curly
@test length(sumexpr.args) == 4
@test sumexpr.args[1] == :sum
@test sumexpr.args[2] == :(x[i,j] * y[i,j])
@test sumexpr.args[3].head == :(=)
@test sumexpr.args[4].head == :(=)

# test JuMP's macros

let 
    m = Model()
    @defVar(m, w)
    @defVar(m, x)
    @defVar(m, y)
    @defVar(m, z)
    t = 10

    @addConstraint(m, 3x - y == 3.3(w + 2z) + 5) 
    @test conToStr(m.linconstr[end]) == "3 x - y - 3.3 w - 6.6 z $eq 5"
    if VERSION >= v"0.4.0-"
        @addConstraint(m, 3x - y == (w + 2z)*3.3 + 5)
        @test conToStr(m.linconstr[end]) == "3 x - y - 3.3 w - 6.6 z $eq 5"
    end
    @addConstraint(m, (x+y)/2 == 1) 
    @test conToStr(m.linconstr[end]) == "0.5 x + 0.5 y $eq 1"
    @addConstraint(m, -1 <= x-y <= t) 
    @test conToStr(m.linconstr[end]) == "-1 $leq x - y $leq 10"
    @addConstraint(m, -1 <= x+1 <= 1)
    @test conToStr(m.linconstr[end]) == "-2 $leq x $leq 0"
    @test_throws ErrorException @addConstraint(m, x <= t <= y)

    @defExpr(aff, 3x - y - 3.3(w + 2z) + 5)
    @test affToStr(aff) == "3 x - y - 3.3 w - 6.6 z + 5"

    if VERSION >= v"0.4.0-"
        @defExpr(qaff, (w+3)*(2x+1)+10)
        @test quadToStr(qaff) == "2 w*x + 6 x + w + 13"
    end
end

let
    m = Model()
    @defVar(m, 3.2 >= x >= 1)
    @test m.colLower == [1.0]
    @test m.colUpper == [3.2]
end

let
    m = Model()
    @defVar(m, x[1:3,1:3])
    @defVar(m, y)
    C = [1 2 3; 4 5 6; 7 8 9]
    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:2, j = 2:3 } <= 1)
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 5 x[2,2] + 6 x[2,3] $leq 1"
    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:3; i != j} == y)
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:i} == 0);
    @test conToStr(m.linconstr[end]) == "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] $eq 0"
end

let
    m = Model()
    @defVar(m, x[1:3,1:3])
    C = [1 2 3; 4 5 6; 7 8 9]

    if VERSION >= v"0.4.0-"
        @addConstraint(m, sum{ x[i,j]*(C[i,j]-1), i = 1:3, j = 1:3; i != j} == 0)
        @test conToStr(m.linconstr[end]) == "x[1,2] + 2 x[1,3] + 3 x[2,1] + 5 x[2,3] + 6 x[3,1] + 7 x[3,2] $eq 0"
    end

    con = @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:3; i != j} == 0)
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] $eq 0"

    @defVar(m, y, 0, [con], [-1.0])
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

    chgConstrRHS(con, 3)
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 3"
end

let
    m = Model()
    @defVar(m, x)
    @defVar(m, y)
    temp = x + 2y + 1
    @addConstraint(m, 3*temp - x - 2 >= 0)
    @test conToStr(m.linconstr[end]) == "6 y + 2 x $geq -1"
end

# test ranges in @defVar
let
    m = Model()
    @defVar(m, x[1:5])
    @defVar(m, y[3:2:9])
    @defVar(m, z[4:3:8])
    @defVar(m, w[6:5])

    @test length(x) == 5
    @test length(y) == 4
    @test length(z) == 2
    @test length(w) == 0

    @test x[end].col == x[5].col
    @test y[3].m == y[5].m == y[7].m == y[9].m # just make sure indexing works alright
    @test_throws KeyError z[8].col
    @test_throws BoundsError w[end]

end

# unicode comparisons
let
    m = Model()
    @defVar(m, 0 ≤ x ≤ 1)
    @defVar(m, y ≥ 2)
    @defVar(m, z ≤ 3)
    @test m.colUpper == [1.0, Inf,  3.0]
    @test m.colLower == [0.0, 2.0, -Inf]
    @addConstraint(m, 0 ≤ x + y ≤ 1)
    @addConstraint(m, x + z ≤ 2)
    @addConstraint(m, y + z ≥ 3)
    @test m.linconstr[1].lb == 0.0
    @test m.linconstr[1].ub == 1.0
    @test m.linconstr[2].lb == -Inf
    @test m.linconstr[2].ub == 2.0
    @test m.linconstr[3].lb == 3.0
    @test m.linconstr[3].ub == Inf
end

# test @addConstraint(a,b,c)
let
    m = Model()
    @defVar(m, x[1:5])
    @defVar(m, y[2:2:6])

    @addConstraint(m, c, x[4] - y[4] == 1)
    @test conToStr(m.linconstr[c.idx]) == "x[4] - y[4] $eq 1"

    @addConstraint(m, d[i=1:5,j=6:-2:2], x[i] - y[j] == 2)
    @test conToStr(m.linconstr[d[4,4].idx]) == "x[4] - y[4] $eq 2"
end

# test @addConstraints
let
    m = Model()
    @defVar(m, x)
    @defVar(m, y[1:3])

    @addConstraints m begin
        x + y[1] == 1
        ref[i=1:3], y[1] + y[i] >= i
    end

    @test conToStr(m.linconstr[1]) == "x + y[1] $eq 1"
    @test conToStr(m.linconstr[2]) == "2 y[1] $geq 1"
    @test conToStr(m.linconstr[3]) == "y[1] + y[2] $geq 2"
    @test conToStr(m.linconstr[4]) == "y[1] + y[3] $geq 3"
end

# test quadratic objective macro
let
    m = Model()
    @defVar(m, x[1:5])
    @setObjective(m, Max, sum{i*x[i]*x[j], i=1:5, j=5:-1:1; isodd(i) && iseven(j)} + 2x[5])

    @test quadToStr(m.obj) == "x[1]*x[2] + 3 x[2]*x[3] + x[1]*x[4] + 3 x[3]*x[4] + 5 x[2]*x[5] + 5 x[4]*x[5] + 2 x[5]"
end

# test quadratic constraint macro
let 
    m = Model()
    @defVar(m, x[1:5])

    @addConstraint(m, x[3]*x[1] + sum{x[i]*x[5-i+1], i=1:5; 2 <= i <= 4} + 4x[5] == 1)
    @test conToStr(m.quadconstr[end]) == "x[1]*x[3] + x[3]² + 2 x[2]*x[4] + 4 x[5] - 1 $eq 0"

    @addConstraint(m, sum{sum{(x[i] - 2)*x[j],j=4:5},i=2:3} >= -3*x[2]*2*x[4])
    @test conToStr(m.quadconstr[end]) == "7 x[2]*x[4] + x[3]*x[4] + x[2]*x[5] + x[3]*x[5] - 4 x[4] - 4 x[5] $geq 0"

    if VERSION > v"0.4.0-"
        @addConstraint(m, sum{x[i],i=1:2}*sum{x[i],i=2:3} >= 0)
        @test conToStr(m.quadconstr[end]) == "x[1]*x[2] + x[2]² + x[1]*x[3] + x[2]*x[3] $geq 0"
        @addConstraint(m, x[1]^2 + x[2]*x[3] >= 0)
        @test conToStr(m.quadconstr[end]) == "x[1]² + x[2]*x[3] + 0 $geq 0"
        @addConstraint(m, x[1]^2 + (x[2]+3)*(x[3]-1) >= 0)
        @test conToStr(m.quadconstr[end]) == "x[1]² + x[2]*x[3] + 3 x[3] - x[2] - 3 $geq 0"
    end

    myquadexpr = x[1]*x[2]
    @addConstraint(m, sum{i*myquadexpr + x[i], i=1:3} + sum{x[i] + myquadexpr*i, i=1:3} == 0)
    @test conToStr(m.quadconstr[end]) == "12 x[1]*x[2] + 2 x[1] + 2 x[2] + 2 x[3] $eq 0"
end

# Test "triangular indexing"
n = 10
trimod = Model()
@defVar(trimod, x[i=1:n,j=i:n])
@defVar(trimod, y[i=3:2:7,j=-i])
@test getNumVars(trimod) == n*(n+1)/2 + 3
S = {(i,i+2) for i in 1:5}
@defVar(trimod, z[(i,j)=S,k=i:j])
@test length(z.tupledict) == 15
@addConstraint(trimod, cref[i=1:n,j=i:n], x[i,j] + y[5,-5] == 1)
@test getNumConstraints(trimod) == n*(n+1)/2

# Test iteration over JuMPDicts
cntr = zeros(Bool, n, n)
for (i,j,var) in x
    @test isequal(x[i,j], var)
    cntr[i,j] = true
end
for i in 1:n, j in 1:n
    if j >= i
        @test cntr[i,j]
    else
        @test !cntr[i,j]
    end
end

# test @defExpr
let
    model = Model()
    @defVar(model, x[1:3,1:3])
    @defExpr(expr, sum{i*x[i,j] + j, i=1:3,j=1:3})
    @test affToStr(expr) == "x[1,1] + x[1,2] + x[1,3] + 2 x[2,1] + 2 x[2,2] + 2 x[2,3] + 3 x[3,1] + 3 x[3,2] + 3 x[3,3] + 18"

    @test_throws ErrorException @defExpr(blah[i=1:3], x[i,1]^2)
end
