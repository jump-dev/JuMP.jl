# macros.jl
# Tests for macros

# Check for changes in Julia's expression parsing
sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j = 1:M; i != j})
if VERSION >= v"0.3.0-"
    @test string(sumexpr) == "sum{\$(Expr(:parameters, :(i != j))),x[i,j] * y[i,j],i = 1:N,j = 1:M}"
end
@test sumexpr.head == :curly
@test length(sumexpr.args) == 5
@test sumexpr.args[1] == :sum
@test sumexpr.args[2].head == :parameters
@test sumexpr.args[3] == :(x[i,j] * y[i,j])
@test sumexpr.args[4].head == :(=)
@test sumexpr.args[5].head == :(=)


sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j = 1:M})
if VERSION >= v"0.3.0-"
    @test string(sumexpr) == "sum{x[i,j] * y[i,j],i = 1:N,j = 1:M}"
end
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
    @test conToStr(m.linconstr[end]) == "3 x - y - 3.3 w - 6.6 z == 5"
    @addConstraint(m, (x+y)/2 == 1) 
    @test conToStr(m.linconstr[end]) == "0.5 x + 0.5 y == 1"
    @addConstraint(m, -1 <= x-y <= t) 
    @test conToStr(m.linconstr[end]) == "-1 <= x - y <= 10"
    @addConstraint(m, -1 <= x+1 <= 1)
    @test conToStr(m.linconstr[end]) == "-2 <= x <= 0"
    @test_throws @addConstraint(m, x <= t <= y)
end

let
    m = Model()
    @defVar(m, x[1:3,1:3])
    @defVar(m, y)
    C = [1 2 3; 4 5 6; 7 8 9]
    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:2, j = 2:3 } <= 1)
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 5 x[2,2] + 6 x[2,3] <= 1"
    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:3; i != j} == y)
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y == 0"

    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:i} == 0);
    @test conToStr(m.linconstr[end]) == "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] == 0"
end

let
    m = Model()
    @defVar(m, x[1:3,1:3])
    C = [1 2 3; 4 5 6; 7 8 9]
    con = @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:3; i != j} == 0)
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] == 0"

    @defVar(m, y, 0, [con], [-1.0])
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y == 0"

    chgConstrRHS(con, 3)
    @test conToStr(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y == 3"
end

let
    m = Model()
    @defVar(m, x)
    @defVar(m, y)
    temp = x + 2y + 1
    @addConstraint(m, 3*temp - x - 2 >= 0)
    @test conToStr(m.linconstr[end]) == "6 y + 2 x >= -1"
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
    @test_throws z[8].col
    @test_throws w[end]

end

# unicode comparisons
if VERSION > v"0.3.0-"
    let
        eval(parse("m = Model(); @defVar(m, 0 ≤ x ≤ 1); @defVar(m, y ≥ 2); @defVar(m, z ≤ 3); @test m.colUpper == [1.0, Inf,  3.0]; @test m.colLower == [0.0, 2.0, -Inf]; @addConstraint(m, 0 ≤ x + y ≤ 1); @addConstraint(m, x + z ≤ 2); @addConstraint(m, y + z ≥ 3); @test m.linconstr[1].lb == 0.0; @test m.linconstr[1].ub == 1.0; @test m.linconstr[2].lb == -Inf; @test m.linconstr[2].ub == 2.0; @test m.linconstr[3].lb == 3.0; @test m.linconstr[3].ub == Inf"))
    end
end

# test @addConstraint(a,b,c)
let
    m = Model()
    @defVar(m, x[1:5])
    @defVar(m, y[2:2:6])

    @addConstraint(m, c, x[4] - y[4] == 1)
    @test conToStr(m.linconstr[c.idx]) == "x[4] - y[4] == 1"

    @addConstraint(m, d[i=1:5,j=6:-2:2], x[i] - y[j] == 2)
    @test conToStr(m.linconstr[d[4,4].idx]) == "x[4] - y[4] == 2"
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

    @test conToStr(m.linconstr[1]) == "x + y[1] == 1"
    @test conToStr(m.linconstr[2]) == "2 y[1] >= 1"
    @test conToStr(m.linconstr[3]) == "y[1] + y[2] >= 2"
    @test conToStr(m.linconstr[4]) == "y[1] + y[3] >= 3"
end
