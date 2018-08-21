# test/macros.jl
# Testing macros work correctly
#############################################################################
using JuMP, Compat
using Compat.Test
import MathProgBase
# To ensure the tests work on Windows and Linux/OSX, we need
# to use the correct comparison operators
const leq = JuMP.repl[:leq]
const geq = JuMP.repl[:geq]
const  eq = JuMP.repl[:eq]
const Vert = JuMP.repl[:Vert]
const sub2 = JuMP.repl[:sub2]

struct __Cone__ end

mutable struct MyVariable
    lowerbound
    upperbound
    category
    basename::String
    start
    test_kw::Int
end

@testset "Macros" begin


    @testset "Check Julia curly expression parsing" begin
        sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j in 1:M; i != j})
        @test sumexpr.head == :curly
        @test length(sumexpr.args) == 5
        @test sumexpr.args[1] == :sum
        @test sumexpr.args[2].head == :parameters
        @test sumexpr.args[3] == :(x[i,j] * y[i,j])
        @test sumexpr.args[4].head == :(=)
        @test sumexpr.args[5].head == :call
        @test sumexpr.args[5].args[1] == :in

        sumexpr = :(sum{x[i,j] * y[i,j], i in 1:N, j = 1:M})
        @test sumexpr.head == :curly
        @test length(sumexpr.args) == 4
        @test sumexpr.args[1] == :sum
        @test sumexpr.args[2] == :(x[i,j] * y[i,j])
        @test sumexpr.args[3].head == :call
        @test sumexpr.args[3].args[1] == :in

        @test sumexpr.args[4].head == :(=)
    end

    @testset "Check Julia generator expression parsing" begin
        sumexpr = :(sum(x[i,j] * y[i,j] for i = 1:N, j in 1:M if i != j))
        @test sumexpr.head == :call
        @test sumexpr.args[1] == :sum
        @test sumexpr.args[2].head == :generator
        @test sumexpr.args[2].args[1] == :(x[i,j] * y[i,j])
        @test sumexpr.args[2].args[2].head == :filter
        @test sumexpr.args[2].args[2].args[1] == :(i != j)
        @test sumexpr.args[2].args[2].args[2] == :(i = 1:N)
        @test sumexpr.args[2].args[2].args[3] == :(j = 1:M)

        sumexpr = :(sum(x[i,j] * y[i,j] for i = 1:N, j in 1:M))
        @test sumexpr.head == :call
        @test sumexpr.args[1] == :sum
        @test sumexpr.args[2].head == :generator
        @test sumexpr.args[2].args[1] == :(x[i,j] * y[i,j])
        @test sumexpr.args[2].args[2] == :(i = 1:N)
        @test sumexpr.args[2].args[3] == :(j = 1:M)
    end

    @testset "Check @constraint basics" begin
        m = Model()
        @variable(m, w)
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        t = 10

        @constraint(m, 3x - y == 3.3(w + 2z) + 5)
        @test string(m.linconstr[end]) == "3 x - y - 3.3 w - 6.6 z $eq 5"
        @constraint(m, 3x - y == (w + 2z)*3.3 + 5)
        @test string(m.linconstr[end]) == "3 x - y - 3.3 w - 6.6 z $eq 5"
        @constraint(m, (x+y)/2 == 1)
        @test string(m.linconstr[end]) == "0.5 x + 0.5 y $eq 1"
        @constraint(m, -1 <= x-y <= t)
        @test string(m.linconstr[end]) == "-1 $leq x - y $leq 10"
        @constraint(m, -1 <= x+1 <= 1)
        @test string(m.linconstr[end]) == "-2 $leq x $leq 0"
        @constraint(m, -1 <= x <= 1)
        @test string(m.linconstr[end]) == "-1 $leq x $leq 1"
        @constraint(m, -1 <= x <= sum(0.5 for i = 1:2))
        @test string(m.linconstr[end]) == "-1 $leq x $leq 1"
        @test_throws ErrorException @constraint(m, x <= t <= y)
        if VERSION < v"0.7"
            @test macroexpand(:(@constraint(m, 1 >= x >= 0))).head == :error
            @test macroexpand(:(@constraint(1 <= x <= 2, foo=:bar))).head == :error
        else
            @test_throws LoadError @macroexpand @constraint(m, 1 ≥ x ≥ 0)
            @test_throws LoadError @macroexpand @constraint(1 ≤ x ≤ 2, foo=:bar)
        end

        @expression(m, aff, 3x - y - 3.3(w + 2z) + 5)
        @test string(aff) == "3 x - y - 3.3 w - 6.6 z + 5"

        @constraint(m, 3 + 5*7 <= 0)
        @test string(m.linconstr[end]) == "0 $leq -38"

        @expression(m, qaff, (w+3)*(2x+1)+10)
        @test string(qaff) == "2 w*x + 6 x + w + 13"
    end

    @testset "Checking @variable with reverse direction bounds" begin
        m = Model()
        @variable(m, 3.2 >= x >= 1)
        @test m.colLower == [1.0]
        @test m.colUpper == [3.2]
    end

    @testset "sum{} (deprecated)" begin
        m = Model()
        @variable(m, x[1:3,1:3])
        @variable(m, y)
        C = [1 2 3; 4 5 6; 7 8 9]
        @constraint(m, sum{ C[i,j]*x[i,j], i in 1:2, j = 2:3 } <= 1)
        @test string(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 5 x[2,2] + 6 x[2,3] $leq 1"
        @constraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j in 1:3; i != j} == y)
        @test string(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

        @constraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:i} == 0);
        @test string(m.linconstr[end]) == "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] $eq 0"

        @constraint(m, sum{ 0*x[i,1], i=1:3} == 0)
        @test string(m.linconstr[end]) == "0 $eq 0"

        @constraint(m, sum{ 0*x[i,1] + y, i=1:3} == 0)
        @test string(m.linconstr[end]) == "3 y $eq 0"

    end

    @testset "[macros] sum(generator)" begin
        m = Model()
        @variable(m, x[1:3,1:3])
        @variable(m, y)
        C = [1 2 3; 4 5 6; 7 8 9]
        @constraint(m, sum( C[i,j]*x[i,j] for i in 1:2, j = 2:3 ) <= 1)
        @test string(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 5 x[2,2] + 6 x[2,3] $leq 1"
        @constraint(m, sum( C[i,j]*x[i,j] for i = 1:3, j in 1:3 if i != j) == y)
        @test string(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

        @constraint(m, sum( C[i,j]*x[i,j] for i = 1:3, j = 1:i) == 0);
        @test string(m.linconstr[end]) == "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] $eq 0"

        @constraint(m, sum( C[i,j]*x[i,j] for i = 1:3 for j = 1:i) == 0);
        @test string(m.linconstr[end]) == "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] $eq 0"

        @constraint(m, sum( C[i,j]*x[i,j] for i = 1:3 if true for j = 1:i) == 0);
        @test string(m.linconstr[end]) == "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] $eq 0"

        @constraint(m, sum( C[i,j]*x[i,j] for i = 1:3 if true for j = 1:i if true) == 0);
        @test string(m.linconstr[end]) == "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] $eq 0"

        @constraint(m, sum( 0*x[i,1] for i=1:3) == 0)
        @test string(m.linconstr[end]) == "0 $eq 0"

        @constraint(m, sum( 0*x[i,1] + y for i=1:3) == 0)
        @test string(m.linconstr[end]) == "3 y $eq 0"

        #@test isexpr(macroexpand(:(@constraint(m, sum( 0*x[i,1] + y for i=1:3 for j in 1:3) == 0))),:error) == true

    end

    @testset "Problem modification" begin
        m = Model()
        @variable(m, x[1:3,1:3])
        C = [1 2 3; 4 5 6; 7 8 9]

        @constraint(m, sum( x[i,j]*(C[i,j]-1) for i in 1:3, j = 1:3 if i != j) == 0)
        @test string(m.linconstr[end]) == "x[1,2] + 2 x[1,3] + 3 x[2,1] + 5 x[2,3] + 6 x[3,1] + 7 x[3,2] $eq 0"

        con = @constraint(m, sum( C[i,j]*x[i,j] for i = 1:3, j = 1:3 if i != j) == 0)
        @test string(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] $eq 0"

        @variable(m, y, objective = 0, inconstraints = [con], coefficients = [-1.0])
        @test string(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

        JuMP.setRHS(con, 3)
        @test string(m.linconstr[end]) == "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 3"
    end

    @testset "Using pre-built affine is OK in macro" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        temp = x + 2y + 1
        @constraint(m, 3*temp - x - 2 >= 0)
        @test string(m.linconstr[end]) == "6 y + 2 x $geq -1"
        # More complex expression
        a = 1.0*x
        @constraint(m, (2+2)*((3+4)*(1+a)) == 0)
        @test string(m.linconstr[end]) == "28 x $eq -28"
        @test string(a) == "x"

        @test string(@LinearConstraint(1 + 0*temp == 0)) == "0 $eq -1"
    end

    @testset "Test ranges in @variable" begin
        m = Model()
        @variable(m, x[1:5])
        @variable(m, y[3:2:9])
        @variable(m, z[4:3:8])
        @variable(m, w[6:5])

        @test length(x) == 5
        @test length(y) == 4
        @test length(z) == 2
        @test length(w) == 0

        @test x[end].col == x[5].col
        @test y[3].m == y[5].m
        @test y[3].m == y[7].m
        @test y[3].m == y[9].m
        @test_throws ErrorException z[8].col
        @test_throws ErrorException w[end]
    end

    @testset "Unicode comparisons" begin
        m = Model()
        @variable(m, 0 ≤ x ≤ 1)
        @variable(m, y ≥ 2)
        @variable(m, z ≤ 3)
        @test m.colUpper == [1.0, Inf,  3.0]
        @test m.colLower == [0.0, 2.0, -Inf]
        @constraint(m, 0 ≤ x + y ≤ 1)
        @constraint(m, x + z ≤ 2)
        @constraint(m, y + z ≥ 3)
        @constraint(m, y*z ≤ 1)
        @test m.linconstr[1].lb == 0.0
        @test m.linconstr[1].ub == 1.0
        @test m.linconstr[2].lb == -Inf
        @test m.linconstr[2].ub == 2.0
        @test m.linconstr[3].lb == 3.0
        @test m.linconstr[3].ub == Inf
        @test m.quadconstr[1].sense == :(<=)
    end

    @testset "Three argument @constraint" begin
        m = Model()
        @variable(m, x[1:5])
        @variable(m, y[2:2:6])

        @constraint(m, c, x[4] - y[4] == 1)
        @test string(m.linconstr[c.idx]) == "x[4] - y[4] $eq 1"

        @constraint(m, d[i in 1:5,j=6:-2:2], x[i] - y[j] == 2)
        @test string(m.linconstr[d[4,4].idx]) == "x[4] - y[4] $eq 2"

        @constraint(m, q[i=1:5], x[i]^2 == 1)
        @test string(m.quadconstr[q[5].idx]) == "x[5]² - 1 $eq 0"
    end

    @testset "@constraints" begin
        m = Model()
        @variable(m, x)
        @variable(m, y[1:3])

        @constraints(m, begin
            x + y[1] == 1
            ref[i=1:3], y[1] + y[i] >= i
        end)

        @test string(m.linconstr[1]) == "x + y[1] $eq 1"
        @test string(m.linconstr[2]) == "2 y[1] $geq 1"
        @test string(m.linconstr[3]) == "y[1] + y[2] $geq 2"
        @test string(m.linconstr[4]) == "y[1] + y[3] $geq 3"
    end

    @testset "@NLconstraints" begin
        m = Model()
        @variable(m, 0 <= x <= 1)
        @variable(m, y[1:3])
        @objective(m, Max, x)

        @NLconstraints(m, begin
            ref[i=1:3], y[i] == 0
            x + y[1] * y[2] * y[3] <= 0.5
        end)

        @test length(m.nlpdata.nlconstr) == 4
        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:ExprGraph])

        @test MathProgBase.constr_expr(d,1) == :(x[2] - 0.0 == 0.0)
        @test MathProgBase.constr_expr(d,2) == :(x[3] - 0.0 == 0.0)
        @test MathProgBase.constr_expr(d,3) == :(x[4] - 0.0 == 0.0)
        @test MathProgBase.constr_expr(d,4) == :((x[1] + x[2] * x[3] * x[4]) - 0.5 <= 0.0)

    end

    @testset "Vectors in nonlinear expressions" begin
        m = Model()
        @variable(m, x[1:3])
        @test_throws ErrorException @NLobjective(m, Min, x)
        @test_throws ErrorException @NLobjective(m, Min, [1,2,3])
    end

    @testset "@objective with quadratic" begin
        m = Model()
        @variable(m, x[1:5])
        @objective(m, Max, sum(i*x[i]*x[j] for i=1:5, j=5:-1:1 if isodd(i) && iseven(j)) + 2x[5])

        @test string(m.obj) == "x[1]*x[2] + 3 x[2]*x[3] + x[1]*x[4] + 3 x[3]*x[4] + 5 x[2]*x[5] + 5 x[4]*x[5] + 2 x[5]"
    end

    @testset "@constraint with quadratic" begin
        m = Model()
        @variable(m, x[1:5])

        @constraint(m, x[3]*x[1] + sum(x[i]*x[5-i+1] for i=1:5 if 2 <= i <= 4) + 4x[5] == 1)
        @test string(m.quadconstr[end]) == "x[1]*x[3] + x[3]² + 2 x[2]*x[4] + 4 x[5] - 1 $eq 0"

        @constraint(m, sum(sum((x[i] - 2)*x[j] for j=4:5) for i=2:3) >= -3*x[2]*2*x[4])
        @test string(m.quadconstr[end]) == "7 x[2]*x[4] + x[3]*x[4] + x[2]*x[5] + x[3]*x[5] - 4 x[4] - 4 x[5] $geq 0"

        foo(x) = x
        @constraint(m, x[1] ≤ foo(x[1])^2)
        @test string(m.quadconstr[end]) == "-x[1]² + x[1] $leq 0"

        @constraint(m, sum(x[i] for i=1:2)*sum(x[i] for i=2:3) >= 0)
        @test string(m.quadconstr[end]) == "x[1]*x[2] + x[2]² + x[1]*x[3] + x[2]*x[3] $geq 0"
        @constraint(m, x[1]^2 + x[2]*x[3] >= 0)
        @test string(m.quadconstr[end]) == "x[1]² + x[2]*x[3] $geq 0"
        @constraint(m, x[1]^2 + (x[2]+3)*(x[3]-1) >= 0)
        @test string(m.quadconstr[end]) == "x[1]² + x[2]*x[3] + 3 x[3] - x[2] - 3 $geq 0"
        @constraint(m, sum(x[i] for i=1:2)^2 >= 0)
        @test string(m.quadconstr[end]) == "x[1]² + 2 x[1]*x[2] + x[2]² $geq 0"

        myquadexpr = x[1]*x[2]
        @constraint(m, sum(i*myquadexpr + x[i] for i=1:3) + sum(x[i] + myquadexpr*i for i=1:3) == 0)
        @test string(m.quadconstr[end]) == "12 x[1]*x[2] + 2 x[1] + 2 x[2] + 2 x[3] $eq 0"

        @constraint(m, (x[1] + x[2])*sum( 0*x[i] + x[3] for i=1:3) == 0)
        @test string(m.quadconstr[end]) == "3 x[1]*x[3] + 3 x[2]*x[3] $eq 0"

        @test string(@QuadConstraint(1 + 0*myquadexpr == 0)) == "1 $eq 0"

        @variable(m, y)
        @test string(@QuadConstraint(1 + (2y)*y   == 0)) == "2 y² + 1 $eq 0"
        @test string(@QuadConstraint(1 +   y *y*2 == 0)) == "2 y² + 1 $eq 0"
        z = 2y
        @test string(@QuadConstraint(y*y + y*z == 0)) == "3 y² $eq 0"
    end

    @testset "Triangular indexing, iteration" begin
        n = 10
        trimod = Model()
        @variable(trimod, x[i=1:n,j=i:n])
        @variable(trimod, y[i=3:2:7,j=-i])
        @test MathProgBase.numvar(trimod) == n*(n+1)/2 + 3
        S = Any[(i,i+2) for i in 1:5]
        @variable(trimod, z[(i,j)=S,k=i:j])
        @test length(z.tupledict) == 15
        @constraint(trimod, cref[i=1:n,j=i:n], x[i,j] + y[5,-5] == 1)
        @test MathProgBase.numconstr(trimod) == n*(n+1)/2

        cntr = zeros(Bool, n, n)
        for ((i,j),var) in zip(keys(x),values(x))
            @test x[i,j] === var
            cntr[i,j] = true
        end
        for i in 1:n, j in 1:n
            @test cntr[i,j] == (j >= i)
        end
    end

    @testset "Multidimensional indexing" begin
        model = Model()
        I1 = 1:5
        I2 = 2:8
        I3 = 5:6
        @variable(model, x[1:5,2:8,5:6])
        coll = Int[]
        for v in values(x)
            push!(coll, v.col)
        end
        p = 1
        match = true
        for v in x.innerArray
            match &= (coll[p] == v.col)
            p += 1
        end
        @test match
    end

    @testset "@expression" begin
        model = Model()
        @variable(model, x[1:3,1:3])
        @expression(model, expr, sum(i*x[i,j] + j for i=1:3,j in 1:3))
        @test string(expr) == "x[1,1] + x[1,2] + x[1,3] + 2 x[2,1] + 2 x[2,2] + 2 x[2,3] + 3 x[3,1] + 3 x[3,2] + 3 x[3,3] + 18"

        @test_throws ErrorException @expression(model, blah[i=1:3], x[i,1]^2)

        @expression(model, y[i=1:2], sum(x[i,1] for _ in 1 if i == 1))
        @test string(y[1]) == "x[1,1]"
        @test string(y[2]) == "0"

        t = @expression(model, sum(x[i,j] for i in 1:3 for j in 1:3 if isodd(i+j)))
        @test string(t) == "x[1,2] + x[2,1] + x[2,3] + x[3,2]"
    end

    @testset "Conditions in constraint indexing" begin
        model = Model()
        @variable(model, x[1:10])
        @constraint(model, c1[i=1:9;isodd(i)], x[i] + x[i+1] <= 1)
        @NLconstraint(model, c2[i=["red","blue","green"], k=9:-2:2; (i == "red" && isodd(k)) || (k >=4 && (i == "blue" || i == "green"))], x[k]^3 <= 1)
        @test length(model.linconstr) == 5
        @test string(model.linconstr[1]) == "x[1] + x[2] $leq 1"
        @test string(model.linconstr[2]) == "x[3] + x[4] $leq 1"
        @test string(model.linconstr[3]) == "x[5] + x[6] $leq 1"
        @test string(model.linconstr[4]) == "x[7] + x[8] $leq 1"
        @test string(model.linconstr[5]) == "x[9] + x[10] $leq 1"
        @test length(model.nlpdata.nlconstr) == 10
    end

    @testset "Test changes in condition parsing" begin
        ex = :(x[12;3])
        @test ex.head == :typed_vcat
        @test ex.args == [:x, 12, 3]

        ex = :(x[i=1:3,j=S;isodd(i) && i+j>=2])
        if VERSION ≤ v"0.7-"
            @test ex.head == :typed_vcat
            @test ex.args == [:x,
                              Expr(:parameters, Expr(:&&, :(isodd(i)), :(i+j>=2))),
                              Expr(:(=), :i, :(1:3)),
                              Expr(:(=), :j, :S)]
        elseif VERSION < v"1.0-"
            @test ex.head == :ref
            @test ex.args == [:x,
                              Expr(:parameters, Expr(:&&, :(isodd(i)), :(i+j>=2))),
                              Expr(:(=), :i, :(1:3)),
                              Expr(:(=), :j, :S)]
        else
            @test ex.head == :ref
            @test ex.args == [:x,
                              Expr(:parameters, Expr(:&&, :(isodd(i)), :(i+j>=2))),
                              Expr(:kw, :i, :(1:3)),
                              Expr(:kw, :j, :S)]
        end
    end

    @testset "Curly norm parsing (deprecated)" begin
        model = Model()
        @variable(model, x[1:2,1:2])
        @constraint(model, -2norm2{x[i,j], i in 1:2, j=1:2} + x[1,2] >= -1)
        @constraint(model, -2norm2{x[i,j], i=1:2, j in 1:2; iseven(i+j)} + x[1,2] >= -1)
        @constraint(model, 1 >= 2*norm2{x[i,1], i in 1:2})
        @test string(model.socconstr[1]) == "2.0 $Vert[x[1,1],x[1,2],x[2,1],x[2,2]]$Vert$sub2 $leq x[1,2] + 1"
        @test string(model.socconstr[2]) == "2.0 $Vert[x[1,1],x[2,2]]$Vert$sub2 $leq x[1,2] + 1"
        @test string(model.socconstr[3]) == "2.0 $Vert[x[1,1],x[2,1]]$Vert$sub2 $leq 1"
        @test_throws MethodError @constraint(model, (x[1,1]+1)*norm2{x[i,j], i=1:2, j=1:2} + x[1,2] >= -1)
        @test_throws ErrorException @constraint(model, norm2{x[i,j], i=1:2, j=1:2} + x[1,2] >= -1)
    end

    @testset "Generator norm parsing" begin
        model = Model()
        @variable(model, x[1:2,1:2])
        @constraint(model, -2norm(x[i,j] for i in 1:2, j=1:2) + x[1,2] >= -1)
        @constraint(model, -2norm(x[i,j] for i=1:2, j in 1:2 if iseven(i+j)) + x[1,2] >= -1)
        @constraint(model, 1 >= 2*norm(x[i,1] for i in 1:2))
        @test string(model.socconstr[1]) == "2.0 $Vert[x[1,1],x[1,2],x[2,1],x[2,2]]$Vert$sub2 $leq x[1,2] + 1"
        @test string(model.socconstr[2]) == "2.0 $Vert[x[1,1],x[2,2]]$Vert$sub2 $leq x[1,2] + 1"
        @test string(model.socconstr[3]) == "2.0 $Vert[x[1,1],x[2,1]]$Vert$sub2 $leq 1"
        @test_throws MethodError @constraint(model, (x[1,1]+1)*norm(x[i,j] for i=1:2, j=1:2) + x[1,2] >= -1)
        @test_throws ErrorException @constraint(model, norm(x[i,j] for i=1:2, j=1:2) + x[1,2] >= -1)
    end

    @testset "Extraneous terms in QuadExpr (#535)" begin
        model = Model()
        @variable(model, x)
        @variable(model, y)
        @constraint(model, x*x <= y*y)
        @test string(model.quadconstr[1]) == "x² - y² $leq 0"
    end

    @testset "Special-case binary multiplication in addtoexpr_reorder (#537)" begin
        dual = Model()
        @variable(dual, α[1:3] <= 0)
        @variable(dual, γ >= 0)
        @constraint(dual, ones(3)*γ .<= α)
        @test string(dual.linconstr[1]) == "γ - α[1] $leq 0"
        @test string(dual.linconstr[2]) == "γ - α[2] $leq 0"
        @test string(dual.linconstr[3]) == "γ - α[3] $leq 0"
    end

    @testset "Indices in macros don't leak out of scope (#582)" begin
        m = Model()
        cnt = 4
        for i in 5:8
            x = @variable(m, [i=1:3,j=1:3], upperbound = i)
            cnt += 1
            @test i == cnt
        end
        cnt = 4
        for i in 5:8
            y = @variable(m, [i=2:4,j=1:3], upperbound = i)
            cnt += 1
            @test i == cnt
        end
        cnt = 4
        for i in 5:8
            z = @variable(m, [i=[1:3;],j=1:3], upperbound = i)
            cnt += 1
            @test i == cnt
        end
        @test m.colUpper == vcat(repeat([1.0,2.0,3.0], inner=[3], outer=[4]),
                                  repeat([2.0,3.0,4.0], inner=[3], outer=[4]),
                                  repeat([1.0,2.0,3.0], inner=[3], outer=[4]))
        @variable(m, x[i=1:3] ≤ i)
        cnt = 4
        for i in 5:8
            @constraint(m, sum(x[i] for i=1, j=1, k=1) == 1)
            cnt += 1
            @test i == cnt
        end
        cnt = 4
        for i in 5:8
            @constraint(m, norm(x[i] for i=1, j=1, k=1) <= 1)
            cnt += 1
            @test i == cnt
        end
        cnt = 4
        for i in 5:8
            @constraint(m, [i=1:3,j=1:3], x[i] == 1)
            cnt += 1
            @test i == cnt
        end
    end

    @testset "Issue #621" begin
        m = Model()
        @variable(m, x)

        q = x^2
        a = x+1
        con = @QuadConstraint(a+q*3 <= 0)
        @test string(con) == "3 x² + x + 1 $leq 0"
    end

    @testset "@variables and @constraints" begin
        m = Model()
        @variables m begin
            0 ≤ x[i=1:2] ≤ i
            y ≥ 2, Int, (start = 0.7)
            z ≤ 3, (start=10)
            q, (Bin, start=0.5)
        end
        @constraints m begin
            0 ≤ x[1] + y ≤ 1
            x[1] + z ≤ 2
            y + z ≥ 3
            y*z ≤ 1
        end
        @test m.colUpper == [1.0, 2.0, Inf,  3.0, 1.0]
        @test m.colLower == [0.0, 0.0, 2.0, -Inf, 0.0]
        @test m.colCat == [:Cont, :Cont, :Int, :Cont, :Bin]
        @test getvalue(y) == 0.7
        @test getvalue(z) == 10
        @test getvalue(q) == 0.5
        @test m.linconstr[1].lb == 0.0
        @test m.linconstr[1].ub == 1.0
        @test m.linconstr[2].lb == -Inf
        @test m.linconstr[2].ub == 2.0
        @test m.linconstr[3].lb == 3.0
        @test m.linconstr[3].ub == Inf
        @test m.quadconstr[1].sense == :(<=)
    end

    @testset "@expressions and @NLexpressions" begin
        m = Model()
        @variable(m, x)

        @expressions(m, begin
            myex[i=1:2], x + i
            myex2,       x + 3
        end)
        setvalue(x, 1)
        @test getvalue(myex[1]) == 2
        @test getvalue(myex[2]) == 3
        @test getvalue(myex2)   == 4
        @test_throws BoundsError getvalue(myex[3])

        @NLexpressions(m, begin
            nlex,         x^2
            nlex2[i=1:2], x^i + sqrt(x)
        end)
        setvalue(x, 2)
        @test getvalue(nlex) == 4
        @test getvalue(nlex2[1]) == 2^1 + sqrt(2)
        @test getvalue(nlex2[2]) == 2^2 + sqrt(2)
        @test_throws BoundsError getvalue(nlex2[3])
    end

    @testset "No bare symbols in constraint macros" begin
        m = Model()
        @variable(m, x)
        if VERSION < v"0.7-"
            @test macroexpand(:(@constraint(m, x))).head == :error
            @test macroexpand(:(@constraint(m, :foo))).head == :error
            @test macroexpand(:(@SDconstraint(m, x))).head == :error
            @test macroexpand(:(@SDconstraint(m, :foo))).head == :error
            @test macroexpand(:(@NLconstraint(m, x))).head == :error
            @test macroexpand(:(@NLconstraint(m, :foo))).head == :error
        else
            @test_throws LoadError @macroexpand @constraint(m, x)
            @test_throws LoadError @macroexpand @constraint(m, :foo)
            @test_throws LoadError @macroexpand @SDconstraint(m, x)
            @test_throws LoadError @macroexpand @SDconstraint(m, :foo)
            @test_throws LoadError @macroexpand @NLconstraint(m, x)
            @test_throws LoadError @macroexpand @NLconstraint(m, :foo)
        end
    end

    @testset "LB/UB kwargs" begin
        m = Model()
        @variable(m, a, lowerbound=0)
        @variable(m, b, lowerbound=0, upperbound=1)
        @variable(m, c, upperbound=1, lowerbound=0)
        @variable(m, d, upperbound=1)
        @variable(m, e >= 0, upperbound=1)
        @variable(m, f <= 1, lowerbound=0)
        @test m.colLower == [0,0,0,-Inf,0,0]
        @test m.colUpper == [Inf,1,1,1,1,1]

        if VERSION < v"0.7-"
            @test macroexpand(:(@variable(m, g >= 0, lowerbound=1))).head == :error
            @test macroexpand(:(@variable(m, h <= 1, upperbound=1))).head == :error
            @test macroexpand(:(@variable(m, 0 <= i <= 1, lowerbound=1))).head == :error
            @test macroexpand(:(@variable(m, 0 <= j <= 1, upperbound=1))).head == :error
        else
            @test_throws LoadError @macroexpand @variable(m, g >= 0, lowerbound=1)
            @test_throws LoadError @macroexpand @variable(m, h <= 1, upperbound=1)
            @test_throws LoadError @macroexpand @variable(m, 0 <= i <= 1, lowerbound=1)
            @test_throws LoadError @macroexpand @variable(m, 0 <= j <= 1, upperbound=1)
        end
    end

    @testset "Anonymous versions of macros" begin
        m = Model()
        x = @variable(m, [1:3], lowerbound=1.0)
        y = @variable(m, [2:3,1:3], upperbound=1.0, lowerbound=0.0, Bin)
        cat = Dict(:red => :Int, :blue => :Bin)
        z = @variable(m, [s=[:red,:blue]], category = cat[s])
        w = @variable(m, [i=1:4,j=1:4;isodd(i+j)], SemiCont)
        # v = @variable(m, [i=1:3,j=1:3], Symmetric, lowerbound = eye(3)[i,j])
        u = @variable(m, [1:4,1:4], SDP)
        if VERSION < v"0.7-"
            @test macroexpand(:(@variable(m, [1:3] <= 1))).head == :error
        else
            @test_throws LoadError @macroexpand @variable(m, [1:3] <= 1)
        end

        @test getlowerbound(x[1]) == 1.0
        @test getupperbound(x[1]) == Inf
        @test getcategory(x[1]) == :Cont
        @test getlowerbound(y[2,1]) == 0.0
        @test getupperbound(y[2,1]) == 1.0
        @test getcategory(y[2,1]) == :Bin
        @test getlowerbound(z[:red]) == -Inf
        @test getupperbound(z[:red]) == Inf
        @test getcategory(z[:red]) == :Int
        @test getcategory(z[:blue]) == :Bin
        @test getlowerbound(w[1,2]) == -Inf
        @test getupperbound(w[1,2]) == Inf
        @test getcategory(w[1,2]) == :SemiCont
        # @test getlowerbound(v[1,2]) == 0.0
        # @test getupperbound(v[1,2]) == Inf
        # @test getcategory(v[1,2]) == :Cont
        @test getlowerbound(u[1,2]) == -Inf
        @test getupperbound(u[1,2]) == Inf
        @test getcategory(u[1,2]) == :Cont

        c = @constraint(m, [i=1:3], x[i] <= z[:red])
        d = @NLconstraint(m, [i=1:3], x[i]^3 == 1)
        e = @NLexpression(m, [i=2:3], y[2,i]^3)
        f = @expression(m, [i=[:red,:blue]], u[1,2] + 2z[i])

        # not sure how else to test this
        @test c[1] == c[1]
        @test d[1] == d[1]
        @test e[2] == e[2]
        @test f[:red] == f[:red]
    end

    @testset "Colons in index sets" begin
        m = Model()
        S = [:]
        @test_throws ErrorException @variable(m, x[S])
    end

    @testset "getindex constraint" begin
        m = Model()
        @variable(m, x)
        _c1 = @constraint(m, c1, x == 1)
        _c2 = @constraint(m, c2[i=1:3], x == i)
        _c3 = @NLconstraint(m, c3, x^3 == 1)
        _c4 = @NLconstraint(m, c4[i=1:3], x^i == 1)
        @test _c1 == m[:c1]
        @test _c2 == m[:c2]
        @test _c3 == m[:c3]
        @test _c4 == m[:c4]
    end

    @testset "Anonymous singleton variables" begin
        m = Model()
        x = @variable(m)
        y = @variable(m, lowerbound=0, upperbound=1, category = :Int)
        @test x == Variable(m, 1)
        @test y == Variable(m, 2)
        @test getcategory(y) == :Int
    end

    @testset "Invalid variable names" begin
        m = Model()
        if VERSION < v"0.7-"
            @test macroexpand(:(@variable(m, Bin))).head == :error
            @test macroexpand(:(@variable(m, Int))).head == :error
            @test macroexpand(:(@variable(m, Cont))).head == :error
            @test macroexpand(:(@variable(m, SemiCont))).head == :error
            @test macroexpand(:(@variable(m, SemiInt))).head == :error
        else
            @test_throws LoadError @macroexpand @variable(m, Bin)
            @test_throws LoadError @macroexpand @variable(m, Int)
            @test_throws LoadError @macroexpand @variable(m, Cont)
            @test_throws LoadError @macroexpand @variable(m, SemiCont)
            @test_throws LoadError @macroexpand @variable(m, SemiInt)
        end
    end

    @testset "Invalid lb/ub in ranged row" begin
        m = Model()
        @variable m x
        @variable m y
        @test_throws ErrorException @constraint m 0 <= x <= y
        @test_throws ErrorException @constraint m 2*x <= y <= 1
        @test_throws ErrorException @constraint m x+y <= x <= x*y
    end

    @testset "Adding vector constraints" begin
        m = Model()
        @variable(m, x[1:2, 1:2])
        @variable(m, y[1:2])
        u = [2, 3]
        v = [4, 5]
        expr_base = zero(JuMP.AffExpr)
        @constraint(m, x*u .<= x*v)
        @test string(m.linconstr[1]) == "-2 x[1,1] - 2 x[1,2] $leq 0"
        @test string(m.linconstr[2]) == "-2 x[2,1] - 2 x[2,2] $leq 0"
        @constraint(m, x*u + y .<= v)
        @test string(m.linconstr[3]) == "2 x[1,1] + 3 x[1,2] + y[1] $leq 4"
        @test string(m.linconstr[4]) == "2 x[2,1] + 3 x[2,2] + y[2] $leq 5"
        expr = JuMP.addtoexpr(expr_base, x, u)
        expr2 = JuMP.addtoexpr(expr, y, 1.0)
        @test expr2 == [2.0x[1,1] + 3.0x[1,2] + y[1];
                        2.0x[2,1] + 3.0x[2,2] + y[2]]
    end

    @testset "Filling @expression with a number" begin
        m = Model()
        @expression(m, y[1:3], 0.0)
        sites = [:A,:B,:C]
        @expression(m, z[i=sites], 0.0)
        @test y[1] == y[2] == y[3] == AffExpr(0.0)
        @test z[:A] == z[:B] == z[:C] == AffExpr(0.0)
    end

    @testset "@LinearConstraints" begin
        m = Model()
        @variable(m, x[1:3])
        lc = @LinearConstraints(begin
           x[1] + x[2] + x[3] ≥ 1
           x[1] + x[2] + x[3] ≤ 2
        end)
       @test lc[1].terms == lc[2].terms == x[1] + x[2] + x[3]
       @test lc[1].lb == 1
       @test lc[1].ub == Inf
       @test lc[2].lb == -Inf
       @test lc[2].ub == 2
    end

    @testset "Set constraint syntax" begin
        m = Model()
        @variable(m, x)

        JuMP.constructconstraint!(aff, ::__Cone__) = (m.ext[:ConeTest] = 1; LinearConstraint(aff, 0, 0))

        @constraint(m, 2x in __Cone__())
        @test m.ext[:ConeTest] == 1
    end

    @testset "Binary variable with invalid bounds" begin
        m = Model()
        @test_throws ArgumentError @variable(m, 2 <= x <= 3, Bin)
        @test_throws ArgumentError @variable(m, 2 <= x <= 3, category=:Bin)
        cat = :Bin
        @test_throws ArgumentError @variable(m, 2 <= x <= 3, category=cat)
        @test m.numCols == 0
        l = 0
        u = 1
        @variable(m, l <= x <= u, category=:Bin)
        @test getcategory(x) == :Bin
        @test getlowerbound(x) == 0
        @test getupperbound(x) == 1
        @variable(m, -1 <= y <= Inf, category=:Bin)
        @test getcategory(y) == :Bin
        @test getlowerbound(y) == 0
        @test getupperbound(y) == 1
        @variable(m, -Inf <= z <= Inf, category=:Bin)
        @test getcategory(z) == :Bin
        @test getlowerbound(z) == 0
        @test getupperbound(z) == 1
    end

    @testset "Extension of @variable with constructvariable! #1029" begin
        JuMP.variabletype(m::Model, ::Type{MyVariable}) = MyVariable
        function JuMP.constructvariable!(m::Model, ::Type{MyVariable}, _error::Function, lowerbound::Number, upperbound::Number, category::Symbol, basename::AbstractString, start::Number; test_kw::Int = 0)
            MyVariable(lowerbound, upperbound, category, basename, start, test_kw)
        end
        m = Model()
        @variable(m, 1 <= x <= 2, MyVariable, category = :Bin, test_kw = 1, start = 3)
        @test isa(x, MyVariable)
        @test x.lowerbound == 1
        @test x.upperbound == 2
        @test x.category == :Bin
        @test x.basename == "x"
        @test x.start == 3
        @test x.test_kw == 1
        @variable(m, y[1:3] >= 0, MyVariable, test_kw = 2)
        @test isa(y, Vector{MyVariable})
        for i in 1:3
            @test y[i].lowerbound == 0
            @test y[i].upperbound == Inf
            @test y[i].category == :Default
            @test isempty(y[i].basename)
            @test isnan(y[i].start)
            @test y[i].test_kw == 2
        end
    end

    @testset "constructconstraint! on variable" begin
        m = Model()
        @variable(m, x)
        @test string(JuMP.constructconstraint!(x, :(>=))) == "x $geq 0"
        @test string(JuMP.constructconstraint!(x, :(<=))) == "x $leq 0"
        @test string(JuMP.constructconstraint!(x, :(==))) == "x $eq 0"
    end

    @testset "Nested tuple destructuring" begin
    m = Model()
    @variable(m, x)
    d = Dict((1,2) => 3)
    @constraint(m, x == sum(i+j+k for ((i,j),k) in d))
    @test m.linconstr[1].lb == m.linconstr[1].ub == 6
    @test m.linconstr[1].terms == 1x
    end
end
