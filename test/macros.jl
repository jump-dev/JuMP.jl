# test/macros.jl
# Testing macros work correctly
#############################################################################
using JuMP, FactCheck
using Base.Test

# To ensure the tests work on Windows and Linux/OSX, we need
# to use the correct comparison operators
const leq = JuMP.repl[:leq]
const geq = JuMP.repl[:geq]
const  eq = JuMP.repl[:eq]
const Vert = JuMP.repl[:Vert]
const sub2 = JuMP.repl[:sub2]

facts("[macros] Check Julia expression parsing") do
    sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j in 1:M; i != j})
    @fact sumexpr.head --> :curly
    @fact length(sumexpr.args) --> 5
    @fact sumexpr.args[1] --> :sum
    @fact sumexpr.args[2].head --> :parameters
    @fact sumexpr.args[3] --> :(x[i,j] * y[i,j])
    @fact sumexpr.args[4].head --> :(=)
    if VERSION < v"0.5.0-dev+3231"
        @fact sumexpr.args[5].head --> :in
    else
        @fact sumexpr.args[5].head --> :call
        @fact sumexpr.args[5].args[1] --> :in
    end

    sumexpr = :(sum{x[i,j] * y[i,j], i in 1:N, j = 1:M})
    @fact sumexpr.head --> :curly
    @fact length(sumexpr.args) --> 4
    @fact sumexpr.args[1] --> :sum
    @fact sumexpr.args[2] --> :(x[i,j] * y[i,j])
    if VERSION < v"0.5.0-dev+3231"
        @fact sumexpr.args[3].head --> :in
    else
        @fact sumexpr.args[3].head --> :call
        @fact sumexpr.args[3].args[1] --> :in
    end

    @fact sumexpr.args[4].head --> :(=)
end

# generator syntax only parses on 0.5
if VERSION >= v"0.5-dev+5475"
    eval("""
facts("[macros] Check Julia expression parsing (0.5)") do
    sumexpr = :(sum(x[i,j] * y[i,j] for i = 1:N, j in 1:M if i != j))
    @fact sumexpr.head --> :call
    @fact sumexpr.args[1] --> :sum
    @fact sumexpr.args[2].head --> :generator
    @fact sumexpr.args[2].args[1] --> :(x[i,j] * y[i,j])
    @fact sumexpr.args[2].args[2].head --> :filter
    @fact sumexpr.args[2].args[2].args[1] --> :(i != j)
    @fact sumexpr.args[2].args[2].args[2] --> :(i = 1:N)
    @fact sumexpr.args[2].args[2].args[3] --> :(j = 1:M)

    sumexpr = :(sum(x[i,j] * y[i,j] for i = 1:N, j in 1:M))
    @fact sumexpr.head --> :call
    @fact sumexpr.args[1] --> :sum
    @fact sumexpr.args[2].head --> :generator
    @fact sumexpr.args[2].args[1] --> :(x[i,j] * y[i,j])
    @fact sumexpr.args[2].args[2] --> :(i = 1:N)
    @fact sumexpr.args[2].args[3] --> :(j = 1:M)
end
"""); end

facts("[macros] Check @constraint basics") do
    m = Model()
    @variable(m, w)
    @variable(m, x)
    @variable(m, y)
    @variable(m, z)
    t = 10

    @constraint(m, 3x - y == 3.3(w + 2z) + 5)
    @fact string(m.linconstr[end]) --> "3 x - y - 3.3 w - 6.6 z $eq 5"
    @constraint(m, 3x - y == (w + 2z)*3.3 + 5)
    @fact string(m.linconstr[end]) --> "3 x - y - 3.3 w - 6.6 z $eq 5"
    @constraint(m, (x+y)/2 == 1)
    @fact string(m.linconstr[end]) --> "0.5 x + 0.5 y $eq 1"
    @constraint(m, -1 <= x-y <= t)
    @fact string(m.linconstr[end]) --> "-1 $leq x - y $leq 10"
    @constraint(m, -1 <= x+1 <= 1)
    @fact string(m.linconstr[end]) --> "-2 $leq x $leq 0"
    @constraint(m, -1 <= x <= 1)
    @fact string(m.linconstr[end]) --> "-1 $leq x $leq 1"
    @constraint(m, -1 <= x <= sum{0.5, i = 1:2})
    @fact string(m.linconstr[end]) --> "-1 $leq x $leq 1"
    @fact_throws @constraint(m, x <= t <= y)
    @fact macroexpand(:(@constraint(m, 1 >= x >= 0))).head --> :error
    @fact macroexpand(:(@constraint(1 <= x <= 2, foo=:bar))).head --> :error

    @expression(m, aff, 3x - y - 3.3(w + 2z) + 5)
    @fact string(aff) --> "3 x - y - 3.3 w - 6.6 z + 5"

    @constraint(m, 3 + 5*7 <= 0)
    @fact string(m.linconstr[end]) --> "0 $leq -38"

    @expression(m, qaff, (w+3)*(2x+1)+10)
    @fact string(qaff) --> "2 w*x + 6 x + w + 13"
end

facts("[macros] Checking @variable with reverse direction bounds") do
    m = Model()
    @variable(m, 3.2 >= x >= 1)
    @fact m.colLower --> [1.0]
    @fact m.colUpper --> [3.2]
end

facts("[macros] sum{}") do
    m = Model()
    @variable(m, x[1:3,1:3])
    @variable(m, y)
    C = [1 2 3; 4 5 6; 7 8 9]
    @constraint(m, sum{ C[i,j]*x[i,j], i in 1:2, j = 2:3 } <= 1)
    @fact string(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 5 x[2,2] + 6 x[2,3] $leq 1"
    @constraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j in 1:3; i != j} == y)
    @fact string(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

    @constraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:i} == 0);
    @fact string(m.linconstr[end]) --> "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] $eq 0"

    @constraint(m, sum{ 0*x[i,1], i=1:3} == 0)
    @fact string(m.linconstr[end]) --> "0 $eq 0"

    @constraint(m, sum{ 0*x[i,1] + y, i=1:3} == 0)
    @fact string(m.linconstr[end]) --> "3 y $eq 0"

end

if VERSION >= v"0.5-dev+5475"
eval("""
facts("[macros] sum(generator)") do
    m = Model()
    @variable(m, x[1:3,1:3])
    @variable(m, y)
    C = [1 2 3; 4 5 6; 7 8 9]
    @constraint(m, sum( C[i,j]*x[i,j] for i in 1:2, j = 2:3 ) <= 1)
    @fact string(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 5 x[2,2] + 6 x[2,3] $leq 1"
    @constraint(m, sum( C[i,j]*x[i,j] for i = 1:3, j in 1:3 if i != j) == y)
    @fact string(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

    @constraint(m, sum( C[i,j]*x[i,j] for i = 1:3, j = 1:i) == 0);
    @fact string(m.linconstr[end]) --> "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] $eq 0"

    @constraint(m, sum( 0*x[i,1] for i=1:3) == 0)
    @fact string(m.linconstr[end]) --> "0 $eq 0"

    @constraint(m, sum( 0*x[i,1] + y for i=1:3) == 0)
    @fact string(m.linconstr[end]) --> "3 y $eq 0"

    @fact isexpr(macroexpand(:(@constraint(m, sum( 0*x[i,1] + y for i=1:3 for j in 1:3) == 0))),:error) --> true

end"""); end

facts("[macros] Problem modification") do
    m = Model()
    @variable(m, x[1:3,1:3])
    C = [1 2 3; 4 5 6; 7 8 9]

    @constraint(m, sum{ x[i,j]*(C[i,j]-1), i in 1:3, j = 1:3; i != j} == 0)
    @fact string(m.linconstr[end]) --> "x[1,2] + 2 x[1,3] + 3 x[2,1] + 5 x[2,3] + 6 x[3,1] + 7 x[3,2] $eq 0"

    con = @constraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:3; i != j} == 0)
    @fact string(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] $eq 0"

    @variable(m, y, objective = 0, inconstraints = [con], coefficients = [-1.0])
    @fact string(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

    JuMP.setRHS(con, 3)
    @fact string(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 3"
end

facts("[macros] Using pre-built affine is OK in macro") do
    m = Model()
    @variable(m, x)
    @variable(m, y)
    temp = x + 2y + 1
    @constraint(m, 3*temp - x - 2 >= 0)
    @fact string(m.linconstr[end]) --> "6 y + 2 x $geq -1"
    # More complex expression
    a = 1.0*x
    @constraint(m, (2+2)*((3+4)*(1+a)) == 0)
    @fact string(m.linconstr[end]) --> "28 x $eq -28"
    @fact string(a) --> "x"

    @fact string(@LinearConstraint(1 + 0*temp == 0)) --> "0 $eq -1"
end

facts("[macros] Test ranges in @variable") do
    m = Model()
    @variable(m, x[1:5])
    @variable(m, y[3:2:9])
    @variable(m, z[4:3:8])
    @variable(m, w[6:5])

    @fact length(x) --> 5
    @fact length(y) --> 4
    @fact length(z) --> 2
    @fact length(w) --> 0

    @fact x[end].col --> x[5].col
    @fact y[3].m --> y[5].m
    @fact y[3].m --> y[7].m
    @fact y[3].m --> y[9].m
    @fact_throws z[8].col  # KeyError
    @fact_throws w[end]  # BoundsError
end

facts("[macros] Unicode comparisons") do
    m = Model()
    @variable(m, 0 ≤ x ≤ 1)
    @variable(m, y ≥ 2)
    @variable(m, z ≤ 3)
    @fact m.colUpper --> [1.0, Inf,  3.0]
    @fact m.colLower --> [0.0, 2.0, -Inf]
    @constraint(m, 0 ≤ x + y ≤ 1)
    @constraint(m, x + z ≤ 2)
    @constraint(m, y + z ≥ 3)
    @constraint(m, y*z ≤ 1)
    @fact m.linconstr[1].lb --> 0.0
    @fact m.linconstr[1].ub --> 1.0
    @fact m.linconstr[2].lb --> -Inf
    @fact m.linconstr[2].ub --> 2.0
    @fact m.linconstr[3].lb --> 3.0
    @fact m.linconstr[3].ub --> Inf
    @fact m.quadconstr[1].sense --> :(<=)
end

facts("[macros] Three argument @constraint") do
    m = Model()
    @variable(m, x[1:5])
    @variable(m, y[2:2:6])

    @constraint(m, c, x[4] - y[4] == 1)
    @fact string(m.linconstr[c.idx]) --> "x[4] - y[4] $eq 1"

    @constraint(m, d[i in 1:5,j=6:-2:2], x[i] - y[j] == 2)
    @fact string(m.linconstr[d[4,4].idx]) --> "x[4] - y[4] $eq 2"

    @constraint(m, q[i=1:5], x[i]^2 == 1)
    @fact string(m.quadconstr[q[5].idx]) --> "x[5]² - 1 $eq 0"
end

facts("[macros] @constraints") do
    m = Model()
    @variable(m, x)
    @variable(m, y[1:3])

    @constraints(m, begin
        x + y[1] == 1
        ref[i=1:3], y[1] + y[i] >= i
    end)

    @fact string(m.linconstr[1]) --> "x + y[1] $eq 1"
    @fact string(m.linconstr[2]) --> "2 y[1] $geq 1"
    @fact string(m.linconstr[3]) --> "y[1] + y[2] $geq 2"
    @fact string(m.linconstr[4]) --> "y[1] + y[3] $geq 3"
end

facts("[macros] @NLconstraints") do
    m = Model()
    @variable(m, 0 <= x <= 1)
    @variable(m, y[1:3])
    @objective(m, Max, x)

    @NLconstraints(m, begin
        ref[i=1:3], y[i] == 0
        x + y[1] * y[2] * y[3] <= 0.5
    end)

    @fact length(m.nlpdata.nlconstr) --> 4
    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:ExprGraph])

    @fact MathProgBase.constr_expr(d,1) --> :(x[2] - 0.0 == 0.0)
    @fact MathProgBase.constr_expr(d,2) --> :(x[3] - 0.0 == 0.0)
    @fact MathProgBase.constr_expr(d,3) --> :(x[4] - 0.0 == 0.0)
    @fact MathProgBase.constr_expr(d,4) --> :((x[1] + x[2] * x[3] * x[4]) - 0.5 <= 0.0)

end

facts("[macros] Vectors in nonlinear expressions") do
    m = Model()
    @variable(m, x[1:3])
    @fact_throws ErrorException @NLobjective(m, Min, x)
    @fact_throws ErrorException @NLobjective(m, Min, [1,2,3])
end

facts("[macros] @objective with quadratic") do
    m = Model()
    @variable(m, x[1:5])
    @objective(m, Max, sum{i*x[i]*x[j], i=1:5, j=5:-1:1; isodd(i) && iseven(j)} + 2x[5])

    @fact string(m.obj) --> "x[1]*x[2] + 3 x[2]*x[3] + x[1]*x[4] + 3 x[3]*x[4] + 5 x[2]*x[5] + 5 x[4]*x[5] + 2 x[5]"
end

facts("[macros] @constraint with quadratic") do
    m = Model()
    @variable(m, x[1:5])

    @constraint(m, x[3]*x[1] + sum{x[i]*x[5-i+1], i=1:5; 2 <= i <= 4} + 4x[5] == 1)
    @fact string(m.quadconstr[end]) --> "x[1]*x[3] + x[3]² + 2 x[2]*x[4] + 4 x[5] - 1 $eq 0"

    @constraint(m, sum{sum{(x[i] - 2)*x[j],j=4:5},i=2:3} >= -3*x[2]*2*x[4])
    @fact string(m.quadconstr[end]) --> "7 x[2]*x[4] + x[3]*x[4] + x[2]*x[5] + x[3]*x[5] - 4 x[4] - 4 x[5] $geq 0"

    foo(x) = x
    @constraint(m, x[1] ≤ foo(x[1])^2)
    @fact string(m.quadconstr[end]) --> "-x[1]² + x[1] $leq 0"

    @constraint(m, sum{x[i],i=1:2}*sum{x[i],i=2:3} >= 0)
    @fact string(m.quadconstr[end]) --> "x[1]*x[2] + x[2]² + x[1]*x[3] + x[2]*x[3] $geq 0"
    @constraint(m, x[1]^2 + x[2]*x[3] >= 0)
    @fact string(m.quadconstr[end]) --> "x[1]² + x[2]*x[3] $geq 0"
    @constraint(m, x[1]^2 + (x[2]+3)*(x[3]-1) >= 0)
    @fact string(m.quadconstr[end]) --> "x[1]² + x[2]*x[3] + 3 x[3] - x[2] - 3 $geq 0"
    @constraint(m, sum{x[i],i=1:2}^2 >= 0)
    @fact string(m.quadconstr[end]) --> "x[1]² + 2 x[1]*x[2] + x[2]² $geq 0"

    myquadexpr = x[1]*x[2]
    @constraint(m, sum{i*myquadexpr + x[i], i=1:3} + sum{x[i] + myquadexpr*i, i=1:3} == 0)
    @fact string(m.quadconstr[end]) --> "12 x[1]*x[2] + 2 x[1] + 2 x[2] + 2 x[3] $eq 0"

    @constraint(m, (x[1] + x[2])*sum{ 0*x[i] + x[3], i=1:3} == 0)
    @fact string(m.quadconstr[end]) --> "3 x[1]*x[3] + 3 x[2]*x[3] $eq 0"

    @fact string(@QuadConstraint(1 + 0*myquadexpr == 0)) --> "1 $eq 0"

    @variable(m, y)
    @fact string(@QuadConstraint(1 + (2y)*y   == 0)) --> "2 y² + 1 $eq 0"
    @fact string(@QuadConstraint(1 +   y *y*2 == 0)) --> "2 y² + 1 $eq 0"
    z = 2y
    @fact string(@QuadConstraint(y*y + y*z == 0)) --> "3 y² $eq 0"
end

facts("[macros] Triangular indexing, iteration") do
    n = 10
    trimod = Model()
    @variable(trimod, x[i=1:n,j=i:n])
    @variable(trimod, y[i=3:2:7,j=-i])
    @fact MathProgBase.numvar(trimod) --> n*(n+1)/2 + 3
    S = Any[(i,i+2) for i in 1:5]
    @variable(trimod, z[(i,j)=S,k=i:j])
    @fact length(z.tupledict) --> 15
    @constraint(trimod, cref[i=1:n,j=i:n], x[i,j] + y[5,-5] == 1)
    @fact MathProgBase.numconstr(trimod) --> n*(n+1)/2

    cntr = zeros(Bool, n, n)
    for ((i,j),var) in zip(keys(x),values(x))
        @fact x[i,j] --> exactly(var)
        cntr[i,j] = true
    end
    for i in 1:n, j in 1:n
        @fact cntr[i,j] --> (j >= i)
    end
end

facts("[macros] Multidimensional indexing") do
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
    @fact match --> true
end

facts("[macros] @expression") do
    model = Model()
    @variable(model, x[1:3,1:3])
    @expression(model, expr, sum{i*x[i,j] + j, i=1:3,j in 1:3})
    @fact string(expr) --> "x[1,1] + x[1,2] + x[1,3] + 2 x[2,1] + 2 x[2,2] + 2 x[2,3] + 3 x[3,1] + 3 x[3,2] + 3 x[3,3] + 18"

    @fact_throws @expression(model, blah[i=1:3], x[i,1]^2)

    @expression(model, y[i=1:2], sum{x[i,1]; i == 1})
    @fact string(y[1]) --> "x[1,1]"
    @fact string(y[2]) --> "0"

    # deprecated versions
    @expression(expr2, sum{i*x[i,j] + j, i=1:3,j in 1:3})
    @fact string(expr2) --> "x[1,1] + x[1,2] + x[1,3] + 2 x[2,1] + 2 x[2,2] + 2 x[2,3] + 3 x[3,1] + 3 x[3,2] + 3 x[3,3] + 18"
    expr2 = @expression(sum{i*x[i,j] + j, i=1:3,j in 1:3})
    @fact string(expr2) --> "x[1,1] + x[1,2] + x[1,3] + 2 x[2,1] + 2 x[2,2] + 2 x[2,3] + 3 x[3,1] + 3 x[3,2] + 3 x[3,3] + 18"
end

facts("[macros] Conditions in constraint indexing") do
    model = Model()
    @variable(model, x[1:10])
    @constraint(model, c1[i=1:9;isodd(i)], x[i] + x[i+1] <= 1)
    @NLconstraint(model, c2[i=["red","blue","green"], k=9:-2:2; (i == "red" && isodd(k)) || (k >=4 && (i == "blue" || i == "green"))], x[k]^3 <= 1)
    @fact length(model.linconstr) --> 5
    @fact string(model.linconstr[1]) --> "x[1] + x[2] $leq 1"
    @fact string(model.linconstr[2]) --> "x[3] + x[4] $leq 1"
    @fact string(model.linconstr[3]) --> "x[5] + x[6] $leq 1"
    @fact string(model.linconstr[4]) --> "x[7] + x[8] $leq 1"
    @fact string(model.linconstr[5]) --> "x[9] + x[10] $leq 1"
    @fact length(model.nlpdata.nlconstr) --> 10
end

facts("[macros] Test changes in condition parsing") do
    ex = :(x[12;3])
    @fact ex.head --> :typed_vcat
    @fact ex.args --> [:x, 12, 3]

    ex = :(x[i=1:3,j=S;isodd(i) && i+j>=2])
    @fact ex.head --> :typed_vcat
    @fact ex.args --> [:x,
                      Expr(:parameters, Expr(:&&, :(isodd(i)), :(i+j>=2))),
                      Expr(:(=), :i, :(1:3)),
                      Expr(:(=), :j, :S)]
end

facts("[macros] Norm parsing") do
    model = Model()
    @variable(model, x[1:2,1:2])
    @constraint(model, -2norm2{x[i,j], i in 1:2, j=1:2} + x[1,2] >= -1)
    @constraint(model, -2norm2{x[i,j], i=1:2, j in 1:2; iseven(i+j)} + x[1,2] >= -1)
    @constraint(model, 1 >= 2*norm2{x[i,1], i in 1:2})
    @fact string(model.socconstr[1]) --> "2.0 $Vert[x[1,1],x[1,2],x[2,1],x[2,2]]$Vert$sub2 $leq x[1,2] + 1"
    @fact string(model.socconstr[2]) --> "2.0 $Vert[x[1,1],x[2,2]]$Vert$sub2 $leq x[1,2] + 1"
    @fact string(model.socconstr[3]) --> "2.0 $Vert[x[1,1],x[2,1]]$Vert$sub2 $leq 1"
    @fact_throws @constraint(model, (x[1,1]+1)*norm2{x[i,j], i=1:2, j=1:2} + x[1,2] >= -1)
    @fact_throws @constraint(model, norm2{x[i,j], i=1:2, j=1:2} + x[1,2] >= -1)
end

if VERSION >= v"0.5-dev+5475"
eval("""
facts("[macros] Generator norm parsing") do
    model = Model()
    @variable(model, x[1:2,1:2])
    @constraint(model, -2norm2(x[i,j] for i in 1:2, j=1:2) + x[1,2] >= -1)
    @constraint(model, -2norm2(x[i,j] for i=1:2, j in 1:2 if iseven(i+j)) + x[1,2] >= -1)
    @constraint(model, 1 >= 2*norm2(x[i,1] for i in 1:2))
    @fact string(model.socconstr[1]) --> "2.0 $Vert[x[1,1],x[1,2],x[2,1],x[2,2]]$Vert$sub2 $leq x[1,2] + 1"
    @fact string(model.socconstr[2]) --> "2.0 $Vert[x[1,1],x[2,2]]$Vert$sub2 $leq x[1,2] + 1"
    @fact string(model.socconstr[3]) --> "2.0 $Vert[x[1,1],x[2,1]]$Vert$sub2 $leq 1"
    @fact_throws @constraint(model, (x[1,1]+1)*norm2(x[i,j] for i=1:2, j=1:2) + x[1,2] >= -1)
    @fact_throws @constraint(model, norm2(x[i,j] for i=1:2, j=1:2) + x[1,2] >= -1)
end""");end

facts("[macros] Extraneous terms in QuadExpr (#535)") do
    model = Model()
    @variable(model, x)
    @variable(model, y)
    @constraint(model, x*x <= y*y)
    @fact string(model.quadconstr[1]) --> "x² - y² $leq 0"
end

facts("[macros] Special-case binary multiplication in addtoexpr_reorder (#537)") do
    dual = Model()
    @variable(dual, α[1:3] <= 0)
    @variable(dual, γ >= 0)
    @constraint(dual, ones(3)*γ .<= α)
    @fact string(dual.linconstr[1]) --> "γ - α[1] $leq 0"
    @fact string(dual.linconstr[2]) --> "γ - α[2] $leq 0"
    @fact string(dual.linconstr[3]) --> "γ - α[3] $leq 0"
end

facts("[macros] Indices in macros don't leak out of scope (#582)") do
    m = Model()
    cnt = 4
    for i in 5:8
        @variable(m, x[i=1:3,j=1:3] ≤ i)
        cnt += 1
        @fact i --> cnt
    end
    cnt = 4
    for i in 5:8
        @variable(m, y[i=2:4,j=1:3] ≤ i)
        cnt += 1
        @fact i --> cnt
    end
    cnt = 4
    for i in 5:8
        @variable(m, z[i=[1:3;],j=1:3] ≤ i)
        cnt += 1
        @fact i --> cnt
    end
    @fact m.colUpper --> vcat(repeat([1.0,2.0,3.0], inner=[3], outer=[4]),
                              repeat([2.0,3.0,4.0], inner=[3], outer=[4]),
                              repeat([1.0,2.0,3.0], inner=[3], outer=[4]))
    @variable(m, x[i=1:3] ≤ i)
    cnt = 4
    for i in 5:8
        @constraint(m, sum{x[i], i=1, j=1, k=1} == 1)
        cnt += 1
        @fact i --> cnt
    end
    cnt = 4
    for i in 5:8
        @constraint(m, norm2{x[i], i=1, j=1, k=1} <= 1)
        cnt += 1
        @fact i --> cnt
    end
    cnt = 4
    for i in 5:8
        @constraint(m, c[i=1:3,j=1:3], x[i] == 1)
        cnt += 1
        @fact i --> cnt
    end
end

facts("[macros] Issue #621") do
    m = Model()
    @variable(m, x)

    q = x^2
    a = x+1
    con = @QuadConstraint(a+q*3 <= 0)
    @fact string(con) --> "3 x² + x + 1 $leq 0"
end

facts("[macros] @variables and @constraints") do
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
    @fact m.colUpper --> [1.0, 2.0, Inf,  3.0, 1.0]
    @fact m.colLower --> [0.0, 0.0, 2.0, -Inf, 0.0]
    @fact m.colCat --> [:Cont, :Cont, :Int, :Cont, :Bin]
    @fact getvalue(y) --> 0.7
    @fact getvalue(z) --> 10
    @fact getvalue(q) --> 0.5
    @fact m.linconstr[1].lb --> 0.0
    @fact m.linconstr[1].ub --> 1.0
    @fact m.linconstr[2].lb --> -Inf
    @fact m.linconstr[2].ub --> 2.0
    @fact m.linconstr[3].lb --> 3.0
    @fact m.linconstr[3].ub --> Inf
    @fact m.quadconstr[1].sense --> :(<=)
end

facts("[macros] @expressions and @NLexpressions") do
    m = Model()
    @variable(m, x)

    @expressions(m, begin
        myex[i=1:2], x + i
        myex2,       x + 3
    end)
    setvalue(x, 1)
    @fact getvalue(myex[1]) --> 2
    @fact getvalue(myex[2]) --> 3
    @fact getvalue(myex2)   --> 4
    @fact_throws getvalue(myex[3])

    @NLexpressions(m, begin
        nlex,         x^2
        nlex2[i=1:2], x^i + sqrt(x)
    end)
    setvalue(x, 2)
    @fact getvalue(nlex) --> 4
    @fact getvalue(nlex2[1]) --> 2^1 + sqrt(2)
    @fact getvalue(nlex2[2]) --> 2^2 + sqrt(2)
    @fact_throws getvalue(nlex2[3])
end

facts("[macros] No bare symbols in constraint macros") do
    m = Model()
    @variable(m, x)
    @fact macroexpand(:(@constraint(m, x))).head --> :error
    @fact macroexpand(:(@constraint(m, :foo))).head --> :error
    @fact macroexpand(:(@SDconstraint(m, x))).head --> :error
    @fact macroexpand(:(@SDconstraint(m, :foo))).head --> :error
    @fact macroexpand(:(@NLconstraint(m, x))).head --> :error
    @fact macroexpand(:(@NLconstraint(m, :foo))).head --> :error
end

facts("[macros] LB/UB kwargs") do
    m = Model()
    @variable(m, a, lowerbound=0)
    @variable(m, b, lowerbound=0, upperbound=1)
    @variable(m, c, upperbound=1, lowerbound=0)
    @variable(m, d, upperbound=1)
    @variable(m, e >= 0, upperbound=1)
    @variable(m, f <= 1, lowerbound=0)
    @fact m.colLower --> [0,0,0,-Inf,0,0]
    @fact m.colUpper --> [Inf,1,1,1,1,1]

    @fact macroexpand(:(@variable(m, g >= 0, lowerbound=1))).head --> :error
    @fact macroexpand(:(@variable(m, h <= 1, upperbound=1))).head --> :error
    @fact macroexpand(:(@variable(m, 0 <= i <= 1, lowerbound=1))).head --> :error
    @fact macroexpand(:(@variable(m, 0 <= j <= 1, upperbound=1))).head --> :error
end

facts("[macros] Anonymous versions of macros") do
    m = Model()
    x = @variable(m, [1:3], lowerbound=1.0)
    y = @variable(m, [2:3,1:3], upperbound=1.0, lowerbound=0.0, Bin)
    z = @variable(m, [[:red,:blue]], Int)
    w = @variable(m, [i=1:4,j=1:4;isodd(i+j)], SemiCont)
    # v = @variable(m, [i=1:3,j=1:3], Symmetric, lowerbound = eye(3)[i,j])
    u = @variable(m, [1:4,1:4], SDP)
    @fact macroexpand(:(@variable(m, [1:3] <= 1))).head --> :error

    @fact getlowerbound(x[1]) --> 1.0
    @fact getupperbound(x[1]) --> Inf
    @fact getcategory(x[1]) --> :Cont
    @fact getlowerbound(y[2,1]) --> 0.0
    @fact getupperbound(y[2,1]) --> 1.0
    @fact getcategory(y[2,1]) --> :Bin
    @fact getlowerbound(z[:red]) --> -Inf
    @fact getupperbound(z[:red]) --> Inf
    @fact getcategory(z[:red]) --> :Int
    @fact getlowerbound(w[1,2]) --> -Inf
    @fact getupperbound(w[1,2]) --> Inf
    @fact getcategory(w[1,2]) --> :SemiCont
    # @fact getlowerbound(v[1,2]) --> 0.0
    # @fact getupperbound(v[1,2]) --> Inf
    # @fact getcategory(v[1,2]) --> :Cont
    @fact getlowerbound(u[1,2]) --> -Inf
    @fact getupperbound(u[1,2]) --> Inf
    @fact getcategory(u[1,2]) --> :Cont

    c = @constraint(m, [i=1:3], x[i] <= z[:red])
    d = @NLconstraint(m, [i=1:3], x[i]^3 == 1)
    e = @NLexpression(m, [i=2:3], y[2,i]^3)
    f = @expression(m, [i=[:red,:blue]], u[1,2] + 2z[i])

    # not sure how else to test this
    @fact c[1] --> c[1]
    @fact d[1] --> d[1]
    @fact e[2] --> e[2]
    @fact f[:red] --> f[:red]
end

facts("[macros] Colons in index sets") do
    m = Model()
    S = [:]
    @fact_throws ErrorException @variable(m, x[S])
end

facts("[macros] getconstraint") do
    m = Model()
    @variable(m, x)
    _c1 = @constraint(m, c1, x == 1)
    _c2 = @constraint(m, c2[i=1:3], x == i)
    _c3 = @NLconstraint(m, c3, x^3 == 1)
    _c4 = @NLconstraint(m, c4[i=1:3], x^i == 1)
    @fact _c1 --> getconstraint(m, :c1)
    @fact _c2 --> getconstraint(m, :c2)
    @fact _c3 --> getconstraint(m, :c3)
    @fact _c4 --> getconstraint(m, :c4)
end

facts("[macros] Anonymous singleton variables") do
    m = Model()
    x = @variable(m)
    y = @variable(m, lowerbound=0, upperbound=1)
    @fact x --> Variable(m, 1)
    @fact y --> Variable(m, 2)
end

facts("[macros] Invalid variable names") do
    m = Model()
    @fact macroexpand(:(@variable(m, Bin))).head --> :error
    @fact macroexpand(:(@variable(m, Int))).head --> :error
    @fact macroexpand(:(@variable(m, Cont))).head --> :error
    @fact macroexpand(:(@variable(m, SemiCont))).head --> :error
    @fact macroexpand(:(@variable(m, SemiInt))).head --> :error
end
