# test/macros.jl
# Testing macros work correctly
#############################################################################
using JuMP, FactCheck
using Base.Test

# To ensure the tests work on Windows and Linux/OSX, we need
# to use the correct comparison operators
const leq = JuMP.repl_leq
const geq = JuMP.repl_geq
const  eq = JuMP.repl_eq

facts("[macros] Check Julia expression parsing") do
    sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j = 1:M; i != j})
    @fact sumexpr.head --> :curly
    @fact length(sumexpr.args) --> 5
    @fact sumexpr.args[1] --> :sum
    @fact sumexpr.args[2].head --> :parameters
    @fact sumexpr.args[3] --> :(x[i,j] * y[i,j])
    @fact sumexpr.args[4].head --> :(=)
    @fact sumexpr.args[5].head --> :(=)

    sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j = 1:M})
    @fact sumexpr.head --> :curly
    @fact length(sumexpr.args) --> 4
    @fact sumexpr.args[1] --> :sum
    @fact sumexpr.args[2] --> :(x[i,j] * y[i,j])
    @fact sumexpr.args[3].head --> :(=)
    @fact sumexpr.args[4].head --> :(=)
end

facts("[macros] Check @addConstraint basics") do
    m = Model()
    @defVar(m, w)
    @defVar(m, x)
    @defVar(m, y)
    @defVar(m, z)
    t = 10

    @addConstraint(m, 3x - y == 3.3(w + 2z) + 5)
    @fact conToStr(m.linconstr[end]) --> "3 x - y - 3.3 w - 6.6 z $eq 5"
    if VERSION >= v"0.4.0-"
        @addConstraint(m, 3x - y == (w + 2z)*3.3 + 5)
        @fact conToStr(m.linconstr[end]) --> "3 x - y - 3.3 w - 6.6 z $eq 5"
    end
    @addConstraint(m, (x+y)/2 == 1)
    @fact conToStr(m.linconstr[end]) --> "0.5 x + 0.5 y $eq 1"
    @addConstraint(m, -1 <= x-y <= t)
    @fact conToStr(m.linconstr[end]) --> "-1 $leq x - y $leq 10"
    @addConstraint(m, -1 <= x+1 <= 1)
    @fact conToStr(m.linconstr[end]) --> "-2 $leq x $leq 0"
    @addConstraint(m, -1 <= x <= 1)
    @fact conToStr(m.linconstr[end]) --> "-1 $leq x $leq 1"
    if VERSION > v"0.4-"
        @addConstraint(m, -1 <= x <= sum{0.5, i = 1:2})
        @fact conToStr(m.linconstr[end]) --> "-1 $leq x $leq 1"
    end
    @fact_throws @addConstraint(m, x <= t <= y)

    @defExpr(aff, 3x - y - 3.3(w + 2z) + 5)
    @fact affToStr(aff) --> "3 x - y - 3.3 w - 6.6 z + 5"

    if VERSION >= v"0.4.0-"
        @defExpr(qaff, (w+3)*(2x+1)+10)
        @fact quadToStr(qaff) --> "2 w*x + 6 x + w + 13"
    end
end

facts("[macros] Checking @defVar with reverse direction bounds") do
    m = Model()
    @defVar(m, 3.2 >= x >= 1)
    @fact m.colLower --> [1.0]
    @fact m.colUpper --> [3.2]
end

facts("[macros] sum{}") do
    m = Model()
    @defVar(m, x[1:3,1:3])
    @defVar(m, y)
    C = [1 2 3; 4 5 6; 7 8 9]
    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:2, j = 2:3 } <= 1)
    @fact conToStr(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 5 x[2,2] + 6 x[2,3] $leq 1"
    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:3; i != j} == y)
    @fact conToStr(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:i} == 0);
    @fact conToStr(m.linconstr[end]) --> "x[1,1] + 4 x[2,1] + 5 x[2,2] + 7 x[3,1] + 8 x[3,2] + 9 x[3,3] $eq 0"
end

facts("[macros] Problem modification") do
    m = Model()
    @defVar(m, x[1:3,1:3])
    C = [1 2 3; 4 5 6; 7 8 9]

    if VERSION >= v"0.4.0-"
        @addConstraint(m, sum{ x[i,j]*(C[i,j]-1), i = 1:3, j = 1:3; i != j} == 0)
        @fact conToStr(m.linconstr[end]) --> "x[1,2] + 2 x[1,3] + 3 x[2,1] + 5 x[2,3] + 6 x[3,1] + 7 x[3,2] $eq 0"
    end

    con = @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:3; i != j} == 0)
    @fact conToStr(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] $eq 0"

    @defVar(m, y, objective = 0, inconstraints = [con], coefficients = [-1.0])
    @fact conToStr(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 0"

    chgConstrRHS(con, 3)
    @fact conToStr(m.linconstr[end]) --> "2 x[1,2] + 3 x[1,3] + 4 x[2,1] + 6 x[2,3] + 7 x[3,1] + 8 x[3,2] - y $eq 3"
end

facts("[macros] Using pre-built affine is OK in macro") do
    m = Model()
    @defVar(m, x)
    @defVar(m, y)
    temp = x + 2y + 1
    @addConstraint(m, 3*temp - x - 2 >= 0)
    @fact conToStr(m.linconstr[end]) --> "6 y + 2 x $geq -1"
    # More complex expression
    a = 1.0*x
    @addConstraint(m, (2+2)*((3+4)*(1+a)) == 0)
    @fact conToStr(m.linconstr[end]) --> "28 x $eq -28"
    @fact affToStr(a) --> "x"

end

facts("[macros] Test ranges in @defVar") do
    m = Model()
    @defVar(m, x[1:5])
    @defVar(m, y[3:2:9])
    @defVar(m, z[4:3:8])
    @defVar(m, w[6:5])

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
    @defVar(m, 0 ≤ x ≤ 1)
    @defVar(m, y ≥ 2)
    @defVar(m, z ≤ 3)
    @fact m.colUpper --> [1.0, Inf,  3.0]
    @fact m.colLower --> [0.0, 2.0, -Inf]
    @addConstraint(m, 0 ≤ x + y ≤ 1)
    @addConstraint(m, x + z ≤ 2)
    @addConstraint(m, y + z ≥ 3)
    @addConstraint(m, y*z ≤ 1)
    @fact m.linconstr[1].lb --> 0.0
    @fact m.linconstr[1].ub --> 1.0
    @fact m.linconstr[2].lb --> -Inf
    @fact m.linconstr[2].ub --> 2.0
    @fact m.linconstr[3].lb --> 3.0
    @fact m.linconstr[3].ub --> Inf
    @fact m.quadconstr[1].sense --> :(<=)
end

facts("[macros] Three argument @addConstraint") do
    m = Model()
    @defVar(m, x[1:5])
    @defVar(m, y[2:2:6])

    @addConstraint(m, c, x[4] - y[4] == 1)
    @fact conToStr(m.linconstr[c.idx]) --> "x[4] - y[4] $eq 1"

    @addConstraint(m, d[i=1:5,j=6:-2:2], x[i] - y[j] == 2)
    @fact conToStr(m.linconstr[d[4,4].idx]) --> "x[4] - y[4] $eq 2"

    @addConstraint(m, q[i=1:5], x[i]^2 == 1)
    @fact conToStr(m.quadconstr[q[5].idx]) --> "x[5]² - 1 $eq 0"
end

facts("[macros] @addConstraints") do
    m = Model()
    @defVar(m, x)
    @defVar(m, y[1:3])

    @addConstraints(m, begin
        x + y[1] == 1
        ref[i=1:3], y[1] + y[i] >= i
    end)

    @fact conToStr(m.linconstr[1]) --> "x + y[1] $eq 1"
    @fact conToStr(m.linconstr[2]) --> "2 y[1] $geq 1"
    @fact conToStr(m.linconstr[3]) --> "y[1] + y[2] $geq 2"
    @fact conToStr(m.linconstr[4]) --> "y[1] + y[3] $geq 3"
end

facts("[macros] @addNLConstraints") do
    m = Model()
    @defVar(m, 0 <= x <= 1)
    @defVar(m, y[1:3])
    @setObjective(m, Max, x)

    @addNLConstraints(m, begin
        ref[i=1:3], y[i] == 0
        x + y[1] * y[2] * y[3] <= 0.5
    end)

    @fact length(m.nlpdata.nlconstr) --> 4
    @fact "$(ReverseDiffSparse.base_expression(m.nlpdata.nlconstr[1].terms))" --> "y[i] - 0"
    @fact "$(ReverseDiffSparse.base_expression(m.nlpdata.nlconstr[2].terms))" --> "y[i] - 0"
    @fact "$(ReverseDiffSparse.base_expression(m.nlpdata.nlconstr[3].terms))" --> "y[i] - 0"
    @fact "$(ReverseDiffSparse.base_expression(m.nlpdata.nlconstr[4].terms))" --> "(x + y[1] * y[2] * y[3]) - 0.5"

end

facts("[macros] @setObjective with quadratic") do
    m = Model()
    @defVar(m, x[1:5])
    @setObjective(m, Max, sum{i*x[i]*x[j], i=1:5, j=5:-1:1; isodd(i) && iseven(j)} + 2x[5])

    @fact quadToStr(m.obj) --> "x[1]*x[2] + 3 x[2]*x[3] + x[1]*x[4] + 3 x[3]*x[4] + 5 x[2]*x[5] + 5 x[4]*x[5] + 2 x[5]"
end

facts("[macros] @addConstraint with quadratic") do
    m = Model()
    @defVar(m, x[1:5])

    @addConstraint(m, x[3]*x[1] + sum{x[i]*x[5-i+1], i=1:5; 2 <= i <= 4} + 4x[5] == 1)
    @fact conToStr(m.quadconstr[end]) --> "x[1]*x[3] + x[3]² + 2 x[2]*x[4] + 4 x[5] - 1 $eq 0"

    @addConstraint(m, sum{sum{(x[i] - 2)*x[j],j=4:5},i=2:3} >= -3*x[2]*2*x[4])
    @fact conToStr(m.quadconstr[end]) --> "7 x[2]*x[4] + x[3]*x[4] + x[2]*x[5] + x[3]*x[5] - 4 x[4] - 4 x[5] $geq 0"

    if VERSION > v"0.4.0-"
        @addConstraint(m, sum{x[i],i=1:2}*sum{x[i],i=2:3} >= 0)
        @fact conToStr(m.quadconstr[end]) --> "x[1]*x[2] + x[2]² + x[1]*x[3] + x[2]*x[3] $geq 0"
        @addConstraint(m, x[1]^2 + x[2]*x[3] >= 0)
        @fact conToStr(m.quadconstr[end]) --> "x[1]² + x[2]*x[3] $geq 0"
        @addConstraint(m, x[1]^2 + (x[2]+3)*(x[3]-1) >= 0)
        @fact conToStr(m.quadconstr[end]) --> "x[1]² + x[2]*x[3] + 3 x[3] - x[2] - 3 $geq 0"
        @addConstraint(m, sum{x[i],i=1:2}^2 >= 0)
        @fact conToStr(m.quadconstr[end]) --> "x[1]² + 2 x[1]*x[2] + x[2]² $geq 0"
    end

    myquadexpr = x[1]*x[2]
    @addConstraint(m, sum{i*myquadexpr + x[i], i=1:3} + sum{x[i] + myquadexpr*i, i=1:3} == 0)
    @fact conToStr(m.quadconstr[end]) --> "12 x[1]*x[2] + 2 x[1] + 2 x[2] + 2 x[3] $eq 0"
end

facts("[macros] Triangular indexing, iteration") do
    n = 10
    trimod = Model()
    @defVar(trimod, x[i=1:n,j=i:n])
    @defVar(trimod, y[i=3:2:7,j=-i])
    @fact MathProgBase.numvar(trimod) --> n*(n+1)/2 + 3
    S = Any[(i,i+2) for i in 1:5]
    @defVar(trimod, z[(i,j)=S,k=i:j])
    @fact length(z.tupledict) --> 15
    @addConstraint(trimod, cref[i=1:n,j=i:n], x[i,j] + y[5,-5] == 1)
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
    @defVar(model, x[1:5,2:8,5:6])
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

facts("[macros] @defExpr") do
    model = Model()
    @defVar(model, x[1:3,1:3])
    @defExpr(expr, sum{i*x[i,j] + j, i=1:3,j=1:3})
    @fact affToStr(expr) --> "x[1,1] + x[1,2] + x[1,3] + 2 x[2,1] + 2 x[2,2] + 2 x[2,3] + 3 x[3,1] + 3 x[3,2] + 3 x[3,3] + 18"

    @fact_throws @defExpr(blah[i=1:3], x[i,1]^2)

    @defExpr(y[i=1:2], sum{x[i,1]; i == 1})
    @fact affToStr(y[1]) --> "x[1,1]"
    @fact affToStr(y[2]) --> "0"
end

if VERSION >= v"0.4-"
    facts("[macros] Conditions in constraint indexing") do
        model = Model()
        @defVar(model, x[1:10])
        @addConstraint(model, c1[i=1:9;isodd(i)], x[i] + x[i+1] <= 1)
        @addNLConstraint(model, c2[i=["red","blue","green"], k=9:-2:2; (i == "red" && isodd(k)) || (k >=4 && (i == "blue" || i == "green"))], x[k]^3 <= 1)
        @fact length(model.linconstr) --> 5
        @fact conToStr(model.linconstr[1]) --> "x[1] + x[2] $leq 1"
        @fact conToStr(model.linconstr[2]) --> "x[3] + x[4] $leq 1"
        @fact conToStr(model.linconstr[3]) --> "x[5] + x[6] $leq 1"
        @fact conToStr(model.linconstr[4]) --> "x[7] + x[8] $leq 1"
        @fact conToStr(model.linconstr[5]) --> "x[9] + x[10] $leq 1"
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
end

facts("[macros] Norm parsing") do
    model = Model()
    @defVar(model, x[1:2,1:2])
    @addConstraint(model, -2norm2{x[i,j], i=1:2, j=1:2} + x[1,2] >= -1)
    @addConstraint(model, -2norm2{x[i,j], i=1:2, j=1:2; iseven(i+j)} + x[1,2] >= -1)
    @addConstraint(model, 1 >= 2*norm2{x[i,1],i=1:2})
    @fact conToStr(model.socconstr[1]) --> "2.0 √(x[1,1]² + x[1,2]² + x[2,1]² + x[2,2]²) $leq x[1,2] + 1"
    @fact conToStr(model.socconstr[2]) --> "2.0 √(x[1,1]² + x[2,2]²) $leq x[1,2] + 1"
    @fact conToStr(model.socconstr[3]) --> "2.0 √(x[1,1]² + x[2,1]²) $leq 1"
    @fact_throws @addConstraint(model, (x[1,1]+1)*norm2{x[i,j], i=1:2, j=1:2} + x[1,2] >= -1)
    @fact_throws @addConstraint(model, norm2{x[i,j], i=1:2, j=1:2} + x[1,2] >= -1)
end

facts("[macros] Extraneous terms in QuadExpr (#535)") do
    model = Model()
    @defVar(model, x)
    @defVar(model, y)
    @addConstraint(model, x*x <= y*y)
    @fact conToStr(model.quadconstr[1]) --> "x² - y² $leq 0"
end

if VERSION > v"0.4-"
facts("[macros] Special-case binary multiplication in addToExpression_reorder (#537)") do
    dual = Model()
    @defVar(dual, α[1:3] <= 0)
    @defVar(dual, γ >= 0)
    @addConstraint(dual, ones(3)*γ .<= α)
    @fact conToStr(dual.linconstr[1]) --> "γ - α[1] $leq 0"
    @fact conToStr(dual.linconstr[2]) --> "γ - α[2] $leq 0"
    @fact conToStr(dual.linconstr[3]) --> "γ - α[3] $leq 0"
end
end
