# expr.jl
# Test coverage for AffExpr and QuadExpr

maff = Model()
@defVar(maff, 0 <= x[1:5] <= 1)
@defVar(maff, 0 <= LongName <= 99)
# Test affToStr
a1 = x[1] + LongName + 5
@test affToStr(a1) == "x[1] + LongName + 5"
# Test like term collection
a2 = 2*(x[2] + LongName + x[2]) + 0
@test affToStr(a2) == "4 x[2] + 2 LongName"
# Test appending functionality
push!(a1, 5.0, x[2])
@test affToStr(a1) == "x[1] + LongName + 5 x[2] + 5"
append!(a1, a2)
@test affToStr(a1) == "x[1] + 3 LongName + 9 x[2] + 5"

# Test quadToStr
q1 = x[1]*x[2] + 27.2*LongName + 5
@test quadToStr(q1) == "x[1]*x[2] + 27.2 LongName + 5"
# Test like term collection
q2 = x[1]*x[2] + x[2]*x[1]
@test quadToStr(q2) == "2 x[1]*x[2]"

# Test getValue(AffExpr), getValue(QuadExpr)
let
    m = Model()
    @defVar(m, 1 <= x[1:3] <= 2)
    @addConstraint(m, x[1] + x[3] == 3)
    @setObjective(m, Max, x[2]+x[3])
    stat = solve(m)
    @test stat == :Optimal
    @test getValue(x[1]-x[2]+2x[3]-1.0) == 2.0
    @test getValue(x[1]*x[1]-2x[2]*x[1]+3x[2]+1) == 4.0
end
