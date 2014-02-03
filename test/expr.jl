# expr.jl
# Test coverage for AffExpr and QuadExpr

maff = Model()
@defVar(maff, 0 <= x[1:5] <= 1)
@defVar(maff, 0 <= LongName <= 99)
# Test affToStr
a1 = x[1] + LongName + 5
@test affToStr(a1) == "1.0 x[1] + 1.0 LongName + 5.0"
# Test like term collection
a2 = 2*(x[2] + LongName + x[2]) + 0
@test affToStr(a2) == "4.0 x[2] + 2.0 LongName"

# Test quadToStr
q1 = x[1]*x[2] + 27.2*LongName + 5
@test quadToStr(q1) == "1.0 x[1]*x[2] + 27.2 LongName + 5.0"
# Test like term collection
q2 = x[1]*x[2] + x[2]*x[1]
@test quadToStr(q2) == "2.0 x[1]*x[2]"
