# affexpr.jl
# Test coverage for AffExpr

maff = Model("max")
@defVar(maff, 0 <= x[1:5] <= 1)
@defVar(maff, 0 <= LongName <= 99)
a1 = x[1] + LongName + 5
@test affToStr(a1) == "1.0 _col1 + 1.0 LongName + 5.0"
a2 = 2*(x[2] + LongName + x[2]) + 0
println("  TODO: Collect like terms before print")
println("        e.g. $a2")
#@test affToStr(a2) == "2.0 _col2 + 2.0 LongName"
