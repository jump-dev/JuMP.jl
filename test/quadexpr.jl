# quadexpr.jl
# Test coverage for QuadExpr

mq = Model("max")
@defVar(mq, 0 <= x[1:5] <= 1)
@defVar(mq, 0 <= LongName <= 99)
println("  TODO: Finish overloads, this is all brittle!")
#a1 = x[1]*x[2] + 27.2*LongName + 5
#@test quadToStr(a1) == "1.0 _col1*_col2 + 27.2 LongName + 5.0"
#a2 = x[1]*x[2] + x[2]*x[1]
println("  TODO: Collect like terms before print")
#println("        e.g. $a2")
#@test quadToStr(a2) == "2.0 _col1*_col2"
