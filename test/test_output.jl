using Base.Test
using ReverseDiffSparse

y = placeholders(3)

ex = @processNLExpr sin(y[1])
@test to_flat_expr(ex) == :(sin(x[1]))

@test base_expression(ex) == :(sin(y[1]))

ex = @parametricExpr i sin(y[i])
@test base_expression(ex) == :(sin(y[i]))


ex = @processNLExpr sin(y[1])^2
@test to_flat_expr(ex) == :(sin(x[1])^2)


ex = @processNLExpr y[1]/y[2]
@test to_flat_expr(ex) == :(x[1]/x[2])

q = 10.0
z = 2
ex = @processNLExpr z*y[1]^q
@test to_flat_expr(ex) == :(2*x[1]^10.0)


x = placeholders(5)
ex = @processNLExpr sum{3x[i], i = 1:5}
@test to_flat_expr(ex) == :(3x[1] + 3x[2] + 3x[3] + 3x[4] + 3x[5])

ex = @processNLExpr sum{3x[i], i = 1:5; iseven(i)}
@test to_flat_expr(ex) == :(3x[2] + 3x[4])

ex = @processNLExpr prod{x[i], i = 1:5}
@test to_flat_expr(ex) == :(x[1]*x[2]*x[3]*x[4]*x[5])

# empty sums and products
S = []
ex = @processNLExpr x[1] + prod{x[i], i in S} + sum{x[i], i in S}
@test to_flat_expr(ex) == :(x[1] + 1 + 0)

# singleton sums and products
S = [2]
ex = @processNLExpr x[1] + prod{x[i], i in S} + sum{x[i], i in S}
@test to_flat_expr(ex) == :(x[1] + x[2] + x[2])

# ifelse
ex = @processNLExpr x[1]*ifelse(x[1] >= x[2], x[1],x[2])
@test to_flat_expr(ex) == :(x[1] * ifelse(x[1] >= x[2],x[1],x[2]))

ex = @processNLExpr x[1]*ifelse(x[1] >= x[2] && x[1]^2 < 5, x[1],x[2])
@test to_flat_expr(ex) == :(x[1] * ifelse(x[1] >= x[2] && x[1]^2 < 5,x[1],x[2]))
