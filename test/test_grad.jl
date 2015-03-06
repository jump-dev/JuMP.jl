using Base.Test
using Base.Meta
using ReverseDiffSparse

x = placeholders(3)

ex = @processNLExpr sin(x[1])
# genfgrad(ex) # examine symbolic function
fg = genfgrad_simple(ex)

out = zeros(3)
xvals = [1.1,2.2,3.3]
fval = fg(xvals, out)
@test_approx_eq fval sin(xvals[1])
@test_approx_eq out [cos(xvals[1]),0.0,0.0]

ex = @processNLExpr sin(x[1])^2
fg = genfgrad_simple(ex)
fval = fg(xvals, out)
@test_approx_eq fval sin(xvals[1])^2
@test_approx_eq out[1] sin(2xvals[1])

ex = @processNLExpr x[1]*x[1]
fg = genfgrad_simple(ex)
fval = fg(xvals, out)
@test_approx_eq fval xvals[1]^2
@test_approx_eq out[1] 2*xvals[1]

ex = @processNLExpr x[1]/x[2]
fg = genfgrad_simple(ex)
fval = fg(xvals, out)
@test_approx_eq fval xvals[1]/xvals[2]
@test_approx_eq out[1] 1/xvals[2]
@test_approx_eq out[2] -xvals[1]/xvals[2]^2

ex = @processNLExpr exp(sin(x[1]*x[2]/x[3])) 
fg = genfgrad_simple(ex)
xvals = [3.4,2.1,6.7]
fval = fg(xvals, out)
q = xvals[1]*xvals[2]/xvals[3] 
@test_approx_eq fval exp(sin(q)) 
@test_approx_eq out[1] xvals[2]*cos(q)*exp(sin(q))/xvals[3] 
@test_approx_eq out[2] xvals[1]*cos(q)*exp(sin(q))/xvals[3] 
@test_approx_eq out[3] -xvals[1]*xvals[2]*cos(q)*exp(sin(q))/xvals[3]^2

ex = @processNLExpr exp(sin(x[1]*x[2]+x[3]^2)) + 2x[1]*x[1]
fg = genfgrad_simple(ex)
xvals = [3.4,2.1,6.7]
fval = fg(xvals, out)
q = xvals[1]*xvals[2]+xvals[3]^2
@test_approx_eq fval exp(sin(q)) + 2xvals[1]^2
@test_approx_eq out[1] xvals[2]*cos(q)*exp(sin(q)) + 4xvals[1]
@test_approx_eq out[2] xvals[1]*cos(q)*exp(sin(q))
@test_approx_eq out[3] 2xvals[3]*cos(q)*exp(sin(q))

# test reusing the same function
xvals = [35.2,-1.2,3.9]
fval = fg(xvals, out)
q = xvals[1]*xvals[2]+xvals[3]^2
@test_approx_eq fval exp(sin(q)) + 2xvals[1]^2
@test_approx_eq out[1] xvals[2]*cos(q)*exp(sin(q)) + 4xvals[1]
@test_approx_eq out[2] xvals[1]*cos(q)*exp(sin(q))
@test_approx_eq out[3] 2xvals[3]*cos(q)*exp(sin(q))

ex = @processNLExpr exp(sin(x[1]*x[2]+x[3]^2)) - 2x[1]*x[1]
fg = genfgrad_simple(ex)
xvals = [3.4,2.1,6.7]
fval = fg(xvals, out)
q = xvals[1]*xvals[2]+xvals[3]^2
@test_approx_eq fval exp(sin(q)) - 2xvals[1]^2
@test_approx_eq out[1] xvals[2]*cos(q)*exp(sin(q)) - 4xvals[1]
@test_approx_eq out[2] xvals[1]*cos(q)*exp(sin(q))
@test_approx_eq out[3] 2xvals[3]*cos(q)*exp(sin(q))

# other variables present
y = 10.0
z = 2
ex = @processNLExpr z*x[1]^y
fg = genfgrad_simple(ex)
fval = fg(xvals,out)
@test_approx_eq fval z*xvals[1]^y
@test_approx_eq out[1] y*z*xvals[1]^(y-1)

# sum syntax
x = placeholders(10)
out = zeros(10)
ex = @processNLExpr sum{3x[i], i = 1:10}
fg = genfgrad_simple(ex)
xvals = rand(10)
fval = fg(xvals, out)
@test_approx_eq fval 3*sum(xvals)
@test_approx_eq out fill(3.0,10)

# with conditions
out = zeros(10)
ex = @processNLExpr sum{3x[i], i = 1:10; i > 3}
fg = genfgrad_simple(ex)
fval = fg(xvals, out)
@test_approx_eq fval 3*sum(xvals[4:end])
@test_approx_eq out [0.0,0.0,0.0,fill(3.0,7)]

vars = placeholders(20)
x = vars[1:10]
y = vars[11:20]
out = zeros(20)
ex = @processNLExpr sum{x[i]*y[i], i = 1:10} + sin(x[1])
fg = genfgrad_simple(ex)
vals = rand(20)
fval = fg(vals, out)
@test_approx_eq fval sum([vals[i]*vals[i+10] for i in 1:10]) + sin(vals[1])
deriv = [vals[11:20],vals[1:10]]
deriv[1] += cos(vals[1])
@test_approx_eq out deriv


# prod syntax
x = placeholders(5)
ex = @processNLExpr prod{x[i], i = 1:5}
fg = genfgrad_simple(ex)
xvals = rand(5)
out = zeros(5)
fval = fg(xvals, out)
@test_approx_eq fval prod(xvals)
@test_approx_eq out prod(xvals)./xvals

# nested sums
S = Array(Vector{Int},3)
S[1] = [1,2,3]
S[2] = [3,4,5]
S[3] = [5,6,7]
x = placeholders(7)
ex = @processNLExpr sum{ sum{ x[k], k in S[i] }, i in 1:3 }
fg = genfgrad_simple(ex)
xvals = rand(7)
out = zeros(7)
fval = fg(xvals, out)
@test_approx_eq out [1.0,1.0,2.0,1.0,2.0,1.0,1.0]

R = [x[1],x[2]]
ex = @processNLExpr sum{ z, z in R }
out = zeros(7)
fg = genfgrad_simple(ex)
fval = fg(xvals, out)
@test_approx_eq out [1.0,1.0,0.0,0.0,0.0,0.0,0.0]

# dot syntax
type MyType
    x
    y
end
t = MyType([1.0,2.0],x)

ex = @processNLExpr sum{ t.x[i]*t.y[i], i = 1:2 }
out = zeros(2)
fg = genfgrad_simple(ex)
fval = fg(xvals[1:2], out)
@test_approx_eq fval xvals[1]+2xvals[2]
@test_approx_eq out [1.0,2.0]

t = [MyType(1.0,2.0),MyType(2.0,3.0)]
ex = @processNLExpr sum{ t[i].x*x[i], i = 1:2 }
out = zeros(2)
fg = genfgrad_simple(ex)
fval = fg(xvals[1:2], out)
@test_approx_eq fval dot([1.0,2.0],xvals[1:2])
@test_approx_eq out [1.0,2.0]

# ref syntax
x = placeholders(2)
k = 2
ex = @processNLExpr 3x[k-1] + x[k+2-2]
fg = genfgrad_simple(ex)
fval = fg(xvals[1:2], out)
@test_approx_eq fval 3xvals[1] + xvals[2]
@test_approx_eq out [3.0,1.0]

# special variable names
T = 10
ex = @processNLExpr x[1]+T
fg = genfgrad_simple(ex)
fval = fg(xvals[1:2], out)

# very simple expressions
x = placeholders(1)[1]
ex = @processNLExpr x
fg = genfgrad_simple(ex)
fval = fg([2.0], out)
@test_approx_eq fval 2.0
@test_approx_eq out[1] 1.0

# expanded indices using tuples
x = placeholders(2)
I = [(1,2)]
ex = @processNLExpr sum{ x[i]*x[j], (i,j) in I }
fg = genfgrad_simple(ex)
fval = fg([2.0,3.0], out)
@test_approx_eq fval 6.0
@test_approx_eq out[1] 3.0
@test_approx_eq out[2] 2.0

# dummy sum
ex = @processNLExpr x[1] + sum{ 2, i in 1:10 }
fg = genfgrad_simple(ex)
fval = fg([2.0], out)
@test_approx_eq fval 22.0
@test_approx_eq out[1] 1.0

# exponents
y = 2
ex = @processNLExpr x[1]^y
fg = genfgrad_simple(ex)
fval = fg([-3.0],out)
@test_approx_eq fval (-3)^2
@test_approx_eq out[1] 2*-3

y = -2
ex = @processNLExpr y^x[1]
fg = genfgrad_simple(ex)
fval = fg([0.3],out)
@test isnan(fval)


# zeros in products
ex = @processNLExpr prod{ x[i], i = 1:2 }
fg = genfgrad_simple(ex)
fval = fg([2.0,0.0],out)
@test_approx_eq fval 0.0
@test_approx_eq out[1] 0.0
@test_approx_eq out[2] 2.0

ex = @processNLExpr prod{ x[i], i = 1:2 }
fg = genfgrad_simple(ex)
fval = fg([0.0,0.0],out)
@test_approx_eq fval 0.0
@test_approx_eq out[1] 0.0
@test_approx_eq out[2] 0.0

# ifelse
ex = @processNLExpr x[1]*ifelse(x[1] >= x[2], x[1],x[2])
fg = genfgrad_simple(ex)
fval = fg([2.5,1.0],out)
@test_approx_eq fval 2.5^2
@test_approx_eq out[1] 2*2.5
@test_approx_eq out[2] 0.0
fval = fg([0.5,1.0],out)
@test_approx_eq fval 0.5
@test_approx_eq out[1] 1.0
@test_approx_eq out[2] 0.5

ex = @processNLExpr x[1]*ifelse(x[1] >= x[2] || true, x[1],x[2])
fg = genfgrad_simple(ex)
fval = fg([0.5,1.0],out)
@test_approx_eq fval 0.25

# errors from duplicate index variables
@test isexpr(macroexpand(:(@processNLExpr sum{sum{x[i],i=1:2},i=1:2})),:error)
@test isexpr(macroexpand(:(@processNLExpr sum{sum{x[i],(i,j)=A},i=1:2})),:error)
@test isexpr(macroexpand(:(@processNLExpr sum{sum{x[i],i=A},(i,j)=B})),:error)

# nested subexpressions
c = 10
ex = @parametricExpr c*x[1]^2
ex2 = @processNLExpr ex+x[2]
fg = genfgrad_simple(ex2)
fval = fg([2.5,1.0],out)
@test_approx_eq fval c*2.5^2+1
@test_approx_eq out[1] 2*c*2.5
@test_approx_eq out[2] 1

@test_approx_eq ReverseDiffSparse.getvalue(ex, [2.0,3.0]) c*2.0^2

ex = @parametricExpr j c*x[j]^2
ex2 = @processNLExpr ex[1]+x[2]
fg = genfgrad_simple(ex2)
fval = fg([2.5,1.0],out)
@test_approx_eq fval c*2.5^2+1
@test_approx_eq out[1] 2*c*2.5
@test_approx_eq out[2] 1

@test_approx_eq ReverseDiffSparse.getvalue(ex[1], [2.0,3.0]) c*2.0^2
@test_approx_eq ReverseDiffSparse.getvalue(ex[2], [2.0,3.0]) c*3.0^2

ex2 = @processNLExpr sum{ex[i],i=1:1}+x[2]
fg = genfgrad_simple(ex2)
fval = fg([2.5,1.0],out)
@test_approx_eq fval c*2.5^2+1
@test_approx_eq out[1] 2*c*2.5
@test_approx_eq out[2] 1

ex = @parametricExpr i j ifelse(i==j,1,0)
ex2 = @processNLExpr ex[1,1] + ex[1,2]
f = genfval_simple(ex2)
@test_approx_eq f([3.0]) 1.0

ex2 = @processNLExpr ifelse(ex[1,1] == ex[2,2],3,0)
f = genfval_simple(ex2)
@test_approx_eq f([3.0]) 3.0

ex = @parametricExpr i sum{x[k],k=1:2;k==i}
ex2 = @processNLExpr ex[1]
fg = genfgrad_simple(ex2)
fval = fg([2.0,3.0],out)
@test_approx_eq fval 2.0
@test_approx_eq out[1] 1.0
@test_approx_eq out[2] 0.0

# embedded indices
idx = [1]
ex = @processNLExpr sum{x[idx[j]], j = 1:1}
fg = genfgrad_simple(ex)
fval = fg([2.5,1.0],out)
@test_approx_eq fval 2.5
@test_approx_eq out[1] 1
@test_approx_eq out[2] 0

# JuMP issue #406
first = @parametricExpr x[1]
second = @parametricExpr x[1] + first
ex = @processNLExpr first + second
fg = genfgrad_simple(ex)
fval = fg([3.0,1.0], out)
@test_approx_eq fval 9.0
@test_approx_eq out[1] 3
@test_approx_eq out[2] 0

# embedded tuples
d = Dict()
d[(1,1)] = [1,2]
ex = @processNLExpr sum{x[i], i = d[(1,1)]}
fg = genfgrad_simple(ex)
fval = fg([2.5,1.0],out)
@test_approx_eq fval 3.5
@test_approx_eq out[1] 1
@test_approx_eq out[2] 1

ex = @parametricExpr x[1] + x[1]^2 + x[1]*x[2] + x[2]^2
ex2 = @processNLExpr ex - 1
fg = genfgrad_simple(ex2)
fval = fg([1.3,2.4],out)
@test_approx_eq fval 1.3 + 1.3^2 + 1.3*2.4 + 2.4^2 - 1.0
@test_approx_eq out[1] 1+2*1.3+2.4
@test_approx_eq out[2] 1.3+2*2.4


println("Passed tests")
