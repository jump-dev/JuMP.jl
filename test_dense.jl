# test_dense.jl
# Make a big dense matrix
require("julp.jl")

N = 1000
M = 100
A = randi(    10,M,N)
b = N+randi(N*10,M)
c = randi(    10,N)

tic()
m = Julp.Model("max")
v = Array(Julp.Variable,0)
for j = 1:N
  push(v,Julp.Variable(m,"x$j",0,1))
end

# Need a default constructor for AffExpr!
m.objective = c[1]*v[1]
for j = 2:N
  m.objective = m.objective + c[j]*v[j]
end

for i = 1:M
  # Need a default constructor for AffExpr!
  lhs = Julp.AffExpr([(v[1],convert(Float64,A[i,1]))],  0.0)
  for j = 1:N
    # Need to overload +=
    lhs = lhs + v[j] * A[i,j]
  end
  Julp.AddConstraint(m, lhs <= b[i])
end
Julp.WriteLP(m,"dense.lp")
toc()
