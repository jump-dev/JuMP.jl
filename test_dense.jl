# test_dense.jl
# Make a big dense matrix
require("julp.jl")

using Julp

function doTest()
  srand(10)
  N = 10000
  M = 1000
  A = randi(    10,M,N)
  b = N+randi(N*10,M)
  c = randi(    10,N)

  tic()
  m = Model("max")
  v = [ Variable(m,0,1,0,"x$j") for j=1:N ]

  #m.objective = Julp.AffExpr([(v[j],convert(Float64,c[j])) for j=1:N])
  #m.objective = sum([ c[j]*v[j] for j=1:N ])
  #m.objective = Julp.lpSum([ v[j]*c[j] for j=1:N ])
  m.objective = @SumExpr([c[j]*v[j] for j=1:N])

  for i = 1:M
    #lhs = sum([v[j]*A[i,j] for j=1:N])
    #lhs = Julp.AffExpr([(v[j],convert(Float64,A[i,j])) for j=1:N])
    #lhs = Julp.lpSum([ v[j] * A[i,j] for j = 1:N ])
    lhs = @SumExpr([ A[i,j]*v[j] for j = 1:N ])
    addConstraint(m, lhs <= b[i])
  end
  toc()
  println("In building the model in memory")
  tic()
  writeLP(m,"dense.lp")
  toc()
  println("Total")
end

doTest()
