# test_dense.jl
# Make a big dense matrix
require("julp.jl")

# Failed attempt at a macro
macro lpsum(coefs,vars,counter, range)
  :(Julp.AffExpr([($vars,$coefs) for $counter in $range]))
end

function doTest()
  N = 1000
  M = 100
  A = randi(    10,M,N)
  b = N+randi(N*10,M)
  c = randi(    10,N)

  tic()
  m = Julp.Model("max")
  v = [ Julp.Variable(m,"x$j",0,1) for j=1:N ]

  m.objective = Julp.AffExpr([(v[j],convert(Float64,c[j])) for j=1:N])
  #m.objective = sum([ c[j]*v[j] for j=1:N ])
  

  for i = 1:M
    #lhs = sum([v[j]*A[i,j] for j=1:N])
    lhs = Julp.AffExpr([(v[j],convert(Float64,A[i,j])) for j=1:N])
    Julp.AddConstraint(m, lhs <= b[i])
  end
  toc()
  println("In building the model in memory")
  tic()
  Julp.WriteLP(m,"dense.lp")
  toc()
  println("Total")
end

doTest()
