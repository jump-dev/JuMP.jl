using JuMP
m = Model()
n = 100
a = rand(n,n)
@defVar(m,x[1:n,1:n])
b = rand(n,n,n)
@defVar(m,y[1:n,1:n,1:n])
c = rand(n)
@defVar(m,z[1:n])

#initialize
@addConstraint(m,sum{c[i]*z[i],i=1:n}<=0)

#Vector
tic()
@addConstraint(m,sum{c[i]*z[i],i=1:n}<=1)
println("Vector with sum{}: $(toq())")
tic()
@addConstraint(m,c*z <= 1)
println("Vector with * : $(toq())")

#2D Matrix
tic()
@addConstraint(m,sum{a[i,j]*x[i,j],i=1:n,j=1:n}<=1)
println("2D Matrix with sum{}: $(toq())")
tic()
@addConstraint(m,bigdot(a,x)<=1)
println("2D Matrix with bigdot(): $(toq())")

#3D Matrix
tic()
@addConstraint(m,sum{b[i,j,k]*y[i,j,k],i=1:n,j=1:n,k=1:n}<=1)
println("3D Matrix with sum{}: $(toq())")
tic()
@addConstraint(m,bigdot(b,y)<=1)
println("3D Matrix with bigdot(): $(toq())")

