import Compat
using JuMP

function test(n::Int)
    m = Model()
    a = rand(n,n)
    @variable(m,x[1:n,1:n])
    b = rand(n,n,n)
    @variable(m,y[1:n,1:n,1:n])
    c = rand(n)
    @variable(m,z[1:n])

    #initialize
    @constraint(m,sum(c[i]*z[i] for i=1:n)<=0)

    #Vector
    elapsed = @elapsed @constraint(m,sum(c[i]*z[i] for i=1:n)<=1)
    println("Vector with sum(): $(elapsed)")

    elapsed = @elapsed @constraint(m,Compat.dot(c,z) <= 1)
    println("Vector with vecdot() : $(elapsed)")

    #2D Matrix
    elapsed = @elapsed @constraint(m,sum(a[i,j]*x[i,j] for i=1:n,j=1:n)<=1)
    println("2D Matrix with sum(): $(elapsed)")

    elapsed = @elapsed @constraint(m,Compat.dot(a,x)<=1)
    println("2D Matrix with bigvecdot(): $(elapsed)")

    #3D Matrix
    elapsed = @elapsed @constraint(m,sum(b[i,j,k]*y[i,j,k] for i=1:n,j=1:n,k=1:n)<=1)
    println("3D Matrix with sum(): $(elapsed)")

    elapsed = @elapsed @constraint(m,Compat.dot(b,y)<=1)
    println("3D Matrix with vecdot(): $(elapsed)")
    return 0
end

for n in [10, 50, 100, 200, 300, 400]
    println("n = $n")
    test(n)
end
