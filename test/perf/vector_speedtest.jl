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
    @constraint(m,sum{c[i]*z[i],i=1:n}<=0)

    #Vector
    tic()
    @constraint(m,sum{c[i]*z[i],i=1:n}<=1)
    println("Vector with sum{}: $(toq())")
    tic()
    @constraint(m,vecdot(c,z) <= 1)
    println("Vector with vecdot() : $(toq())")

    #2D Matrix
    tic()
    @constraint(m,sum{a[i,j]*x[i,j],i=1:n,j=1:n}<=1)
    println("2D Matrix with sum{}: $(toq())")
    tic()
    @constraint(m,vecdot(a,x)<=1)
    println("2D Matrix with bigvecdot(): $(toq())")

    #3D Matrix
    tic()
    @constraint(m,sum{b[i,j,k]*y[i,j,k],i=1:n,j=1:n,k=1:n}<=1)
    println("3D Matrix with sum{}: $(toq())")
    tic()
    @constraint(m,vecdot(b,y)<=1)
    println("3D Matrix with vecdot(): $(toq())")
    return 0
end

for n in [10, 50, 100, 200, 300, 400]
    println("n = $n")
    test(n)
end
