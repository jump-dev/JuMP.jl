using JuMP

function test(n::Int)
    m = Model()
    a = rand(n,n)
    @defVar(m,x[1:n,1:n])
    b = rand(n,n,n)
    @defVar(m,y[1:n,1:n,1:n])
    c = rand(n)
    @defVar(m,z[1:n])

    #initialize
    @constraint(m,sum{c[i]*z[i],i=1:n}<=0)

    #Vector
    tic()
    @constraint(m,sum{c[i]*z[i],i=1:n}<=1)
    println("Vector with sum{}: $(toq())")
    tic()
    @constraint(m,dot(c,z) <= 1)
    println("Vector with dot() : $(toq())")

    #2D Matrix
    tic()
    @constraint(m,sum{a[i,j]*x[i,j],i=1:n,j=1:n}<=1)
    println("2D Matrix with sum{}: $(toq())")
    tic()
    @constraint(m,dot(a,x)<=1)
    println("2D Matrix with bigdot(): $(toq())")

    #3D Matrix
    tic()
    @constraint(m,sum{b[i,j,k]*y[i,j,k],i=1:n,j=1:n,k=1:n}<=1)
    println("3D Matrix with sum{}: $(toq())")
    tic()
    @constraint(m,dot(b,y)<=1)
    println("3D Matrix with dot(): $(toq())")
    return 0
end

for n in [10, 50, 100, 200, 300, 400]
    println("n = $n")
    test(n)
end
