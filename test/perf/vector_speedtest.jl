using LinearAlgebra
using JuMP

function test(n::Int)
    model = Model()
    a = rand(n, n)
    @variable(model, x[1:n, 1:n])
    b = rand(n, n, n)
    @variable(model, y[1:n, 1:n, 1:n])
    c = rand(n)
    @variable(model, z[1:n])

    #initialize
    @constraint(model, sum(c[i] * z[i] for i in 1:n) <= 0)

    #Vector
    elapsed = @elapsed @constraint(model, sum(c[i] * z[i] for i in 1:n) <= 1)
    println("Vector with sum(): $(elapsed)")

    elapsed = @elapsed @constraint(model, dot(c, z) <= 1)
    println("Vector with vecdot() : $(elapsed)")

    #2D Matrix
    elapsed = @elapsed @constraint(
        model,
        sum(a[i, j] * x[i, j] for i in 1:n, j in 1:n) <= 1
    )
    println("2D Matrix with sum(): $(elapsed)")

    elapsed = @elapsed @constraint(model, dot(a, x) <= 1)
    println("2D Matrix with bigvecdot(): $(elapsed)")

    #3D Matrix
    elapsed = @elapsed @constraint(
        model,
        sum(b[i, j, k] * y[i, j, k] for i in 1:n, j in 1:n, k in 1:n) <= 1
    )
    println("3D Matrix with sum(): $(elapsed)")

    elapsed = @elapsed @constraint(model, dot(b, y) <= 1)
    println("3D Matrix with vecdot(): $(elapsed)")
    return 0
end

for n in [10, 50, 100, 200, 300]
    println("n = $n")
    test(n)
end
