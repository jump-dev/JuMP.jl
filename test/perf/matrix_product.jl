using LinearAlgebra
using JuMP

function test(n::Int, verbose::Bool)
    model = Model()
    a = rand(n, n)
    @variable(model, x[1:n, 1:n])

    # Compile it
    x * a
    a * x * a

    elapsed = @elapsed x * a
    if verbose
        println("2D Matrix product `x * a`: $(elapsed)")
    end

    elapsed = @elapsed a * x * a
    if verbose
        println("2D Matrix product `a * x * a`: $(elapsed)")
    end

    return
end

# Compile the methods needed.
test(2, false)

for n in [10, 20, 50, 100]
    println("n = $n")
    test(n, true)
end
