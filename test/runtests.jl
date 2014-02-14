
tests = ["test_grad.jl",
        "test_hessian.jl"]

for t in tests
    println("$t:")
    include(t)
end
