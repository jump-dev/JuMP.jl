# Compile speed tests
# Tests time to get from nothing to certain tasks
using JuMP  # Ensure compiled before running the tests

# Test 1
# Size 100, aff_str into model_str
test1_time = @elapsed begin
    run(`$(ENV["_"]) -e """
    using JuMP
    m = Model()
    N = 100
    @variable(m, x[1:N])
    @constraint(m, sum(i*x[i] for i=1:N) >= N)
    @time JuMP.aff_str(JuMP.REPLMode, m.linconstr[end].terms)
    @time JuMP.model_str(JuMP.REPLMode, m)"""
    `)
end

# Test 2
# Size 100, model_str
test2_time = @elapsed begin
    run(`$(ENV["_"]) -e """
    using JuMP
    m = Model()
    N = 100
    @variable(m, x[1:N])
    @constraint(m, sum(i*x[i] for i=1:N) >= N)
    @time JuMP.model_str(JuMP.REPLMode, m)"""
    `)
end

# Test 3
# Size 10,000, print model
test3_time = @elapsed begin
    run(`$(ENV["_"]) -e """
    using JuMP
    m = Model()
    N = 10000
    @variable(m, x[1:N])
    @constraint(m, sum(i*x[i] for i=1:N) >= N)
    @time sprint(print, m)"""
    `)
end

# Test 4
# Size 10, 4 variables, print model
test4_time = @elapsed begin
    run(`$(ENV["_"]) -e """
    using JuMP
    m = Model()
    N = 10
    @variable(m, x1[1:N])
    @variable(m, x2[1:N,f=1:N])
    @variable(m, x3[1:N,f=1:2:N])
    @variable(m, x4[[:a,:b,:c]])
    @variable(m, x5[[:a,:b,:c],[:d,"e",4]])
    @constraint(m,
        sum(i*x1[i] for i=1:N) +
        sum(i*f*x2[i,f] for i=1:N,f=1:N) +
        sum(i*f*x3[i,f] for i=1:N,f=1:2:N) +
        sum(x4) >= N)
    @time sprint(print, m)"""
    `)
end

println("Julia ", VERSION)
println(
    "Current branch: $(chomp(read(`git rev-parse --abbrev-ref HEAD`, String)))",
)
println("Current commit: $(chomp(read(`git rev-parse HEAD`, String)))")
@printf("Test 1: %7.2f s\n", test1_time)
@printf("Test 2: %7.2f s\n", test2_time)
@printf("Test 3: %7.2f s\n", test3_time)
@printf("Test 4: %7.2f s\n", test4_time)
