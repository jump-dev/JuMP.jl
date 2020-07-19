using BenchmarkTools
using JuMP

function model0()
    eval(quote
        m = Model()
        @variable(m, a)
        @variable(m, b)
        @variable(m, c)
        @variable(m, x[1:3])
        @constraint(m, a + b + c <= 1)
        @constraint(m, x[1] + x[2] + x[3] <= 1)
        return m
    end)
end

function model1()
    eval(quote
        m = Model()
        @variable(m, x[1:10])
        @variable(m, y[1:10])
        @constraint(m, sum(x) <= 1)
        @constraint(m, sum(x .* y) <= sum(y))
        return m
    end)
end

suite = BenchmarkGroup()
suite["model0"] = @benchmarkable model0()
suite["model1"] = @benchmarkable model1()
# tune!(suite)
results = run(suite, verbose=true)
