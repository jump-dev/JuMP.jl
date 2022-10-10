module NonlinearBenchmark

using JuMP
import BenchmarkTools
import InfiniteOpt
import Ipopt
import Random
import Symbolics

function benchmark_group()
    lookup = Dict(
        "perf_nl_" => "@NL",
        "perf_nlexpr_" => "NonlinearExpr",
        "perf_infopt_" => "InfiniteOpt",
        "perf_symbolics_" => "Symbolics",
    )
    suite = BenchmarkTools.BenchmarkGroup()
    for v in values(lookup)
        suite[v] = BenchmarkTools.BenchmarkGroup()
    end
    for name in names(@__MODULE__; all = true)
        f = getfield(@__MODULE__, name)
        for (k, v) in lookup
            if startswith("$name", k)
                fname = replace("$name", k => "")
                suite[v][fname] = BenchmarkTools.@benchmarkable $f()
                break
            end
        end
    end
    return suite
end

function runbenchmarks()
    suite = benchmark_group()
    return BenchmarkTools.run(suite)
end

# sum
#
# nlexpr is slower because it builds up the product via operator overloading,
# creating a lot of temporary objects. @NL gets to see the full +(args...) to it
# builds the expression in-place.
#
# We could fix this by implementing a n-argy method for +, but that gets
# difficult with method ambiguities.

function perf_nl_sum()
    model = Model()
    @variable(model, x)
    @NLobjective(model, Min, sum(x^i for i in 1:10_000))
    return
end

function perf_nlexpr_sum()
    model = Model()
    @variable(model, x)
    @objective(model, Min, sum(x^Float64(i) for i in 1:10_000))
    return
end

function perf_infopt_sum()
    model = InfiniteOpt.InfiniteModel()
    @variable(model, x)
    @objective(model, Min, sum(x^i for i in 1:10_000))
    return
end

function perf_symbolics_sum()
    Symbolics.@variables x
    sum(x^i for i in 1:10_000)
    return
end

# prod
#
# nlexpr is slower because it builds up the product via operator overloading,
# creating a lot of temporary objects. @NL gets to see the full *(args...) to it
# builds the expression in-place.
#
# We could fix this by implementing a n-argy method for *, but that gets
# difficult with method ambiguities.

function perf_nl_prod()
    model = Model()
    @variable(model, x)
    @NLobjective(model, Min, prod(x^i for i in 1:10_000))
    return
end

function perf_nlexpr_prod()
    model = Model()
    @variable(model, x)
    @objective(model, Min, prod(x^Float64(i) for i in 1:10_000))
    return
end

function perf_infopt_prod()
    model = InfiniteOpt.InfiniteModel()
    @variable(model, x)
    @objective(model, Min, prod(x^i for i in 1:10_000))
    return
end

function perf_symbolics_prod()
    Symbolics.@variables x
    prod(x^i for i in 1:10_000)
    return
end

# many_constraints

function perf_nl_many_constraints()
    model = Model()
    @variable(model, x[1:10_000])
    @NLconstraint(model, [i = 1:10_000], sin(x[i]) <= cos(i))
    return
end

function perf_nlexpr_many_constraints()
    model = Model()
    @variable(model, x[1:10_000])
    @constraint(model, [i = 1:10_000], sin(x[i]) <= cos(i))
    return
end

function perf_infopt_many_constraints()
    model = InfiniteOpt.InfiniteModel()
    @variable(model, x[1:10_000])
    @constraint(model, [i = 1:10_000], sin(x[i]) <= cos(i))
    return
end

function perf_symbolics_many_constraints()
    Symbolics.@variables x[1:10_000]
    [sin(x[i]) - cos(i) for i in 1:10_000]
    return
end

# mle

function perf_nl_mle()
    Random.seed!(1234)
    n = 1_000
    data = randn(n)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, μ, start = 0.0)
    @variable(model, σ >= 0.0, start = 1.0)
    @NLobjective(
        model,
        Max,
        n / 2 * log(1 / (2 * π * σ^2)) -
        sum((data[i] - μ)^2 for i in 1:n) / (2 * σ^2)
    )
    optimize!(model)
    return
end

function perf_nlexpr_mle()
    Random.seed!(1234)
    n = 1_000
    data = randn(n)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, μ, start = 0.0)
    @variable(model, σ >= 0.0, start = 1.0)
    @objective(
        model,
        Max,
        n / 2 * log(1 / (2 * π * σ^2)) -
        sum((data[i] - μ)^2 for i in 1:n) / (2 * σ^2)
    )
    optimize!(model)
    return
end

function perf_infopt_mle()
    Random.seed!(1234)
    n = 1_000
    data = randn(n)
    model = InfiniteOpt.InfiniteModel(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, μ, start = 0.0)
    @variable(model, σ >= 0.0, start = 1.0)
    @objective(
        model,
        Max,
        n / 2 * log(1 / (2 * π * σ^2)) -
        sum((data[i] - μ)^2 for i in 1:n) / (2 * σ^2)
    )
    optimize!(model)
    return
end

function perf_symbolics_mle()
    Random.seed!(1234)
    n = 1_000
    data = randn(n)
    Symbolics.@variables μ σ
    n / 2 * log(1 / (2 * π * σ^2)) -
    sum((data[i] - μ)^2 for i in 1:n) / (2 * σ^2)
    return
end

# clnlbeam

function perf_nl_clnlbeam()
    N = 1000
    h = 1 / N
    alpha = 350
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variables(model, begin
        -1 <= t[1:(N+1)] <= 1
        -0.05 <= x[1:(N+1)] <= 0.05
        u[1:(N+1)]
    end)
    @NLobjective(
        model,
        Min,
        sum(
            0.5 * h * (u[i+1]^2 + u[i]^2) +
            0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) for i in 1:N
        ),
    )
    @NLconstraint(
        model,
        [i = 1:N],
        x[i+1] - x[i] - 0.5 * h * (sin(t[i+1]) + sin(t[i])) == 0,
    )
    @constraint(
        model,
        [i = 1:N],
        t[i+1] - t[i] - 0.5 * h * u[i+1] - 0.5 * h * u[i] == 0,
    )
    optimize!(model)
    return
end

function perf_nlexpr_clnlbeam()
    N = 1000
    h = 1 / N
    alpha = 350
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variables(model, begin
        -1 <= t[1:(N+1)] <= 1
        -0.05 <= x[1:(N+1)] <= 0.05
        u[1:(N+1)]
    end)
    @objective(
        model,
        Min,
        sum(
            0.5 * h * (u[i+1]^2 + u[i]^2) +
            0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) for i in 1:N
        ),
    )
    @constraint(
        model,
        [i = 1:N],
        x[i+1] - x[i] - 0.5 * h * (sin(t[i+1]) + sin(t[i])) == 0,
    )
    @constraint(
        model,
        [i = 1:N],
        t[i+1] - t[i] - 0.5 * h * u[i+1] - 0.5 * h * u[i] == 0,
    )
    optimize!(model)
    return
end

function perf_infopt_clnlbeam()
    N = 1000
    h = 1 / N
    alpha = 350
    model = InfiniteOpt.InfiniteModel(Ipopt.Optimizer)
    set_silent(model)
    @variables(model, begin
        -1 <= t[1:(N+1)] <= 1
        -0.05 <= x[1:(N+1)] <= 0.05
        u[1:(N+1)]
    end)
    @objective(
        model,
        Min,
        sum(
            0.5 * h * (u[i+1]^2 + u[i]^2) +
            0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) for i in 1:N
        ),
    )
    @constraint(
        model,
        [i = 1:N],
        x[i+1] - x[i] - 0.5 * h * (sin(t[i+1]) + sin(t[i])) == 0,
    )
    @constraint(
        model,
        [i = 1:N],
        t[i+1] - t[i] - 0.5 * h * u[i+1] - 0.5 * h * u[i] == 0,
    )
    optimize!(model)
    return
end

# rosenbrock

function perf_nl_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    optimize!(model)
    return
end

function perf_nlexpr_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x)
    @variable(model, y)
    @objective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    optimize!(model)
    return
end

function perf_infopt_rosenbrock()
    model = InfiniteOpt.InfiniteModel(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x)
    @variable(model, y)
    @objective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    optimize!(model)
    return
end

function perf_symbolics_rosenbrock()
    Symbolics.@variables x y
    (1 - x)^2 + 100 * (y - x^2)^2
    return
end

# JuMP#2788

function perf_nl_jump_2788()
    N = 400
    Random.seed!(1234)
    k = N
    n = 12
    p = rand(400:700, k, 1)
    c1 = rand(100:200, k, n)
    c2 = 0.9 .* c1
    b = rand(150:250, k, 1)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, 0 <= x[i = 1:n] <= 1)
    @variable(model, 0 <= var1 <= 1)
    @variable(model, 0 <= var2 <= 1)
    @variable(model, 0 <= var3 <= 1)
    @objective(model, Max, var1 - var2 + var3)
    @NLexpression(model, expr, sum(x[i] * p[i] for i in 1:n))
    @NLexpression(model, expr_c1[j = 1:k], sum(x[i] * c1[j, i] for i in 1:n))
    @NLexpression(model, expr_c2[j = 1:k], sum(x[i] * c2[j, i] for i in 1:n))
    @NLconstraint(model, expr == sum(b[j] / (1 + var1)^j for j in 1:k))
    @NLconstraint(model, expr == sum(expr_c1[j] / (1 + var2)^j for j in 1:k))
    @NLconstraint(model, expr == sum(expr_c2[j] / (1 + var3)^j for j in 1:k))
    @NLconstraint(model, [j = 1:k], expr_c1[j] >= b[j])
    optimize!(model)
    return
end

function perf_nlexpr_jump_2788()
    N = 400
    Random.seed!(1234)
    k = N
    n = 12
    p = rand(400:700, k, 1)
    c1 = rand(100:200, k, n)
    c2 = 0.9 .* c1
    b = rand(150:250, k, 1)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, 0 <= x[i = 1:n] <= 1)
    @variable(model, 0 <= var1 <= 1)
    @variable(model, 0 <= var2 <= 1)
    @variable(model, 0 <= var3 <= 1)
    @objective(model, Max, var1 - var2 + var3)
    @expression(model, expr, sum(x[i] * p[i] for i in 1:n))
    @expression(model, expr_c1[j = 1:k], sum(x[i] * c1[j, i] for i in 1:n))
    @expression(model, expr_c2[j = 1:k], sum(x[i] * c2[j, i] for i in 1:n))
    @constraint(model, expr == sum(b[j] / (1 + var1)^Float64(j) for j in 1:k))
    @constraint(
        model,
        expr == sum(expr_c1[j] / (1 + var2)^Float64(j) for j in 1:k),
    )
    @constraint(
        model,
        expr == sum(expr_c2[j] / (1 + var3)^Float64(j) for j in 1:k),
    )
    @constraint(model, [j = 1:k], expr_c1[j] >= b[j])
    optimize!(model)
    return
end

function perf_infopt_jump_2788()
    N = 400
    Random.seed!(1234)
    k = N
    n = 12
    p = rand(400:700, k, 1)
    c1 = rand(100:200, k, n)
    c2 = 0.9 .* c1
    b = rand(150:250, k, 1)
    model = InfiniteOpt.InfiniteModel(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, 0 <= x[i = 1:n] <= 1)
    @variable(model, 0 <= var1 <= 1)
    @variable(model, 0 <= var2 <= 1)
    @variable(model, 0 <= var3 <= 1)
    @objective(model, Max, var1 - var2 + var3)
    @expression(model, expr, sum(x[i] * p[i] for i in 1:n))
    @expression(model, expr_c1[j = 1:k], sum(x[i] * c1[j, i] for i in 1:n))
    @expression(model, expr_c2[j = 1:k], sum(x[i] * c2[j, i] for i in 1:n))
    @constraint(model, expr == sum(b[j] / (1 + var1)^Float64(j) for j in 1:k))
    @constraint(
        model,
        expr == sum(expr_c1[j] / (1 + var2)^Float64(j) for j in 1:k),
    )
    @constraint(
        model,
        expr == sum(expr_c2[j] / (1 + var3)^Float64(j) for j in 1:k),
    )
    @constraint(model, [j = 1:k], expr_c1[j] >= b[j])
    optimize!(model)
    return
end

end  # module
