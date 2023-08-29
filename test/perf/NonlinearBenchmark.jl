module NonlinearBenchmark

using JuMP
import BenchmarkTools
import DataFrames
import Ipopt
import Random

function benchmark_group()
    lookup = Dict("perf_nl_" => "@NL", "perf_nlexpr_" => "NonlinearExpr")
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
    results = BenchmarkTools.run(suite)
    df_time = build_table(x -> minimum(x).time / 1e9, results)
    df_memory = build_table(x -> minimum(x).memory / 1024^2, results)
    @info "minimum(time) [s]"
    display(df_time)
    @info "minimum(memory) [MiB]"
    display(df_memory)
    return results
end

function build_table(f, results)
    tables = map(sort!(collect(keys(results["NonlinearExpr"])))) do b
        old = f(results["@NL"][b])
        new = f(results["NonlinearExpr"][b])
        return (benchmark = b, NL = old, NonlinearExpr = new, ratio = new / old)
    end
    return DataFrames.DataFrame(tables)
end

# sum
#
# nlexpr is slower because it builds up the product via operator overloading,
# creating a lot of temporary objects. @NL gets to see the full +(args...) to it
# builds the expression in-place.
#
# We could fix this by implementing a n-argy method for +, but that gets
# difficult with method ambiguities.

function perf_nl_micro_sum()
    model = Model()
    @variable(model, x)
    @NLobjective(model, Min, sum(x^i for i in 1:10_000))
    return
end

function perf_nlexpr_micro_sum()
    model = Model()
    @variable(model, x)
    @objective(model, Min, sum(x^i for i in 1:10_000))
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

function perf_nl_micro_prod()
    model = Model()
    @variable(model, x)
    @NLobjective(model, Min, prod(x^i for i in 1:10_000))
    return
end

function perf_nlexpr_micro_prod()
    model = Model()
    @variable(model, x)
    @objective(model, Min, prod(x^i for i in 1:10_000))
    return
end

# many_constraints

function perf_nl_micro_many_constraints()
    model = Model()
    @variable(model, x[1:10_000])
    @NLconstraint(model, [i = 1:10_000], sin(x[i]) <= cos(i))
    return
end

function perf_nlexpr_micro_many_constraints()
    model = Model()
    @variable(model, x[1:10_000])
    @constraint(model, [i = 1:10_000], sin(x[i]) <= cos(i))
    return
end

# value_expr_many_small

function perf_nl_micro_value_expr_many_small()
    model = Model()
    @variable(model, x)
    @NLexpression(model, expr[i = 1:10_000], x^i)
    value.(x -> 2.0, expr)
    return
end

function perf_nlexpr_micro_value_expr_many_small()
    model = Model()
    @variable(model, x)
    @expression(model, expr[i = 1:10_000], x^i)
    value.(x -> 2.0, expr)
    return
end

# value_expr_few_large

function perf_nl_micro_value_expr_few_large()
    model = Model()
    @variable(model, x)
    @NLexpression(model, expr, sum(x^i for i in 1:10_000))
    value(x -> 2.0, expr)
    return
end

function perf_nlexpr_micro_value_expr_few_large()
    model = Model()
    @variable(model, x)
    @expression(model, expr, sum(x^i for i in 1:10_000))
    value(x -> 2.0, expr)
    return
end

# mle

function perf_nl_model_mle()
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

function perf_nlexpr_model_mle()
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

# clnlbeam

function perf_nl_model_clnlbeam()
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

function perf_nlexpr_model_clnlbeam()
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

# rosenbrock

function perf_nl_model_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    optimize!(model)
    return
end

function perf_nlexpr_model_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x)
    @variable(model, y)
    @objective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    optimize!(model)
    return
end

# JuMP#2788

function perf_nl_model_jump_2788()
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

function perf_nlexpr_model_jump_2788()
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
    @constraint(model, expr == sum(b[j] / (1 + var1)^j for j in 1:k))
    @constraint(model, expr == sum(expr_c1[j] / (1 + var2)^j for j in 1:k),)
    @constraint(model, expr == sum(expr_c2[j] / (1 + var3)^j for j in 1:k),)
    @constraint(model, [j = 1:k], expr_c1[j] >= b[j])
    optimize!(model)
    return
end

# nested_problems

function perf_nl_model_nested_problems()
    function solve_lower_level(x...)
        model = Model(Ipopt.Optimizer)
        set_silent(model)
        @variable(model, y[1:2])
        @NLobjective(
            model,
            Max,
            x[1]^2 * y[1] + x[2]^2 * y[2] - x[1] * y[1]^4 - 2 * x[2] * y[2]^4,
        )
        @constraint(model, (y[1] - 10)^2 + (y[2] - 10)^2 <= 25)
        optimize!(model)
        @assert termination_status(model) == LOCALLY_SOLVED
        return objective_value(model), value.(y)
    end
    function V(x...)
        f, _ = solve_lower_level(x...)
        return f
    end
    function ∇V(g::AbstractVector, x...)
        _, y = solve_lower_level(x...)
        g[1] = 2 * x[1] * y[1] - y[1]^4
        g[2] = 2 * x[2] * y[2] - 2 * y[2]^4
        return
    end
    function ∇²V(H::AbstractMatrix, x...)
        _, y = solve_lower_level(x...)
        H[1, 1] = 2 * y[1]
        H[2, 2] = 2 * y[2]
        return
    end
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x[1:2] >= 0)
    register(model, :f_V, 2, V, ∇V, ∇²V)
    @NLobjective(model, Min, x[1]^2 + x[2]^2 + f_V(x[1], x[2]))
    optimize!(model)
    solution_summary(model)
    return
end

function perf_nlexpr_model_nested_problems()
    function solve_lower_level(x...)
        model = Model(Ipopt.Optimizer)
        set_silent(model)
        @variable(model, y[1:2])
        @objective(
            model,
            Max,
            x[1]^2 * y[1] + x[2]^2 * y[2] - x[1] * y[1]^4 - 2 * x[2] * y[2]^4,
        )
        @constraint(model, (y[1] - 10)^2 + (y[2] - 10)^2 <= 25)
        optimize!(model)
        @assert termination_status(model) == LOCALLY_SOLVED
        return objective_value(model), value.(y)
    end
    function V(x...)
        f, _ = solve_lower_level(x...)
        return f
    end
    function ∇V(g::AbstractVector, x...)
        _, y = solve_lower_level(x...)
        g[1] = 2 * x[1] * y[1] - y[1]^4
        g[2] = 2 * x[2] * y[2] - 2 * y[2]^4
        return
    end
    function ∇²V(H::AbstractMatrix, x...)
        _, y = solve_lower_level(x...)
        H[1, 1] = 2 * y[1]
        H[2, 2] = 2 * y[2]
        return
    end
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x[1:2] >= 0)
    @operator(model, f_V, 2, V, ∇V, ∇²V)
    @objective(model, Min, x[1]^2 + x[2]^2 + f_V(x[1], x[2]))
    optimize!(model)
    solution_summary(model)
    return
end

# ###

function perf_nl_model_votroto()
    Q = -0.8:0.4:0.8
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, -2 <= p[1:5] <= 2)
    @variable(model, -1 <= w <= 3)
    @variable(model, -1 <= q <= 3)
    @objective(model, Min, w)
    total = Dict(
        _q => @NLexpression(
            model,
            sum(
                _p / sqrt(2π) * exp(-(i - _q)^2 / 2) for
                (i, _p) in enumerate(p)
            )
        ) for _q in Any[Q; q; 0.5]
    )
    l1 = Dict(
        _q => @NLexpression(model, 1 - total[_q] + 0.5 * total[0.5]) for
        _q in Any[Q; q]
    )
    @NLconstraint(
        model,
        [_q in Q],
        w * (l1[q] - l1[_q]) + (1 - w) * (total[q] - 1) <= 0
    )
    optimize!(model)
    return
end

function perf_nlexpr_model_votroto()
    Q = -0.8:0.4:0.8
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, -2 <= p[1:5] <= 2)
    @variable(model, -1 <= w <= 3)
    @variable(model, -1 <= q <= 3)
    @objective(model, Min, w)
    f(p, q) = (1 / sqrt(2π)) * exp(-((p - q)^2) / 2)
    total(p, q) = sum(_p * f(i, q) for (i, _p) in enumerate(p))
    l1(p, q) = 1 - total(p, q) + 0.5 * total(p, 0.5)
    l2(p, q) = total(p, q) - 1
    lhs(p, q, _q) = l1(p, q) - l1(p, _q)
    @constraint(model, [_q in Q], w * lhs(p, q, _q) + (1 - w) * l2(p, q) <= 0)
    optimize!(model)
    return
end

# large_expressions

function perf_nl_model_large_expressions()
    N = 50_000
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_attribute(model, "max_iter", 1)
    @variable(model, y[1:N], start = 1)
    @variable(model, z[1:N])
    @NLobjective(model, Max, sum(2z[i]^2 + sin(1 / y[i]) for i in 1:N))
    @NLconstraint(
        model,
        [i = 1:N],
        ifelse(z[i] <= y[i]^3, log(y[i] / i), z[i] / cos(y[i])) <= 42,
    )
    @NLconstraint(model, sum(z[i]^i + log(y[i]) for i in 1:N) == 0)
    optimize!(model)
    return
end

function perf_nlexpr_model_large_expressions()
    N = 50_000
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_attribute(model, "max_iter", 1)
    @variable(model, y[1:N], start = 1)
    @variable(model, z[1:N])
    @objective(model, Max, sum(2z[i]^2 + sin(1 / y[i]) for i in 1:N))
    @constraint(
        model,
        [i = 1:N],
        ifelse(z[i] <= y[i]^3, log(y[i] / i), z[i] / cos(y[i])) <= 42,
    )
    @constraint(model, sum(z[i]^i + log(y[i]) for i in 1:N) == 0)
    optimize!(model)
    return
end

# large_expressions_2

function perf_nl_model_large_expressions_2()
    N = 100
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_attribute(model, "max_iter", 1)
    @variable(model, y[1:N], start = 1)
    @variable(model, z[1:N])
    @NLobjective(model, Max, sum(2z[i]^2 + sin(1 / y[i]) for i in 1:N))
    @NLconstraint(
        model,
        prod(
            ifelse(z[i] <= y[i]^3, log(y[i] / i), z[i] / cos(y[i])) for i in 1:N
        ) <= 42
    )
    @NLconstraint(model, sum(z[i]^i + log(y[i]) for i in 1:N) == 0)
    optimize!(model)
    return
end

function perf_nlexpr_model_large_expressions_2()
    N = 100
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_attribute(model, "max_iter", 1)
    @variable(model, y[1:N], start = 1)
    @variable(model, z[1:N])
    @objective(model, Max, sum(2z[i]^2 + sin(1 / y[i]) for i in 1:N))
    @constraint(
        model,
        prod(
            ifelse(z[i] <= y[i]^3, log(y[i] / i), z[i] / cos(y[i])) for i in 1:N
        ) <= 42
    )
    @constraint(model, sum(z[i]^i + log(y[i]) for i in 1:N) == 0)
    optimize!(model)
    return
end

end  # module
