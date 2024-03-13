#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module JuMPBenchmarks

using JuMP

import BenchmarkTools
import LinearAlgebra
import Random

function benchmark(;
    baseline::String = "baseline",
    compare_against::Bool = false,
    directory::String = @__DIR__,
    report_filename::String = joinpath(@__DIR__, "report.txt"),
    kwargs...,
)
    group = BenchmarkTools.BenchmarkGroup()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "bench_")
            f = getfield(@__MODULE__, name)
            group[name] = BenchmarkTools.@benchmarkable $f()
        end
    end
    if compare_against
        MOI.Benchmarks.compare_against_baseline(
            group,
            baseline;
            directory = directory,
            report_filename = report_filename,
            kwargs...,
        )
    else
        MOI.Benchmarks.create_baseline(
            group,
            baseline;
            directory = directory,
            kwargs...,
        )
    end
    return
end

"""
    run_microbenchmark(N::Int)

Run each micro benchmark function `N` times.

 * If `N == 0`, uses the `BenchmarkTools.@benchmark` macro.
 * If `N > 0`, runs the function multiple times with `@time`.

!!! warning
    This is not a rigorous benchmark if `N > 0`. It uses `@time`, so it doesn't
    accurately capture top-level compilation and inference. Use `N = 0` for a
    more rigorous evaluation, or use `JuMPBenchmarks.benchmark` for a comparison
    with prior benchmarks.
"""
function run_microbenchmark(N::Int)
    for name in sort!(names(@__MODULE__; all = true))
        if startswith("$(name)", "bench_")
            f = getfield(@__MODULE__, name)
            @info("$(name)")
            if N == 0
                result = BenchmarkTools.@benchmark $f()
                display(result)
                println()
            else
                for _ in 1:N
                    @time f()
                end
            end
        end
    end
    return
end

"""
    run_examples(julia_cmd = "julia")

Run each example in the `/docs/src/tutorials` folder in a fresh Julia instance.
"""
function run_examples(julia_cmd = "julia")
    src_dir = joinpath(dirname(dirname(@__DIR__)), "docs", "src", "tutorials")
    doc_project = joinpath(dirname(dirname(@__DIR__)), "docs")
    log_file = joinpath(@__DIR__, "run_examples.log")
    for (root, dir, files) in walkdir(src_dir)
        for file in files
            if !endswith(file, ".jl")
                continue
            end
            filename = joinpath(root, file)
            code = """
            data_1 = @timed include(\"$(filename)\")
            data_2 = @timed include(\"$(filename)\")
            open("$(log_file)", "a") do io
                println(io, "$(file)")
                println(io, "  Time (s): ", data_1.time)
                println(io, "  Bytes   : ", data_1.bytes)
                println(io, "  -- run 2 --")
                println(io, "  Time (s): ", data_2.time)
                println(io, "  Bytes   : ", data_2.bytes)
            end
            """
            println("Running: ", file)
            @time redirect_stderr(devnull) do
                return redirect_stdout(devnull) do
                    return run(
                        `$(julia_cmd) --project=$(doc_project) -e $(code)`,
                    )
                end
            end
        end
    end
    return
end

###
### Matrix Product
###

function _matrix_affine_product(n::Int)
    model = Model()
    a = rand(n, n)
    @variable(model, x[1:n, 1:n])
    return @expression(model, a * x)
end

function _matrix_quadratic_product(n::Int)
    model = Model()
    a = rand(n, n)
    @variable(model, x[1:n, 1:n])
    return @expression(model, a * x * a)
end

for s in (:affine, :quadratic), N in [10, 50]
    f = getfield(@__MODULE__, Symbol("_matrix_$(s)_product"))
    new_name = Symbol("bench_matrix_$(s)_product_$(lpad(N, 3, '0'))")
    @eval $(new_name)() = $f($N)
end

###
### Vector speed test
###

function _vector_sum(n::Int)
    model = Model()
    c = rand(n)
    @variable(model, z[1:n])
    @constraint(model, sum(c[i] * z[i] for i in 1:n) <= 1)
end

function _vector_dot(n::Int)
    model = Model()
    c = rand(n)
    @variable(model, z[1:n])
    @constraint(model, LinearAlgebra.dot(c, z) <= 1)
end

function _matrix_sum(n::Int)
    model = Model()
    a = rand(n, n)
    @variable(model, x[1:n, 1:n])
    @constraint(model, sum(a[i, j] * x[i, j] for i in 1:n, j in 1:n) <= 1)
end

function _matrix_dot(n::Int)
    model = Model()
    a = rand(n, n)
    @variable(model, x[1:n, 1:n])
    @constraint(model, LinearAlgebra.dot(a, x) <= 1)
end

function _array_sum(n::Int)
    model = Model()
    b = rand(n, n, n)
    @variable(model, y[1:n, 1:n, 1:n])
    @constraint(
        model,
        sum(b[i, j, k] * y[i, j, k] for i in 1:n, j in 1:n, k in 1:n) <= 1,
    )
end

function _array_dot(n::Int)
    model = Model()
    b = rand(n, n, n)
    @variable(model, y[1:n, 1:n, 1:n])
    @constraint(model, LinearAlgebra.dot(b, y) <= 1)
end

for s in (:vector, :matrix, :array), N in [10, 50]
    f = getfield(@__MODULE__, Symbol("_$(s)_sum"))
    new_name = Symbol("bench_$(s)_sum_$(lpad(N, 3, '0'))")
    @eval $(new_name)() = $f($N)
end

###
### Macro
###

function _macro_linear(N::Int)
    m = Model()
    @variable(m, x[1:10N, 1:5N])
    @variable(m, y[1:N, 1:N, 1:N])
    for z in 1:10
        @constraint(
            m,
            9 * y[1, 1, 1] - 5 * y[N, N, N] -
            2 * sum(z * x[j, i*N] for j in ((z-1)*N+1):z*N, i in 3:4) +
            sum(i * (9 * x[i, j] + 3 * x[j, i]) for i in N:2N, j in N:2N) +
            x[1, 1] +
            x[10N, 5N] +
            x[2N, 1] +
            1 * y[1, 1, N] +
            2 * y[1, N, 1] +
            3 * y[N, 1, 1] +
            y[N, N, N] - 2 * y[N, N, N] + 3 * y[N, N, N] <=
            sum(
                sum(
                    sum(
                        N * i * j * k * y[i, j, k] + x[i, j] for
                        k in 1:N if i != j && j != k
                    ) for j in 1:N
                ) for i in 1:N
            ) + sum(
                sum(x[i, j] for j in 1:5N if j % i == 3) for
                i in 1:10N if i <= N * z
            )
        )
    end
end

function _macro_quad(N::Int)
    m = Model()
    @variable(m, x[1:10N, 1:5N])
    @variable(m, y[1:N, 1:N, 1:N])

    for z in 1:10
        @constraint(
            m,
            9 * y[1, 1, 1] - 5 * y[N, N, N] -
            2 * sum(z * x[j, i*N] for j in ((z-1)*N+1):z*N, i in 3:4) +
            sum(i * (9 * x[i, j] + 3 * x[j, i]) for i in N:2N, j in N:2N) +
            x[1, 1] +
            x[10N, 5N] * x[2N, 1] +
            1 * y[1, 1, N] * 2 * y[1, N, 1] +
            3 * y[N, 1, 1] +
            y[N, N, N] - 2 * y[N, N, N] * 3 * y[N, N, N] <=
            sum(
                sum(
                    sum(
                        N * i * j * k * y[i, j, k] * x[i, j] for
                        k in 1:N if i != j && j != k
                    ) for j in 1:N
                ) for i in 1:N
            ) + sum(
                sum(x[i, j] for j in 1:5N if j % i == 3) for
                i in 1:10N if i <= N * z
            )
        )
    end
end

for s in (:linear, :quad), N in [10, 50]
    f = getfield(@__MODULE__, Symbol("_macro_$(s)"))
    new_name = Symbol("bench_macro_$(s)_$(lpad(N, 3, '0'))")
    @eval $(new_name)() = $f($N)
end

###
### Printing
###

function benchmark_print_AffExpr()
    m = Model()
    N = 100
    @variable(m, x[1:N])
    c = @constraint(m, sum(i * x[i] for i in 1:N) >= N)
    return sprint(print, c)
end

function benchmark_print_model()
    m = Model()
    N = 100
    @variable(m, x[1:N])
    @constraint(m, sum(i * x[i] for i in 1:N) >= N)
    return sprint(print, m)
end

function benchmark_print_model_10000()
    m = Model()
    N = 10_000
    @variable(m, x[1:N])
    @constraint(m, sum(i * x[i] for i in 1:N) >= N)
    return sprint(print, m)
end

function benchmark_print_small_model()
    m = Model()
    N = 10
    @variable(m, x1[1:N])
    @variable(m, x2[1:N, f = 1:N])
    @variable(m, x3[1:N, f = 1:2:N])
    @variable(m, x4[[:a, :b, :c]])
    @variable(m, x5[[:a, :b, :c], [:d, "e", 4]])
    @constraint(
        m,
        sum(i * x1[i] for i in 1:N) +
        sum(i * f * x2[i, f] for i in 1:N, f in 1:N) +
        sum(i * f * x3[i, f] for i in 1:N, f in 1:2:N) +
        sum(x4) >= N
    )
    return sprint(print, m)
end

###
### Axis constraints
###

function _sum_iterate(con_refs)
    x = 0.0
    for con_ref in con_refs
        x += dual(con_ref)
    end
    return x
end

function _sum_index(con_refs)
    x = 0.0
    for i in eachindex(con_refs)
        x += dual(con_refs[i])
    end
    return x
end

function _dense_axis_constraints(key = :index)
    n = 1_000
    model = Model()
    mock = MOIU.MockOptimizer(
        MOIU.Model{Float64}();
        eval_variable_constraint_dual = false,
    )
    MOIU.set_mock_optimize!(
        mock,
        mock -> MOIU.mock_optimize!(
            mock,
            zeros(n),
            (MOI.VariableIndex, MOI.EqualTo{Float64}) => ones(n - 1),
        ),
    )
    MOIU.reset_optimizer(model, mock)
    @variable(model, x[1:n])
    set = MOI.EqualTo(0.0)
    con_refs = @constraint(model, [i = 2:n], x[i] in set)
    optimize!(model)
    if key == :index
        _sum_index(con_refs)
    else
        _sum_iterate(con_refs)
    end
    return
end

function _sparse_axis_constraints(key = :index)
    n = 1_000
    model = Model()
    mock = MOIU.MockOptimizer(
        MOIU.Model{Float64}();
        eval_variable_constraint_dual = false,
    )
    MOIU.set_mock_optimize!(
        mock,
        mock -> MOIU.mock_optimize!(
            mock,
            zeros(n),
            (MOI.VariableIndex, MOI.EqualTo{Float64}) => ones(div(n, 2)),
        ),
    )
    MOIU.reset_optimizer(model, mock)
    @variable(model, x[1:n])
    set = MOI.EqualTo(0.0)
    con_refs = @constraint(model, [i = 1:n; iseven(i)], x[i] in set)
    optimize!(model)
    if key == :index
        _sum_index(con_refs)
    else
        _sum_iterate(con_refs)
    end
    return
end

for container in (:dense, :sparse), sum_type in (:iterate, :index)
    f = getfield(@__MODULE__, Symbol("_$(container)_axis_constraints"))
    new_name = Symbol("bench_$(container)_axis_constraints_$(sum_type)")
    @eval $(new_name)() = $f($(sum_type))
end

###
### Model constructors
###

"""
    benchmark_p_median(
        num_facilities = 100,
        num_customers = 100,
        num_locations = 5_000,
    )

Implements the "p-median" facility location problem. We try to locate N
facilities such that we minimize the distance any given customer has to travel
to reach their closest facility. In this simple instance we will operate
in a 1D world with L possible locations for facilities, and customers being
located at random locations along the number line from 1 to D.

We use anonymous variables to remove the cost of name generation from the
benchmark.
"""
function benchmark_p_median(
    num_facilities = 100,
    num_customers = 100,
    num_locations = 5_000,
)
    Random.seed!(10)
    customer_locations = [rand(1:num_locations) for _ in 1:num_customers]
    model = Model()
    has_facility = @variable(model, [1:num_locations], Bin)
    is_closest = @variable(model, [1:num_locations, 1:num_customers], Bin)
    @objective(
        model,
        Min,
        sum(
            abs(customer_locations[customer] - location) *
            is_closest[location, customer] for customer in 1:num_customers,
            location in 1:num_locations
        )
    )
    for customer in 1:num_customers
        # `location` can't be closest for `customer` if there is no facility.
        @constraint(
            model,
            [location in 1:num_locations],
            is_closest[location, customer] <= has_facility[location]
        )
        # One facility must be the closest for `customer`.
        @constraint(
            model,
            sum(
                is_closest[location, customer] for location in 1:num_locations
            ) == 1
        )
    end
    # Must place all facilities.
    @constraint(model, sum(has_facility) == num_facilities)
    return model
end

"""
    benchmark_cont5(n = 500)

Based on a linear-Quadratic control problem (cont5_2_1) from one of Hans
Mittleman's instance collections.

We use anonymous variables to remove the cost of name generation from the
benchmark.
"""
function benchmark_cont5(n = 500)
    m = n
    n1 = n - 1
    m1 = m - 1
    dx = 1 / n
    T = 1.58
    dt = T / m
    h2 = dx^2
    a = 0.001
    yt = [0.5 * (1 - (j * dx)^2) for j in 0:n]
    model = Model()
    y = @variable(model, [0:m, 0:n], lower_bound = 0, upper_bound = 1)
    u = @variable(model, [1:m], lower_bound = -1, upper_bound = 1)
    @objective(model, Min, y[0, 0])
    # PDE
    for i in 0:m1
        for j in 1:n1
            @constraint(
                model,
                h2 * (y[i+1, j] - y[i, j]) ==
                0.5 *
                dt *
                (
                    y[i, j-1] - 2 * y[i, j] + y[i, j+1] + y[i+1, j-1] -
                    2 * y[i+1, j] + y[i+1, j+1]
                )
            )
        end
    end
    # Initial conditions.
    for j in 0:n
        @constraint(model, y[0, j] == 0)
    end
    # Boundary conditions.
    for i in 1:m
        @constraint(model, y[i, 2] - 4 * y[i, 1] + 3 * y[i, 0] == 0)
        @constraint(
            model,
            y[i, n-2] - 4 * y[i, n1] + 3 * y[i, n] ==
            (2 * dx) * (u[i] - y[i, n])
        )
    end
    return model
end

end  # module

function _print_help()
    return println(
        """
julia test/perf/JuMPBenchmarks.jl [-r N] [-f file [--compare]] [-j julia]

Run a script to benchmark various aspects of JuMP.

## Run a simple micro-benchmark

Pass `-r N` to run each benchmark function `N` times.

### Example

```
\$ julia test/perf/JuMPBenchmarks.jl -r 2
```

## Run a rigorous micro-benchmark

To run a more rigorous benchmark, do not pass `-r` and pass `-f file`
instead.

 * If `--compare` is not given, save a new benchmark dataset to `file`.
 * If `--compare`, compare aginst the data in `file`.

### Example

```
\$ julia test/perf/JuMPBenchmarks.jl -f my_benchmark_run
# Make changes to JuMP, then run:
\$ julia test/perf/JuMPBenchmarks.jl -f my_benchmark_run --compare
```

## Run the examples

To run all of the examples, pass `-j julia` where `julia` is the command to
start Julia at the command line. This uses the Project.toml located at
`/docs/src/Project.toml`, which assumes you have dev'd JuMP to it as
appropriate.

### Example

```
\$ julia test/perf/JuMPBenchmarks.jl -j /Users/oscar/julia1.6
```
""",
    )
end

function _run(f, key)
    arg = findfirst(isequal(key), ARGS)
    if arg !== nothing
        f(ARGS[arg+1])
    end
    return
end

if length(ARGS) > 0
    if findfirst(isequal("-h"), ARGS) !== nothing
        _print_help()
        exit(0)
    end
    _run("-r") do N
        return JuMPBenchmarks.run_microbenchmark(parse(Int, N))
    end
    _run("-j") do julia
        return JuMPBenchmarks.run_examples(julia)
    end
    _run("-f") do baseline
        compare_against = findfirst(isequal("--compare"), ARGS) !== nothing
        return JuMPBenchmarks.benchmark(;
            baseline = baseline,
            compare_against = compare_against,
        )
    end
else
    _print_help()
end
