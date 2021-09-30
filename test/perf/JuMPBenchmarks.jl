module JuMPBenchmarks

using JuMP
using LinearAlgebra
using BenchmarkTools

function benchmark(;
    baseline::String = "baseline",
    compare_against::Bool = false,
    directory::String = @__DIR__,
    report_filename::String = joinpath(@__DIR__, "report.txt"),
    kwargs...,
)
    group = BenchmarkTools.BenchmarkGroup()
    for name in names(@__MODULE__, all = true)
        if startswith("$(name)", "bench_")
            f = getfield(@__MODULE__, name)
            group[name] = @benchmarkable $f()
        end
    end
    if compare_against
        MOI.Benchmarks.compare_against_baseline(
            group,
            baseline;
            directory = directory,
            report_filename = report_filename,
            kwargs...
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
    run(N::Int)

Run each benchmark function `N` times.

!!! warning
    This is not a rigorous benchmark. It uses `@time`, so it doesn't accurately
    capture top-level compilation and inference. Use `benchmark` for a more
    rigorous comparison.
"""
function run(N::Int)
    for name in sort!(names(@__MODULE__, all = true))
        if startswith("$(name)", "bench_")
            f = getfield(@__MODULE__, name)
            @info("$(name)")
            for _ in 1:N
                @time f()
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
                        k = 1:N if i != j && j != k
                    ) for j in 1:N
                ) for i in 1:N
            ) + sum(
                sum(x[i, j] for j = 1:5N if j % i == 3) for
                i = 1:10N if i <= N * z
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
                        k = 1:N if i != j && j != k
                    ) for j in 1:N
                ) for i in 1:N
            ) + sum(
                sum(x[i, j] for j = 1:5N if j % i == 3) for
                i = 1:10N if i <= N * z
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
    c = @constraint(m, sum(i*x[i] for i=1:N) >= N)
    return sprint(print, c)
end

function benchmark_print_model()
    m = Model()
    N = 100
    @variable(m, x[1:N])
    @constraint(m, sum(i*x[i] for i=1:N) >= N)
    return sprint(print, m)
end

function benchmark_print_model_10000()
    m = Model()
    N = 10_000
    @variable(m, x[1:N])
    @constraint(m, sum(i*x[i] for i=1:N) >= N)
    return sprint(print, m)
end

function benchmark_print_small_model()
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
    return sprint(print, m)
end

end  # module

function _print_help()
    println("""
    julia test/perf/JuMPBenchmarks.jl [-r N] [-f file] [--compare]

    Run a script to benchmark JuMP.

    ## Run a simple benchmark

     * Pass `-r N` to run each benchmark function `N` times.

    ## Run a rigorous benchmark

    To run a more rigorous benchmark, do not pass `-r` and pass `-f file`
    instead.

     * If `--compare` is not given, save a new benchmark dataset to `file`.
     * If `--compare`, compare aginst the data in `file`.

    ## Example

    Run a simple benchmark
    ```
    \$ julia test/perf/JuMPBenchmarks.jl -r 2
    ```

    Run a complicated benchmark:
    ```
    \$ julia test/perf/JuMPBenchmarks.jl -f my_benchmark_run
    # Make changes to JuMP, then run:
    \$ julia test/perf/JuMPBenchmarks.jl -f my_benchmark_run --compare
    ```
    """)
end

if length(ARGS) > 0
    file_index = findfirst(isequal("-f"), ARGS)
    run = findfirst(isequal("-r"), ARGS)
    if run !== nothing
        JuMPBenchmarks.run(parse(Int, ARGS[run+1]))
        exit(0)
    end
    if file_index === nothing
        _print_help()
    end
    JuMPBenchmarks.benchmark(
        baseline = ARGS[file_index+1],
        compare_against = findfirst(isequal("--compare"), ARGS) !== nothing,
    )
else
    _print_help()
end
