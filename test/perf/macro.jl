# macro.jl
# Macro-exercising speed tests
using JuMP

function test_linear(N)
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

function test_quad(N)
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

# Warmup
println("Test 1")
test_linear(1)
test_quad(1)
for N in [20, 50, 100]
    println("  Running N=$(N)...")
    N1_times = Any[]
    N2_times = Any[]
    for iter in 1:10
        push!(N1_times, @elapsed test_linear(N))
        push!(N2_times, @elapsed test_quad(N))
    end
    println("    N=$(N) min $(minimum(N1_times))")
    println("    N=$(N) min $(minimum(N2_times))")
end
