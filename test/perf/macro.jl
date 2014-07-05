# macro.jl
# Macro-exercising speed tests
using JuMP

function test_1(N)
    m = Model()
    @defVar(m, x[1:10N,1:5N])
    @defVar(m, y[1:N,1:N,1:N])

    @addConstraint(m, cons[z=1:10],
        y[1,1,1]*9 - 5*y[N,N,N] -
        2*sum{ z*x[j,i*N],                j=((z-1)*N+1):z*N, i=3:4} +
          sum{ i*(9*x[i,j] + x[j,i]*3),   i=N:2N,            j=N:2N} + 
        x[1,1] + x[10N,5N] + x[2N,1] + 
        y[1,1,N]*1 + 2*y[1,N,1] + 3*y[N,1,1] +
        y[N,N,N] - 2*y[N,N,N] + 3*y[N,N,N] 
         <=
        sum{sum{sum{N*i*j*k*y[i,j,k] + x[i,j],k=1:N; i!=j && j!=k},j=1:N},i=1:N} +
        sum{sum{x[i,j], j=1:5N; j % i == 3}, i=1:10N; i <= z*N}
        )
end

# Warmup
println("Test 1")
test_1(1)
for N in [20,50,100]
    println("  Running N=$(N)...")
    N_times = {}
    for iter in 1:10
        tic()
        test_1(N)
        push!(N_times, toq())
    end
    println("    N=$(N) min $(minimum(N_times))")
end

