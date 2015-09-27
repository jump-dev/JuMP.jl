using JuMP

function bench(n)
    t1 = @elapsed begin
        m = Model()
        @defVar(m, x[1:n,2:n])
    end
    t2 = @elapsed begin
        cntr = 0
        for (ii,jj,v) in x
            cntr += ii + jj + v.col
        end
    end
    t1, t2
end

bench(1)
bench(100)
bench(1000)
bench(2000)

function bench_runtime(n)
    t1 = @elapsed begin
        m = Model()
        I, J = 1:n, 2:n
        @defVar(m, x[I,J])
    end
    t2 = @elapsed begin
        cntr = 0
        for (ii,jj,v) in x
            cntr += ii + jj + v.col
        end
    end
    t1, t2
end

bench_runtime(1)
bench_runtime(100)
bench_runtime(1000)
bench_runtime(2000)
