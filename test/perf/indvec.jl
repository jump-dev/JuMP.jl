using JuMP

function bench_add(n)
    v = JuMP.IndexedVector(Float64,n)
    for i = 1:2:n
        JuMP.addelt!(v, i, 1.0)
    end
    for i = 2:2:n
        JuMP.addelt!(v, i, 1.0)
    end
    for i = 1:3:n
        JuMP.addelt!(v, i, 1.0)
    end
end

bench_add(10)
gc_disable()
@time bench_add(1000000)
@time bench_add(1000000)
@time bench_add(1000000)
gc_enable()

