# This benchmarks measures the overhead in accessing `moibackend`.
# See https://github.com/JuliaOpt/JuMP.jl/pull/1348

using BenchmarkTools
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
using JuMP

function modif(m::Model, cref)
    MOI.modifyconstraint!(m.moibackend, cref.index, MOI.EqualTo(0.0))
end

function bench(m::Model, cref)
    @btime JuMP.num_variables($m)
    @btime modif($m, $cref)
end

optimizer = MOIU.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
#m = Model(backend=optimizer, mode=JuMP.Direct) # syntax before #1348
m = Model(optimizer) # syntax after #1348
@variable m x
@variable m y
cref = @constraint m x + y == 1
bench(m, cref)
