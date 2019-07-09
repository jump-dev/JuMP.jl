# Based on https://matbesancon.github.io/post/2018-05-25-colgen2/
# and https://github.com/matbesancon/column_generation_jump

using JuMP, Cbc
const MOI = JuMP.MathOptInterface

maxwidth = 100
cost = 500
prices = Float64[167.0, 197.0, 281.0, 212.0, 225.0, 111.0, 93.0, 129.0, 108.0, 106.0, 55.0, 85.0, 66.0, 44.0, 47.0, 15.0, 24.0, 13.0, 16.0, 14.0]
widths = Float64[75.0, 75.0, 75.0, 75.0, 75.0, 53.8, 53.0, 51.0, 50.2, 32.2, 30.8, 29.8, 20.1, 16.2, 14.5, 11.0, 8.6, 8.2, 6.6, 5.1]
demand = Float64[38, 44, 30, 41, 36, 33, 36, 41, 35, 37, 44, 49, 37, 36, 42, 33, 47, 35, 49, 42]
nwidths = length(prices)

"""
    Solve the subproblem: find a new cutting pattern.
"""
function subproblem(reduced_costs, sizes, maxcapacity)
    n = length(reduced_costs)

    subm = Model(with_optimizer(Cbc.Optimizer))
    xs = @variable(subm, xs[1:n] >= 0, Int)
    @constraint(subm, sum(xs.*sizes) <= maxcapacity)
    @objective(subm, Max, sum(xs .* reduced_costs))

    optimize!(subm)
    return round.(Int, value.(xs)), round(Int, objective_value(subm))
end

"""
    Initialize the master problem.
"""
function init_master(maxwidth, widths, rollcost, demand, prices)
    n = length(widths)
    ncols = length(widths)
    patterns = spzeros(UInt16, n, ncols)
    for i in 1:n
        patterns[i, i] = min(floor(Int, maxwidth / widths[i]), round(Int, demand[i]))
    end

    m = Model(with_optimizer(Cbc.Optimizer))
    θ = @variable(m, θ[1:ncols] >= 0)
    @objective(m, Min,
        sum(θ[p] * (rollcost - sum(patterns[j, p] * prices[j] for j=1:n)) for p in 1:ncols)
    )
    @constraint(m, demand_satisfaction[j=1:n], sum(patterns[j, p] * θ[p] for p in 1:ncols) >= demand[j])

    optimize!(m)
    if termination_status(m) != MOI.OPTIMAL
        warn("No optimal")
    end
    return m, value.(θ), demand_satisfaction, patterns
end

function column_generation(maxwidth, widths, rollcost, demand, prices; maxcols = 5000)
    (m, θ, demand_satisfaction, patterns) = init_master(maxwidth, widths, rollcost, demand, prices)
    ncols = nwidths
    while ncols <= maxcols
        reduced_costs = getdual(demand_satisfaction) + prices
        newcol, newobj = subproblem(reduced_costs, widths, maxwidth)
        netcost = cost - sum(newcol[j]*(getdual(demand_satisfaction)[j]+prices[j]) for j in 1:nwidths)
        println("New reduced cost: $netcost")
        if netcost >= 0
            return (MOI.OPTIMAL, patterns, getvalue(θ))
        end
        patterns = hcat(patterns, newcol)
        ncols += 1

        m = Model(with_optimizer(Cbc.Optimizer))
        θ = @variable(m, θ[1:ncols] >= 0)
        @objective(m, Min,
            sum(θ[p] * (rollcost - sum(patterns[j ,p] * prices[j] for j=1:nwidths)) for p in 1:ncols)
        )
        @constraint(m, demand_satisfaction[j=1:nwidths], sum(patterns[j, p] * θ[p] for p in 1:ncols) >= demand[j])

        optimize!(m)
        if termination_status(m) != MOI.OPTIMAL
            warn("No optimal")
            return termination_status(m), patterns, value.(θ)
        end
    end
    return (:NotFound, patterns, :NoVariable)
end

"""
    From patterns built in the column generation phase, find an integer solution
"""
function branched_model(patterns, demand, rollcost, prices; npatts = size(patterns)[2])
    npatts = size(patterns)[2]

    m = Model(solver = CbcSolver())
    θ = @variable(m, θ[p = 1:npatts] >= 0, Int)
    @objective(m, Min,
        sum(θ[p] * (rollcost - sum(patterns[j, p] * prices[j] for j=1:nwidths)) for p in 1:npatts)
    )
    @constraint(m, demand_satisfaction[j=1:nwidths], sum(θ[p] * patterns[j, p] for p in 1:npatts) >= demand[j])

    optimize!(m)
    return termination_status(m), round.(Int, value.(θ))
end

# Actually run the code.
status, θ_final = branched_model(patterns, demand, cost, prices)
