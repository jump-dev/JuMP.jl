# Based on https://matbesancon.github.io/post/2018-05-25-colgen2/
# and https://github.com/matbesancon/column_generation_jump

using JuMP, GLPK
using SparseArrays
const MOI = JuMP.MathOptInterface

maxwidth = 100
rollcost = 500
prices = Float64[167.0, 197.0, 281.0, 212.0, 225.0, 111.0, 93.0, 129.0, 108.0, 106.0, 55.0, 85.0, 66.0, 44.0, 47.0, 15.0, 24.0, 13.0, 16.0, 14.0]
widths = Float64[75.0, 75.0, 75.0, 75.0, 75.0, 53.8, 53.0, 51.0, 50.2, 32.2, 30.8, 29.8, 20.1, 16.2, 14.5, 11.0, 8.6, 8.2, 6.6, 5.1]
demand = Float64[38, 44, 30, 41, 36, 33, 36, 41, 35, 37, 44, 49, 37, 36, 42, 33, 47, 35, 49, 42]
nwidths = length(prices)

"""
    Solve the pricing problem to generate a new pattern, if there is at least one
    that could improve the current solutions.
"""
function solve_pricing(dual_demand_satisfaction, maxwidth, widths, rollcost, demand, prices)
    reduced_costs = dual_demand_satisfaction + prices
    n = length(reduced_costs)

    # The actual pricing model.
    submodel = Model(with_optimizer(GLPK.Optimizer))
    @variable(submodel, xs[1:n] >= 0, Int)
    @constraint(submodel, sum(xs .* widths) <= maxwidth)
    @objective(submodel, Max, sum(xs .* reduced_costs))

    optimize!(submodel)

    # If the net cost of this new pattern is nonnegative, no more patterns to add.
    new_pattern = round.(Int, value.(xs))
    net_cost = rollcost - sum(new_pattern .* (dual_demand_satisfaction .+ prices))

    println(net_cost)
    if net_cost >= 0 # No new pattern to add.
        return nothing
    else
        return new_pattern
    end
end

"""
    Run the full column-generation solving process, starting with an initial basis
    of cutting patterns and adding more and more of them by solving the pricing problem.
"""
function solve_cutting_stock(maxwidth, widths, rollcost, demand, prices; max_gen_cols=5000)
    n = length(widths)
    ncols = length(widths)

    # Initial set of patterns (stored in a sparse matrix: a pattern won't include many different cuts).
    patterns = spzeros(UInt16, n, ncols)
    for i in 1:n
        patterns[i, i] = min(floor(Int, maxwidth / widths[i]), round(Int, demand[i]))
    end

    # Write the master problem with this "reduced" set of patterns.
    # Not yet integer variables: otherwise, the dual values may make no sense
    # (actually, GLPK will yell at you if you're trying to get duals for integer problems)
    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, θ[1:ncols] >= 0)
    @objective(m, Min,
        sum(θ[p] * (rollcost - sum(patterns[j, p] * prices[j] for j in 1:n)) for p in 1:ncols)
    )
    @constraint(m, demand_satisfaction[j=1:n], sum(patterns[j, p] * θ[p] for p in 1:ncols) >= demand[j])

    # First solve of the master problem.
    optimize!(m)
    if termination_status(m) != MOI.OPTIMAL
        warn("Master not optimal ($ncols patterns so far)")
    end

    # Then, generate new patterns, based on the dual information.
    while ncols - n <= max_gen_cols # Generate at most max_gen_cols columns.
        if ! has_duals(m)
            break
        end

        new_pattern = try
            solve_pricing(dual.(demand_satisfaction), maxwidth, widths, rollcost, demand, prices)
        catch
            # At the final iteration, GLPK has dual values, but at least one of them is 0.0, and thus GLPK crashes.
            break
        end

        # No new pattern to add to the formulation: done!
        if new_pattern === nothing
            break
        end

        # Otherwise, add the new pattern to the master problem, recompute the duals,
        # and go on waltzing one more time with the pricing problem.
        ncols += 1
        patterns = hcat(patterns, new_pattern)

        # One new variable.
        new_var = @variable(m, [ncols], base_name="θ", lower_bound=0)
        push!(θ, new_var[ncols])

        # Update the objective function.
        set_objective_function(m, objective_function(m)
            + θ[ncols] * (rollcost - sum(patterns[j, ncols] * prices[j] for j=1:n))
        )

        # Update the constraint number j if the new pattern impacts this production.
        for j in 1:n
            if new_pattern[j] > 0
                set_standard_form_coefficient(demand_satisfaction[j], new_var[ncols], new_pattern[j])
            end
        end

        # Solve the new master problem to update the dual variables.
        optimize!(m)
        if termination_status(m) != MOI.OPTIMAL
            warn("Master not optimal ($ncols patterns so far)")
        end
    end

    # Finally, impose the master variables to be integer and resolve.
    # To be exact, at each node in the branch-and-bound tree, we would need to restart
    # the column generation process (just in case a new column would be interesting
    # to add). This way, we only get an upper bound (a feasible solution).
    set_integer.(θ)
    optimize!(m)

    if termination_status(m) != MOI.OPTIMAL
        warn("Final master not optimal ($ncols patterns)")
    end

    return value.(θ)
end

# Actually run the code.
θ_final = solve_cutting_stock(maxwidth, widths, rollcost, demand, prices)
