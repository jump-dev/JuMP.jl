# Copyright (c) 2022 Daniel Schermer and contributors                            #src
#                                                                                #src
# Permission is hereby granted, free of charge, to any person obtaining a copy   #src
# of this software and associated documentation files (the "Software"), to deal  #src
# in the Software without restriction, including without limitation the rights   #src
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      #src
# copies of the Software, and to permit persons to whom the Software is          #src
# furnished to do so, subject to the following conditions:                       #src
#                                                                                #src
# The above copyright notice and this permission notice shall be included in all #src
# copies or substantial portions of the Software.                                #src
#                                                                                #src
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #src
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       #src
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    #src
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         #src
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  #src
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  #src
# SOFTWARE.

# # [Traveling Salesperson Problem](@id tsp_lazy)

# **Originally Contributed by**: Daniel Schermer

# This tutorial describes how to implement the
# [Traveling Salesperson Problem](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
# in JuMP using solver-independent lazy constraints that dynamically separate
# subtours. To be more precise, we use lazy constraints to cut off infeasible
# subtours only when necessary and not before needed.

# It uses the following packages:

using JuMP
import GLPK
import Random
import Plots

# ## [Mathematical Formulation](@id tsp_model)

# Assume that we are given a complete graph $\mathcal{G}(V,E)$ where $V$ is the
# set of vertices (or cities) and $E$ is the set of edges (or roads). For each
# pair of vertices $i, j \in V, i \neq j$ the edge $(i,j) \in E$ is associated
# with a weight (or distance) $d_{ij} \in \mathbb{R}^+$.

# For this tutorial, we assume the problem to be symmetric, that is,
# $d_{ij} = d_{ji} \, \forall i,j \in V$.

# In the Traveling Salesperson Problem, we are tasked with finding a tour with
# minimal length that visits every vertex exactly once and then returns to the
# point of origin, that is, a *hamiltonian cycle* with minimal weight.
#
# To model the problem, we introduce a binary variable,
# $x_{ij} \in \{0,1\} \; \forall i, j \in V$, that indicates if edge $(i,j)$ is
# part of the tour or not. Using these variables, the Traveling Salesperson
# Problem can be modeled as the following integer linear program.

# ### [Objective Function](@id tsp_objective)

# The objective is to minimize the length of the tour (due to the assumed
# symmetry, the second sum only contains $i<j$):
# ```math
# \text{min\ } \sum_{i \in V}  \sum_{j \in V, i < j} d_{ij} x_{ij}.
# ```

# Note that it is also possible to use the following objective function instead:
# ```math
# \text{min } \sum_{i \in V}  \sum_{j \in V} \dfrac{d_{ij} x_{ij}}{2}.
# ```

# ### [Constraints](@id tsp_constraints)

# There are four classes of constraints in our formulation.

# First, due to the presumed symmetry, the following constraints must hold:
# ```math
# x_{ij} = x_{ji} \quad \forall i,j \in V.
# ```

# Second, for each vertex $i$, exactly two edges must be selected that connect it
# to other vertices $j$ in the graph $G$:
# ```math
# \sum_{j \in V} x_{ij} = 2 \quad \forall i \in V.
# ```

# Third, we do not permit loops to occur:
# ```math
# x_{ii} = 0 \quad \forall i \in V.
# ```

# The fourth constraint is more complicated. A major difficulty of the Traveling
# Salesperson Problem arises from the fact that we need to prevent *subtours*,
# that is, several distinct Hamiltonian cycles existing on subgraphs of $G$.

# Note that the previous constraints *do not* guarantee that the solution will
# be free of subtours. To this end, by $S$ we label a subset of vertices.
# Then, for each proper subset $S \subset V$, the following constraints
# guarantee that no subtour may occur:
# ```math
# \sum_{i \in S} \sum_{j \in S, i < j} x_{ij} \leq \vert S \vert - 1 \quad \forall S \subset V.
# ```

# Problematically, we require exponentially many of these constraints as
# $\vert V \vert$ increases. Therefore, we will add these constraints only when
# necessary.

# ## [Implementation](@id tsp_implementation)

# There are two ways we can eliminate subtours in JuMP, both of which will be
# shown in what follows:
# - iteratively solving a new model that incorporates previously identified
#   subtours,
# - or adding violated subtours as *lazy constraints*.

# ### Data

# The vertices are assumed to be randomly distributed in the Euclidean space;
# thus, the weight (distance) of each edge is defined as follows.

function generate_distance_matrix(n; random_seed = 1)
    Random.seed!(random_seed)
    X = 100 * rand(n)
    Y = 100 * rand(n)
    d = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1:n, j in 1:n]
    return X, Y, d
end

n = 40
X, Y, d = generate_distance_matrix(n)

# For the JuMP model, we first initialize the model object. Then, we create the
# binary decision variables and add the objective function and constraints. By
# defining the `x` matrix as `Symmetric`, we do not need to add explicit
# constraints that `x[i, j] == x[j, i]`.

function build_tsp_model(d, n)
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:n, 1:n], Bin, Symmetric)
    @objective(model, Min, sum(d .* x) / 2)
    @constraint(model, [i in 1:n], sum(x[i, :]) == 2)
    @constraint(model, [i in 1:n], x[i, i] == 0)
    return model
end

# To search for violated constraints, based on the edges that are currently in
# the solution (that is, those that have value $x_{ij} = 1$), we identify the
# shortest cycle through the function `subtour()`. Whenever a subtour has been
# identified, a constraint corresponding to the form above can be added to the
# model.

function subtour(edges::Vector{Tuple{Int,Int}}, n)
    shortest_subtour, unvisited = collect(1:n), Set(collect(1:n))
    while !isempty(unvisited)
        this_cycle, neighbors = Int[], unvisited
        while !isempty(neighbors)
            current = pop!(neighbors)
            push!(this_cycle, current)
            if length(this_cycle) > 1
                pop!(unvisited, current)
            end
            neighbors =
                [j for (i, j) in edges if i == current && j in unvisited]
        end
        if length(this_cycle) < length(shortest_subtour)
            shortest_subtour = this_cycle
        end
    end
    return shortest_subtour
end

# Let us declare a helper function `selected_edges()` that will be repeatedly
# used in what follows.

function selected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] > 0.5]
end

# Other helper functions for computing subtours:

subtour(x::Matrix{Float64}) = subtour(selected_edges(x, size(x, 1)), size(x, 1))
subtour(x::AbstractMatrix{VariableRef}) = subtour(value.(x))

# ### Iterative method

# An iterative way of eliminating subtours is the following.

# However, it is reasonable to assume that this is not the most efficient way:
# Whenever a new subtour elimination constraint is added to the model, the
# optimization has to start from the very beginning.

# That way, the solver will repeatedly discard useful information encountered
# during previous solves (e.g., all cuts, the incumbent solution, or lower
# bounds).

# !!! info
#     Note that, in principle, it would also be feasible to add all subtours
#     (instead of just the shortest one) to the model. However, preventing just
#     the shortest cycle is often sufficient for breaking other subtours and
#     will keep the model size smaller.

iterative_model = build_tsp_model(d, n)
optimize!(iterative_model)
time_iterated = solve_time(iterative_model)
cycle = subtour(iterative_model[:x])
while 1 < length(cycle) < n
    println("Found cycle of length $(length(cycle))")
    S = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i < j]
    @constraint(
        iterative_model,
        sum(iterative_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1,
    )
    optimize!(iterative_model)
    global time_iterated += solve_time(iterative_model)
    global cycle = subtour(iterative_model[:x])
end

objective_value(iterative_model)

#-

time_iterated

# As a quick sanity check, we visualize the optimal tour to verify that no
# subtour is present:

function plot_tour(X, Y, x)
    plt = Plots.plot()
    for (i, j) in selected_edges(x, size(x, 1))
        Plots.plot!([X[i], X[j]], [Y[i], Y[j]], legend = false)
    end
    return plt
end

plot_tour(X, Y, value.(iterative_model[:x]))

# ### Lazy constraint method

# A more sophisticated approach makes use of _lazy constraints_. To be more
# precise, we do this through the `subtour_elimination_callback()` below, which
# is only run whenever we encounter a new integer-feasible solution.

lazy_model = build_tsp_model(d, n)
function subtour_elimination_callback(cb_data)
    status = callback_node_status(cb_data, lazy_model)
    if status != MOI.CALLBACK_NODE_STATUS_INTEGER
        return  # Only run at integer solutions
    end
    cycle = subtour(callback_value.(cb_data, lazy_model[:x]))
    if !(1 < length(cycle) < n)
        return  # Only add a constraint if there is a cycle
    end
    println("Found cycle of length $(length(cycle))")
    S = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i < j]
    con = @build_constraint(
        sum(lazy_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1,
    )
    MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), con)
    return
end
MOI.set(lazy_model, MOI.LazyConstraintCallback(), subtour_elimination_callback)
optimize!(lazy_model)
objective_value(lazy_model)

# This finds the same optimal tour:

plot_tour(X, Y, value.(lazy_model[:x]))

# Surprisingly, for this particular model with GLPK, the solution time is worse
# than the iterative method:

time_lazy = solve_time(lazy_model)

# In most other cases and solvers, however, the lazy time should be faster than
# the iterative method.
