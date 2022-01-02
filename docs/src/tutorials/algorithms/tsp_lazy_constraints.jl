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

# # [Traveling Salesperson Problem (via callbacks)](@id tsp_lazy)
# 
# **Originally Contributed by**: Daniel Schermer
# 
# This tutorial describes how to implement the [Traveling Salesperson Problem](https://en.wikipedia.org/wiki/Travelling_salesman_problem) in JuMP using solver-independent lazy constraints that dynamically separate subtours.
# To be more precise, we use lazy constraints to cut off infeasible subtours only when necessary and not before needed.

using JuMP
import GLPK
import LinearAlgebra
import Random
import Plots

# # [Mathematical Formulation](@id tsp_model)
# 
# Assume that we are given a complete graph $\mathcal{G}(V,E)$ where $V$ is the set of vertices (or cities) and $E$ is the set of edges (or roads). 
# For each pair of vertices $i, j \in V, i \neq j$ the edge $(i,j) \in E$ is associated with a weight (or distance) $d_{ij} \in \mathbb{R}^+$.
# For this tutorial, we assume the problem to be symmetric, i.e., $d_{ij} = d_{ji} \, \forall i,j \in V$.

# In the Traveling Salesperson Problem, we are tasked with finding a tour with minimal length that visits every vertex exactly once and then returns to the point of origin, i.e., a *hamiltonian cycle* with minimal weight.
# 
# In order to model the problem, we introduce a binary variable $x_{ij} \in \{0,1\} \; \forall i, j \in V$ that indicates if edge $(i,j)$ is part of the tour or not.
# Using these variables, the Traveling Salesperson Problem can be modeled through the following Integer Linear Programming Formulation.
# 
# ## [Objective Function](@id tsp_objective)
# The objective consists in minimizing the weighted edges (due to the assumed symmetry, the second sum only contains $i<j$).
#
# ```math
# \text{min } \sum_{i \in V}  \sum_{j \in V, i < j} d_{ij} x_{ij}
# ```
#
# Note that it is also possible to use the following objective function instead.
#
# ```math
# \text{min } \sum_{i \in V}  \sum_{j \in V} \dfrac{d_{ij} x_{ij}}{2}
# ```
#
# ## [Constraints](@id tsp_constraints)

# There are four classes of constraints in our formulation.
# Due to the presumed symmetry, the following constraints must hold.
#
# ```math
# x_{ij} = x_{ji} \quad \forall i,j \in V
# ```
#
# For each vertex $i$, exactly two edges must be selected that connect it to other vertices $j$ in the graph $G$.
#
# ```math
# \sum_{j \in V} x_{ij} = 2 \quad \forall i \in V
# ```
# 
# We do not permit loops to occur.
#
# ```math
# x_{ii} = 0 \quad \forall i \in V
# ```
#
# The number of vertices can be adjusted here.
# The vertices are assumed to be randomly distributed in the Euclidean space;
# thus, the weight (distance) of each edge is defined as follows.

Random.seed!(1)
n = 75
X = 100 * rand(n)
Y = 100 * rand(n)
d = zeros(n, n)

for i in 1:n
    for j in 1:n
        d[i, j] = LinearAlgebra.norm([X[i] - X[j], Y[i] - Y[j]], 2)
    end
end

# For the JuMP model, we first initialize the model object.
# Then, we create the the binary decision variables and add the objective function and constraints (as introduced above).

# !!! warning
#     When using `GLPK` in conjunction with callbacks, the simple rounding heuristic (`sr_heur`) needs to be turned *off*.
#     Otherwise, the final solution will be infeasible, because the lazy constraint callback is not involved when a new incumbent is identified through the rounding heuristic.
#     For more details, refer to: [https://discourse.julialang.org/t/solution-foun-by-heuristic-glpk/38772](https://discourse.julialang.org/t/solution-foun-by-heuristic-glpk/38772).

tsp = Model(GLPK.Optimizer);
set_optimizer_attribute(tsp, "sr_heur", GLPK.GLP_OFF)

@variable(tsp, x[i in 1:n, j in 1:n], Bin);

@objective(tsp, Min, sum(x[i, j] * d[i, j] / 2 for i in 1:n for j in 1:n));

@constraint(tsp, symmetry[i in 1:n, j in i:n, i!=j], x[i, j] == x[j, i]);

@constraint(tsp, two_degree[i in 1:n], sum(x[i, j] for j in 1:n) == 2);

@constraint(tsp, no_loop[i in 1:n], x[i, i] == 0);

# ## [Subtour Elimination](@id tsp_sec)
# A major difficulty of the Traveling Salesperson Problem arises from the fact that we need to prevent *subtours*, i.e., several distinct Hamiltonian cycles existing on subgraphs of $G$.
# Note that the previous parts of the model (listed above) *do not* guarantee that the solution will be free of subtours.
# To this end, by $S$ we label a subset of vertices.
# Then, for each proper subset $S \subset V$, the following constraints guarantee that no subtour may occur.
#
# ```math
# \sum_{i \in S} \sum_{j \in S, i < j} x_{ij} \leq \vert S \vert - 1 \quad \forall S \subset V
# ```
# These constraints have the disadvantage that we would require exponentially many of them as $\vert V \vert$  increases. 
# Therefore, we will add these constraints not before needed and only when necessary.
#
# To search for violated constraints, based on the edges that are currently in the solution (i.e., those that have value $x_{ij} = 1$), we identify the shortest cycle through the function `subtour()`.
# Whenever a subtour has been identified, a constraint corresponding to the form above can be added to the model.

function subtour(edges)
    ## A list of all unvisited vertices
    unvisited = Set(collect(1:n))

    ## Placeholder for the shortest subtour
    cycle = collect(1:n)

    while !(isempty(unvisited))
        thiscycle = []
        neighbors = unvisited
        while !(isempty(neighbors))
            ## Get the first item         
            current = pop!(neighbors)

            ## Add it to the current cycle and remove it from unvisited
            push!(thiscycle, current)

            ## If we are in the first iteration of the inner while loop,
            ## then the previous pop! already removed `current'
            if length(thiscycle) > 1
                pop!(unvisited, current)
            end

            ## Get the index of all edges to which the current node is connected
            index = findall(edges -> edges[1] == current, edges)

            ## Based on the index, add the neighbors
            neighbors = []
            for i in index
                append!(neighbors, edges[i][2])
            end

            ## We only consider neighbors that have not yet been visited
            neighbors = intersect(neighbors, unvisited)
        end
        ## We always store the shortest cycle as subtour
        if length(thiscycle) < length(cycle)
            cycle = thiscycle
        end
    end
    return cycle
end

# Let us declare a helper function `selected_edges()` that will be repeatedly used in what follows.

function selected_edges(x)
    edges = Tuple{Int,Int}[]
    for i in 1:n
        for j in 1:n
            if (value(x[i, j]) > 0.5)
                push!(edges, (i, j))
            end
        end
    end
    return edges
end

# An iterative way of eliminating subtours is the following.
# However, it is reasonable to assume that this is not the most efficient way: Whenever a new subtour elimination constraint is added to the model, the optimization has to start from the very beginning.
# That way, the solver will repeatedly discard useful information encountered during previous solves (e.g., all cuts, the incumbent solution, or lower bounds).

# !!! info
#     Note that, in principle, it would also be feasible to add all subtours (instead of just the shortest one) to the model.
#     However, preventing just the shortest cycle is often sufficient for breaking other subtours and will keep the model size smaller.

optimize!(tsp)
time_iterated = solve_time(tsp);
cycle = subtour(selected_edges(x))

## Check for cycles, add the subtour elimination constraint,
## and optimize until no more cycles exist
while length(cycle) > 1 && length(cycle) < n
    S = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i < j]
    @constraint(tsp, sum(x[t[1], t[2]] for t in S) <= length(cycle) - 1)
    optimize!(tsp)
    global time_iterated += solve_time(tsp)
    global cycle = subtour(selected_edges(x))
end

info = (objective_value(tsp), time_iterated)
@show(info)

# A more sophisticated approach makes use of **lazy constraints**.
# To be more precise, we do this through the callback `subtour_elimination()` below, which is only run whenever we encounter a new integer-feasible solution.
# !!! warning
#     We use `seen` as a `Dict()` to store if we have seen a particular lazy constraint (determined by a cycle) before.
#     For performance reasons, when using `GLPK`, we must check if we have previously added a lazy constraint before, as `GLPK` does not check this before adding the constraint to the model.
#     If we were to ignore this; identical constraints might be added repeatedly, causing performance to drop drastically.

seen = Dict()

function subtour_elimination(cb_data)
    ## We only checkfor subtours when we encounter integer-feasible solutions
    status = callback_node_status(cb_data, tsp)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER

        ## Load the callback data at the current node  
        x_val = callback_value.(Ref(cb_data), x)

        ## Write the current edges in a tuple list and identify the shortest  cycle
        edges = selected_edges(x_val)
        cycle = subtour(edges)
        h = hash(sort(cycle))

        if !(h in keys(seen))
            ## A subtour contains at least 2 locations and at most (n-1)
            if length(cycle) > 1 && length(cycle) < n
                S = collect(
                    Iterators.filter(
                        t -> t[1] < t[2],
                        Iterators.product(cycle, cycle),
                    ),
                )
                con = @build_constraint(
                    sum(x[t[1], t[2]] for t in S) <= length(cycle) - 1
                )

                MOI.submit(tsp, MOI.LazyConstraint(cb_data), con)
                seen[h] = con
            end
        end
    end
end

# For solving the model with lazy constraints, we need to make sure that `LazyConstraints` are enabled and that the proper `CallbackFunction()` passed.
# Before we do this, we are going to quickly rebuild the model, such that previously added subtour elimination constraints are removed.

tsp = Model(GLPK.Optimizer);
set_optimizer_attribute(tsp, "sr_heur", GLPK.GLP_OFF)
## set_optimizer_attribute(tsp, "msg_lev", GLPK.GLP_MSG_ON)

@variable(tsp, x[i in 1:n, j in 1:n], Bin);

@objective(tsp, Min, sum(x[i, j] * d[i, j] / 2 for i in 1:n for j in 1:n));

@constraint(tsp, symmetry[i in 1:n, j in i:n, i!=j], x[i, j] == x[j, i]);

@constraint(tsp, two_degree[i in 1:n], sum(x[i, j] for j in 1:n) == 2);

@constraint(tsp, no_loop[i in 1:n], x[i, i] == 0);

MOI.set(tsp, MOI.LazyConstraintCallback(), subtour_elimination)
optimize!(tsp)
time_lazy = solve_time(tsp);

# !!! warning
#     The following is a bit of a hack, as different solvers appear to treat lazy constraints internally very differently!
#     While the following part is not necessary when using `Gurobi` as an Optimizer, in the case of `GLPK` things are different:
#     in combination with the dict `lazy_constraints` introduced above, some lazy constraints may be ignored, even though they were explicitely added to the model.
#     This might be remedied by removing the dict; however, then the optimization using lazy constraints will be very slow.
#     Therefore, for use with `GLPK`,  after the first solution we will verify that all previously identified lazy constraints are actually respected, by adding them to the model and re-solving it.
#     Note that this fix also works for the `sr_heur` issue discussed at the beginning of this document.

cycle = subtour(selected_edges(x))
while length(cycle) < n && length(cycle) > 1
    for k in keys(seen)
        add_constraint(tsp, seen[k])
    end
    #global lazy_constraints = Dict()
    optimize!(tsp)
    global time_lazy += solve_time(tsp)
    global cycle = subtour(selected_edges(x))
end

info = (objective_value(tsp), time_lazy)
@show(info)

# When we observe the output, we can see that the objective value for both approaches is identical.
# However, for the single instance of $n=75$ and the Optimizer `GLPK`, the solution time to find these solution is very different when using lazy constraints.

# As a quick sanity check, we might visualize the optimal tour to verify that no subtour is present.

plt = Plots.plot()
edges = selected_edges(x)
for e in edges
    Plots.plot!(
        [X[e[1]], X[e[2]]],
        [Y[e[1]], Y[e[2]]],
        legend = false,
        linecolor = :black,
    )
end
Plots.plot!()

# # References
# 
# The mathematical formulation was inspired by the following reference.
# 
# ```@raw html
# <a id='c1'></a>
# ```
# 
# 1. Gurobi Optimization, LLC. Gurobi Optimizer Reference Manual. (2021).
