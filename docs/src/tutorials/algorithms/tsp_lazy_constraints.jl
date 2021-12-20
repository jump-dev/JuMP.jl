
# Copyright (c) 2019 Arpit Bhatia and contributors                               #src
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


# ## [Traveling Salesperson Problem (via callbacks)](@id tsp_lazy)
# 
# **Originally Contributed by**: Daniel Schermer
# 
# This notebook describes how to implement the Traveling Salesperson Problem in JuMP using lazy constraints that dynamically separate subtours.
# To be more precise, we use lazy constraints to cut off infeasible subtours only when necessary and not before needed.
# The model has been tested with Julia Version 1.7.0, JuMP.jl Version 0.22.1, and GLPK Version 0.15.2.


using JuMP
import Combinatorics
import GLPK
import LinearAlgebra
import Random


# # Mathematical Formulation
# 
# Assume that we are given a complete graph $\mathcal{G}(V,E)$ where $V$ is the set of vertices and $E$ is the set of edges. 
# For each pair of vertices $i, j \in V, i \neq j$ the edge $(i,j) \in E$ is associated with a distance (or weight) $d_{ij} \in \mathbb{R}^+$.
# For the purpose of this Notebook, we assume the problem to be symmetric, i.e., $d_{ij} = d_{ji} \, \forall i,j \in V$.
# In the Traveling Salesperson Problem, we are tasked with finding a tour with minimal length that visits every vertex exactly once and then returns to the point of origin, i.e., a hamiltonian cycle with minimal weight.
# 
# In order to model the problem, we introduce a binary variable $x_{ij} \in \{0,1\} \; \forall i, j \in V$ that indicates if edge $(i,j)$ is part of the tour or not.
# Using these variables, the Traveling Salesperson Problem can be modeled through the following Integer Linear Programming Formulation:
# 
# ## Objective Function
# The objective consists in minimizing the weighted edges (due to the assumed symmetry, the second sum only contains $j>i$).
# $$\text{min } \sum_{i \in V}  \sum_{j \in V, j > i} d_{ij} x_{ij}$$
# 
# ## Constraints
# For the formulation of our problem, we need four classes of constraints which we will discuss in what follows.
# 
# ### Symmetry
# Due to the presumed symmetry, the following constraints must hold.
# 
# $$ x_{ij} = x_{ji} \quad \forall (i,j) \in V $$
# 
# ### Degree Constraints
# For each vertex $i$, exactly two edges must be selected that connect it to other vertices $j$.
# $$ \sum_{i,j \in V} x_{ij} = 2 \quad \forall i \in V $$
# 
# ### Loops
# We do not permit loops.
# $$ x_{ii} = 0 \quad \forall i \in V $$


Random.seed!(1)

# Number of vertices
n = 25

# The vertices are assumed to be randomly distributed in the Euclidean space
X = 100 * rand(n)
Y = 100 * rand(n)

# Thus, the distance (weight) of each edge is defined as follows
d = zeros(n, n)
for i = 1:n
    for j = 1:n
        d[i, j] = LinearAlgebra.norm([X[i] - X[j], Y[i] - Y[j]], 2)
    end
end

# Initialize the model object
tsp = Model(GLPK.Optimizer)

# Initialize the binary decision variables
@variable(tsp, x[i in 1:n, j in 1:n], Bin)

@constraint(tsp, symmetry[i in 1:n, j in 1:n], x[i, j] == x[j, i]);

@constraint(tsp, two_degree[i in 1:n], sum(x[i, j] for j = 1:n) == 2);

@constraint(tsp, no_loop[i in 1:n], x[i, i] == 0);

@objective(tsp, Min, sum(x[i, j] * d[i, j] for i = 1:n for j = 1:n))


# ## Subtour Elimination
# A major difficulty of the Traveling Salesperson Problem arises from the fact that we need to prevent *subtours*, i.e., several distinct Hamiltonian cycles existing on distinct subgraphs of $G$.
# Note that the previous parts of the model (listed above) *do not* guarantee that the solution will be free of subtours.
# 
# To this end, by $S$ we label a subset of vertices. Then, for each proper subset $S \subset V$, the following constraints guarantee that no subtour may occur.
# 
# $$ \sum_{i \in S} \sum_{j \in S, j \neq i} x_{ij} \leq \vert S \vert - 1 \quad \forall S \subset V $$
# 
# In general, we would require an exponential number of subtour eliminations constraints as $\vert V \vert$ increases.
# Therefore, in this model, we will add these constraints as **lazy constraints** , i.e., not before needed and only when necessary.
# 
# We do this through the callback **subtour_elimination()** below, which is only run whenever we encounter an integer-feasible solution.
# Based on the edges that are currently in the solution, we identify the shortest cycle through the function **subtour()**.
# Once the subtour has been identified, a corresponding subtour elimination constraint is added to the model.


function subtour_elimination(cb_data)
    # We only checkfor subtours when we encounter integer-feasible solutions
    status = callback_node_status(cb_data, tsp)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER

        # Load the callback data at the current node  
        x_val = zeros(n, n)
        for i = 1:n
            for j = 1:n
                x_val[i, j] = callback_value(cb_data, x[i, j])
            end
        end

        # Write the current edges in a tuple list
        edges = Tuple{Int,Int}[]
        for i = 1:n
            for j = 1:n
                if (x_val[i, j] > 0.5)
                    push!(edges, (i, j))
                end
            end
        end

        # Get the shortest cycle from the list of edges
        cycle = subtour(edges)

        # A subtour contains at least 2 locations and at most (n-1)
        if length(cycle) > 1 && length(cycle) < n
            subtour_edges = subtour_edges_helper(cycle)
            con = @build_constraint(sum(x[e[1], e[2]] for e in subtour_edges) <= length(cycle) - 1)
            MOI.submit(tsp, MOI.LazyConstraint(cb_data), con)
        end
    end
end

# Helper for constraint building
function subtour_edges_helper(cycle)
    subtour_edges = []
    for i = 1:length(cycle)-1
        push!(subtour_edges, [cycle[i], cycle[i+1]])
    end
    push!(subtour_edges, [cycle[length(cycle)], cycle[1]])
    return subtour_edges
end


# Given a list of edges, the function **subtour()**, is programmed to identify the shortest subtour as follows.
# Note that the resulting cycle is passed back to **subtour_elimination()**.


function subtour(edges)
    # A list of all unvisited vertices
    unvisited = Set(collect(1:n))

    # Placeholder for the shortest subtour
    cycle = collect(1:n)

    while !(isempty(unvisited))
        thiscycle = []
        neighbors = unvisited
        while !(isempty(neighbors))
            # Get the first item         
            current = pop!(neighbors)

            # Add it to the current cycle and remove it from unvisited
            push!(thiscycle, current)

            # If we are in the first iteration of the inner while loop,
            # then the previous pop! already removed `current'
            if length(thiscycle) > 1
                pop!(unvisited, current)
            end

            # Get the index of all edges to which the current node is connected
            index = findall(edges -> edges[1] == current, edges)

            # Based on the index, add the neighbors
            neighbors = []
            for i in index
                append!(neighbors, edges[i][2])
            end

            # We only consider neighbors that have not yet been visited
            neighbors = intersect(neighbors, unvisited)
        end
        # We always store the shortest cycle as subtour
        if length(thiscycle) < length(cycle)
            cycle = thiscycle
        end
    end
    return cycle
end


# All that is left to do is to run the model.
# To do this, we need to make sure that **LazyConstraints** are enabled and that the proper **Callback** is passed.


# Lazy constraints need to be set; otherwise, subtours might exist.
MOI.set(tsp, MOI.LazyConstraintCallback(), subtour_elimination)

# CAREFUL: When using GPLK in conjunction with callbacks, the simple rounding heuristic needs to be turned off.
# Otherwise, the final solution will be infeasible, because the lazy constraint callback is not involved.
# For more details, refer to: https://discourse.julialang.org/t/solution-foun-by-heuristic-glpk/38772/7
set_optimizer_attribute(tsp, "sr_heur", GLPK.GLP_OFF)
set_optimizer_attribute(tsp, "msg_lev", GLPK.GLP_MSG_ON)

optimize!(tsp)


# As a quick sanity check, we might visualize the optimal tour to verify that no subtour is present.


import Plots

plt = plot()
for i = 1:n
    for j = i:n
        if (value.(x[i, j]) > 0.8)
            plot!([X[i], X[j]], [Y[i], Y[j]], legend = false, linecolor = :black)
        end
    end
end
plot!()


# # References
# 
# The mathematical formulation was inspired by the following reference.
# 
# ```@raw html
# <a id='c1'></a>
# ```
# 
# 1. Gurobi Optimization, LLC. Gurobi Optimizer Reference Manual. (2021).