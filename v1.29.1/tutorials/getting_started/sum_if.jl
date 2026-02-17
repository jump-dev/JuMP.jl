# Copyright (c) 2024 Oscar Dowson and contributors                               #src
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
# SOFTWARE.                                                                      #src

# # Performance problems with sum-if formulations

# The purpose of this tutorial is to explain a common performance issue that can
# arise with summations like `sum(x[a] for a in list if condition(a))`. This
# issue is particularly common in models with graph or network structures.

# !!! tip
#     This tutorial is more advanced than the other "Getting started" tutorials.
#     It's in the "Getting started" section because it is one of the most common
#     causes of performance problems that users experience when they first start
#     using JuMP to write large scale programs. If you are new to JuMP, you may
#     want to briefly skim the tutorial and come back to it once you have
#     written a few JuMP models.

# ## Required packages

# This tutorial uses the following packages

using JuMP
import Plots

# ## Data

# As a motivating example, we consider a network flow problem, like the examples
# in [Network flow problems](@ref) or [The network multi-commodity flow problem](@ref).

# Here is a function that builds a random graph. The specifics do not matter.

function build_random_graph(num_nodes::Int, num_edges::Int)
    nodes = 1:num_nodes
    edges = Pair{Int,Int}[i - 1 => i for i in 2:num_nodes]
    while length(edges) < num_edges
        edge = rand(nodes) => rand(nodes)
        if !(edge in edges)
            push!(edges, edge)
        end
    end
    function demand(n)
        if n == 1
            return -1
        elseif n == num_nodes
            return 1
        else
            return 0
        end
    end
    return nodes, edges, demand
end

nodes, edges, demand = build_random_graph(4, 8)

# The goal is to decide the flow of a commodity along each edge in `edges` to
# satisfy the `demand(n)` of each node `n` in `nodes`.

# The mathematical formulation is:
#
# ```math
# \begin{aligned}
# s.t. && \sum_{(i,n)\in E} x_{i,n} - \sum_{(n,j)\in E} x_{n,j} = d_n && \forall n \in N\\
# && x_{e} \ge 0 && \forall e \in E
# \end{aligned}
# ```

# ## Naïve model

# The first model you might write down is:

model = Model()
@variable(model, flows[e in edges] >= 0)
@constraint(
    model,
    [n in nodes],
    sum(flows[(i, j)] for (i, j) in edges if j == n) -
    sum(flows[(i, j)] for (i, j) in edges if i == n) == demand(n)
);

# The benefit of this formulation is that it looks very similar to the
# mathematical formulation of a network flow problem.

# The downside to this formulation is subtle. Behind the scenes, the JuMP
# `@constraint` macro expands to something like:

model = Model()
@variable(model, flows[e in edges] >= 0)
for n in nodes
    flow_in = AffExpr(0.0)
    for (i, j) in edges
        if j == n
            add_to_expression!(flow_in, flows[(i, j)])
        end
    end
    flow_out = AffExpr(0.0)
    for (i, j) in edges
        if i == n
            add_to_expression!(flow_out, flows[(i, j)])
        end
    end
    @constraint(model, flow_in - flow_out == demand(n))
end

# This formulation includes two for-loops, with a loop over every edge (twice) for
# every node. The [big-O notation](https://en.wikipedia.org/wiki/Big_O_notation)
# of the runtime is ``O(|nodes| \times |edges|)``. If
# you have a large number of nodes and a large number of edges, the runtime of
# this loop can be large.

# Let's build a function to benchmark our formulation:

function build_naive_model(nodes, edges, demand)
    model = Model()
    @variable(model, flows[e in edges] >= 0)
    @constraint(
        model,
        [n in nodes],
        sum(flows[(i, j)] for (i, j) in edges if j == n) -
        sum(flows[(i, j)] for (i, j) in edges if i == n) == demand(n)
    )
    return model
end

nodes, edges, demand = build_random_graph(1_000, 2_000)
@elapsed build_naive_model(nodes, edges, demand)

# A good way to benchmark is to measure the runtime across a wide range of input
# sizes. From our big-O analysis, we should expect that doubling the number of
# nodes and edges results in a 4x increase in the runtime.

run_times = Float64[]
factors = 1:10
for factor in factors
    GC.gc()  #src
    graph = build_random_graph(1_000 * factor, 5_000 * factor)
    push!(run_times, @elapsed build_naive_model(graph...))
end
Plots.plot(; xlabel = "Factor", ylabel = "Runtime [s]")
Plots.scatter!(factors, run_times; label = "Actual")
a, b = hcat(ones(10), factors .^ 2) \ run_times
Plots.plot!(factors, a .+ b * factors .^ 2; label = "Quadratic fit")

# As expected, the runtimes demonstrate quadratic scaling: if we double the
# number of nodes and edges, the runtime increases by a factor of four.

# ## Caching

# We can improve our formulation by caching the list of incoming and outgoing
# nodes for each node `n`:

out_nodes = Dict(n => Int[] for n in nodes)
in_nodes = Dict(n => Int[] for n in nodes)
for (i, j) in edges
    push!(out_nodes[i], j)
    push!(in_nodes[j], i)
end

# with the corresponding change to our model:

model = Model()
@variable(model, flows[e in edges] >= 0)
@constraint(
    model,
    [n in nodes],
    sum(flows[(i, n)] for i in in_nodes[n]) -
    sum(flows[(n, j)] for j in out_nodes[n]) == demand(n)
);

# The benefit of this formulation is that we now loop over `out_nodes[n]`
# rather than `edges` for each node `n`, and so the runtime is ``O(|edges|)``.

# Let's build a new function to benchmark our formulation:

function build_cached_model(nodes, edges, demand)
    out_nodes = Dict(n => Int[] for n in nodes)
    in_nodes = Dict(n => Int[] for n in nodes)
    for (i, j) in edges
        push!(out_nodes[i], j)
        push!(in_nodes[j], i)
    end
    model = Model()
    @variable(model, flows[e in edges] >= 0)
    @constraint(
        model,
        [n in nodes],
        sum(flows[(i, n)] for i in in_nodes[n]) -
        sum(flows[(n, j)] for j in out_nodes[n]) == demand(n)
    )
    return model
end

nodes, edges, demand = build_random_graph(1_000, 2_000)
@elapsed build_cached_model(nodes, edges, demand)

# ## Analysis

# Now we can analyse the difference in runtime of the two formulations:

run_times_naive = Float64[]
run_times_cached = Float64[]
factors = 1:10
for factor in factors
    GC.gc()  #src
    graph = build_random_graph(1_000 * factor, 5_000 * factor)
    push!(run_times_naive, @elapsed build_naive_model(graph...))
    push!(run_times_cached, @elapsed build_cached_model(graph...))
end
Plots.plot(; xlabel = "Factor", ylabel = "Runtime [s]")
Plots.scatter!(factors, run_times_naive; label = "Actual")
a, b = hcat(ones(10), factors .^ 2) \ run_times_naive
Plots.plot!(factors, a .+ b * factors .^ 2; label = "Quadratic fit")
Plots.scatter!(factors, run_times_cached; label = "Cached")
a, b = hcat(ones(10), factors) \ run_times_cached
Plots.plot!(factors, a .+ b * factors; label = "Linear fit")

# Even though the cached model needs to build `in_nodes` and `out_nodes`, it is
# asymptotically faster than the naïve model, scaling linearly with `factor`
# rather than quadratically.

# ## Lesson

# If you write code with `sum-if` type conditions, for example,
# `@constraint(model, [a in set], sum(x[b] for b in list if condition(a, b))`,
# you can improve the performance by caching the elements for which `condition(a, b)`
# is true.

# Finally, you should understand that this behavior is not specific to JuMP, and
# that it applies more generally to all computer programs you might write.
# (Python programs that use Pyomo or gurobipy would similarly benefit from this
# caching approach.)
#
# Understanding big-O notation and algorithmic complexity is a useful debugging
# skill to have, regardless of the type of program that you are writing.
