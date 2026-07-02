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
# SOFTWARE.                                                                      #src

# # Network Flows

# **Originally Contributed by**: Arpit Bhatia

# In graph theory, a flow network (also known as a transportation network) is a
# directed graph where each edge has a capacity and each edge receives a flow.
# The amount of flow on an edge cannot exceed the capacity of the edge.

# Often in operations research, a directed graph is called a network, the
# vertices are called nodes and the edges are called arcs.

# A flow must satisfy the restriction that the amount of flow into a node equals
# the amount of flow out of it,  unless it is a source, which has only outgoing
# flow, or sink, which has only incoming flow.

# A network can be used to model traffic in a computer network, circulation with
# demands, fluids in pipes,  currents in an electrical circuit, or anything
# similar in which something travels through a network of nodes.

using JuMP
import GLPK
import LinearAlgebra

# ## The Shortest Path Problem

# Suppose that each arc $(i, j)$ of a graph is assigned a scalar cost $a_{i,j}$,
# and suppose that we define the cost of a forward path to be the sum of the
# costs of its arcs.

# Given a pair of nodes, the shortest path problem is to find a forward path
# that connects these nodes and has minimum cost.

# ```math
# \begin{aligned}
# \min && \sum_{\forall e(i,j) \in E} a_{i,j} \times x_{i,j} \\
# s.t. && b(i) = \sum_j x_{ij} - \sum_k x_{ki} = \begin{cases}
# 1 &\mbox{if $i$ is the starting node,} \\
# -1 &\mbox{if $i$ is the ending node,} \\
# 0 &\mbox{otherwise.} \end{cases} \\
# && x_{e} \in \{0,1\} && \forall e \in E
# \end{aligned}
# ```

# ![Flow Network 1](../../assets/g1.svg)

G = [
    0 100 30 0 0
    0 0 20 0 0
    0 0 0 10 60
    0 15 0 0 50
    0 0 0 0 0
]

n = size(G)[1]

shortest_path = Model(GLPK.Optimizer)

@variable(shortest_path, x[1:n, 1:n], Bin)
# Arcs with zero cost are not a part of the path as they do no exist
@constraint(shortest_path, [i = 1:n, j = 1:n; G[i, j] == 0], x[i, j] == 0)
# Flow conservation constraint
@constraint(
    shortest_path,
    [i = 1:n; i != 1 && i != 2],
    sum(x[i, :]) == sum(x[:, i])
)
# Flow coming out of source = 1
@constraint(shortest_path, sum(x[1, :]) - sum(x[:, 1]) == 1)
# Flowing coming out of destination = -1 i.e. Flow entering destination = 1
@constraint(shortest_path, sum(x[2, :]) - sum(x[:, 2]) == -1)
@objective(shortest_path, Min, LinearAlgebra.dot(G, x))

optimize!(shortest_path)
objective_value(shortest_path)
#-
value.(x)

# ## The Assignment Problem

# Suppose that there are $n$ persons and $n$ objects that we have to match on a
# one-to-one basis. There is a benefit or value $a_{i,j}$ for matching person
# $i$ with object $j$, and we want to assign persons to objects so as to
# maximize the total benefit.

# There is also a restriction that person $i$ can be assigned to object $j$ only
# if $(i, j)$ belongs to a given set of pairs $A$.

# Mathematically, we want to find a set of person-object pairs
# $(1, j_{1}),..., (n, j_{n})$ from $A$ such that the objects $j_{1},...,j_{n}$
# are all distinct, and the total benefit $\sum_{i=1}^{y} a_{ij_{i}}$ is
# maximized.

# ```math
# \begin{aligned}
# \max && \sum_{(i,j) \in A} a_{i,j} \times y_{i,j} \\
# s.t. && \sum_{\{j|(i,j) \in A\}} y_{i,j} = 1 && \forall i = \{1,2....n\} \\
# && \sum_{\{i|(i,j) \in A\}} y_{i,j} = 1 && \forall j = \{1,2....n\} \\
# && y_{i,j} \in \{0,1\} && \forall (i,j) \in \{1,2...k\}
# \end{aligned}
# ```

# ![Flow Network 2](../../assets/g2.svg)

G = [
    6 4 5 0
    0 3 6 0
    5 0 4 3
    7 5 5 5
]

n = size(G)[1]

assignment = Model(GLPK.Optimizer)
@variable(assignment, y[1:n, 1:n], Bin)
# One person can only be assigned to one object
@constraint(assignment, [i = 1:n], sum(y[:, i]) == 1)
# One object can only be assigned to one person
@constraint(assignment, [j = 1:n], sum(y[j, :]) == 1)
@objective(assignment, Max, LinearAlgebra.dot(G, y))

optimize!(assignment)
objective_value(assignment)
#-
value.(y)

# ## The Max-Flow Problem

# In the max-flow problem, we have a graph with two special nodes: the $source$,
# denoted by $s$, and the $sink$, denoted by $t$.

# The objective is to move as much flow as possible from $s$ into $t$ while
# observing the capacity constraints.

# ```math
# \begin{aligned}
# \max && \sum_{v:(s,v) \in E} f(s,v) \\
# s.t. && \sum_{u:(u,v) \in E} f(u,v)  = \sum_{w:(v,w) \in E} f(v,w) && \forall v \in V - \{s,t\} \\
# && f(u,v) \leq c(u,v) && \forall (u,v) \in E \\
# && f(u,v) \geq 0 && \forall (u,v) \in E
# \end{aligned}
# ```

# ![Flow Network 3](../../assets/g3.svg)

G = [
    0 3 2 2 0 0 0 0
    0 0 0 0 5 1 0 0
    0 0 0 0 1 3 1 0
    0 0 0 0 0 1 0 0
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 2
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 0
]

n = size(G)[1]

max_flow = Model(GLPK.Optimizer)

@variable(max_flow, f[1:n, 1:n] >= 0)
# Capacity constraints
@constraint(max_flow, [i = 1:n, j = 1:n], f[i, j] <= G[i, j])
# Flow conservation constraints
@constraint(max_flow, [i = 1:n; i != 1 && i != 8], sum(f[i, :]) == sum(f[:, i]))
@objective(max_flow, Max, sum(f[1, :]))

optimize!(max_flow)
objective_value(max_flow)
#-
value.(f)
