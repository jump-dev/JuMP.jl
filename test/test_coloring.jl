using Base.Test
using ReverseDiffSparse
using Graphs

# tests for acyclic coloring

g = simple_graph(10, is_directed=false)
color, numcolors = ReverseDiffSparse.acyclic_coloring(g)
@test numcolors == 1

add_edge!(g, 2, 4)
color, numcolors = ReverseDiffSparse.acyclic_coloring(g)
@test numcolors == 2
@test color[2] != color[4]

add_edge!(g, 2, 3)
color, numcolors = ReverseDiffSparse.acyclic_coloring(g)
@test numcolors == 2
@test color[3] == color[4]

add_edge!(g, 3, 4)
color, numcolors = ReverseDiffSparse.acyclic_coloring(g)
@test numcolors == 3
ReverseDiffSparse.recovery_preprocess(g, color, verify_acyclic=true)

g = simple_graph(3, is_directed=false)
add_edge!(g, 1, 3)
add_edge!(g, 2, 3)
color, numcolors = ReverseDiffSparse.acyclic_coloring(g)
@test numcolors == 2

g = simple_graph(4, is_directed=false)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 1)
color, numcolors = ReverseDiffSparse.acyclic_coloring(g)
@test numcolors == 3

# test our topological sort method
g = simple_graph(6, is_directed=false)
add_edge!(g, 1,2)
add_edge!(g, 1,3)
add_edge!(g, 1,6)
add_edge!(g, 2,4)
add_edge!(g, 2,5)

v = ReverseDiffSparse.reverse_topological_sort_by_dfs(g, zeros(Int,num_vertices(g)))
@test reverse(v[1]) == [1,6,3,2,5,4]
@test v[2] == [0,1,1,2,2,1]

println("Passed tests")
