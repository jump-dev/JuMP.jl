using Base.Test
using ReverseDiffSparse
using Graphs

function to_adjlist(g::SimpleGraph)
    I = Int[]
    J = Int[]
    for e in edges(g)
        push!(I,source(e,g))
        push!(J,target(e,g))
    end
    return gen_adjlist(I,J,num_vertices(g))
end

# tests for acyclic coloring

g = simple_graph(10, is_directed=false)
color, numcolors = ReverseDiffSparse.acyclic_coloring(to_adjlist(g))
@test numcolors == 1

add_edge!(g, 2, 4)
color, numcolors = ReverseDiffSparse.acyclic_coloring(to_adjlist(g))
@test numcolors == 2
@test color[2] != color[4]

add_edge!(g, 2, 3)
color, numcolors = ReverseDiffSparse.acyclic_coloring(to_adjlist(g))
@test numcolors == 2
@test color[3] == color[4]

add_edge!(g, 3, 4)
color, numcolors = ReverseDiffSparse.acyclic_coloring(to_adjlist(g))
@test numcolors == 3
ReverseDiffSparse.recovery_preprocess(to_adjlist(g), color, numcolors)

g = simple_graph(3, is_directed=false)
add_edge!(g, 1, 3)
add_edge!(g, 2, 3)
color, numcolors = ReverseDiffSparse.acyclic_coloring(to_adjlist(g))
@test numcolors == 2

g = simple_graph(4, is_directed=false)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 1)
color, numcolors = ReverseDiffSparse.acyclic_coloring(to_adjlist(g))
@test numcolors == 3

# test our topological sort method
#=
g = simple_graph(6, is_directed=false)
add_edge!(g, 1,2)
add_edge!(g, 1,3)
add_edge!(g, 1,6)
add_edge!(g, 2,4)
add_edge!(g, 2,5)
=#
vec = [3,6,2,1,4,5,1,2,2,1]
offset = [1,4,7,8,9,10,11]

v = ReverseDiffSparse.reverse_topological_sort_by_dfs(vec, offset, 6, zeros(Int,6))
@test v[1] == [3,6,4,5,2,1]
@test v[2] == [0,1,1,2,2,1]

println("Passed tests")
