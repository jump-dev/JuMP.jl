using Base.Test
using LightGraphs

function to_adjlist(g::Graph)
    I = Int[]
    J = Int[]
    for e in edges(g)
        push!(I,src(e))
        push!(J,dst(e))
    end
    return gen_adjlist(I,J,length(vertices(g)))
end

import ReverseDiffSparse.Coloring: acyclic_coloring, recovery_preprocess, reverse_topological_sort_by_dfs, gen_adjlist, hessian_color_preprocess, prepare_seed_matrix!, recover_from_matmat!, seed_matrix

# tests for acyclic coloring

g = Graph(10)
color, numcolors = acyclic_coloring(to_adjlist(g))
@test numcolors == 1

add_edge!(g, 2, 4)
color, numcolors = acyclic_coloring(to_adjlist(g))
@test numcolors == 2
@test color[2] != color[4]

add_edge!(g, 2, 3)
color, numcolors = acyclic_coloring(to_adjlist(g))
@test numcolors == 2
@test color[3] == color[4]

add_edge!(g, 3, 4)
color, numcolors = acyclic_coloring(to_adjlist(g))
@test numcolors == 3
recovery_preprocess(to_adjlist(g), color, numcolors, Int[])

g = Graph(3)
add_edge!(g, 1, 3)
add_edge!(g, 2, 3)
color, numcolors = acyclic_coloring(to_adjlist(g))
@test numcolors == 2

g = Graph(4)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 1)
color, numcolors = acyclic_coloring(to_adjlist(g))
@test numcolors == 3

# test our topological sort method
#=
g = Graph(6)
add_edge!(g, 1,2)
add_edge!(g, 1,3)
add_edge!(g, 1,6)
add_edge!(g, 2,4)
add_edge!(g, 2,5)
=#
vec = [3,6,2,1,4,5,1,2,2,1]
offset = [1,4,7,8,9,10,11]

v = reverse_topological_sort_by_dfs(vec, offset, 6, zeros(Int,6))
@test v[1] == [3,6,4,5,2,1]
@test v[2] == [0,1,1,2,2,1]

I,J,rinfo = hessian_color_preprocess(Set([(1,2)]),2)
num_colors = rinfo.num_colors
N = length(rinfo.color)
R = seed_matrix(rinfo)
prepare_seed_matrix!(R, rinfo)
@test I == [1,2,2]
@test J == [1,2,1]
@test R == [1.0 0.0; 0.0 1.0]
hess = [3.4 2.1; 2.1 1.3]
matmat = hess*R
V = zeros(3)
recover_from_matmat!(V, matmat, rinfo, zeros(3))
@test V == [3.4,1.3,2.1]
