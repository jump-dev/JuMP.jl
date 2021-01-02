using Test

import JuMP._Derivatives.Coloring:
    acyclic_coloring,
    recovery_preprocess,
    reverse_topological_sort_by_dfs,
    gen_adjlist,
    hessian_color_preprocess,
    prepare_seed_matrix!,
    recover_from_matmat!,
    seed_matrix

struct Graph
    num_vertices::Int
    edges::Vector{Tuple{Int,Int}}
end

function to_adjlist(graph::Graph)
    I = [i for (i, j) in graph.edges]
    J = [j for (i, j) in graph.edges]
    return gen_adjlist(I, J, graph.num_vertices)
end

# tests for acyclic coloring
@testset "Derivatives (coloring)" begin
    @testset "Edge-free graph" begin
        graph = Graph(10, [])
        color, numcolors = acyclic_coloring(to_adjlist(graph))
        @test numcolors == 1
    end

    @testset "One-edge graph" begin
        graph = Graph(10, [(2, 4)])
        color, numcolors = acyclic_coloring(to_adjlist(graph))
        @test numcolors == 2
        @test color[2] != color[4]
    end

    @testset "Two-edge graph" begin
        graph = Graph(10, [(2, 4), (2, 3)])
        color, numcolors = acyclic_coloring(to_adjlist(graph))
        @test numcolors == 2
        @test color[3] == color[4]
    end

    @testset "Three-edge graph" begin
        graph = Graph(10, [(2, 4), (2, 3), (3, 4)])
        color, numcolors = acyclic_coloring(to_adjlist(graph))
        @test numcolors == 3
        # TODO: What is this testing?
        recovery_preprocess(to_adjlist(graph), color, numcolors, Int[])
    end

    @testset "Two-edge three-vertex graph" begin
        graph = Graph(3, [(1, 3), (2, 3)])
        color, numcolors = acyclic_coloring(to_adjlist(graph))
        @test numcolors == 2
    end

    @testset "Four-edge four-vertex graph" begin
        graph = Graph(4, [(1, 2), (2, 3), (3, 4), (4, 1)])
        color, numcolors = acyclic_coloring(to_adjlist(graph))
        @test numcolors == 3
    end

    @testset "Topological sort" begin
        # graph = Graph(6, [(1, 2), (1, 3), (1, 6), (2, 4), (2, 5)])
        vec = [3, 6, 2, 1, 4, 5, 1, 2, 2, 1]
        offset = [1, 4, 7, 8, 9, 10, 11]

        v = reverse_topological_sort_by_dfs(vec, offset, 6, zeros(Int, 6))
        @test v[1] == [3, 6, 4, 5, 2, 1]
        @test v[2] == [0, 1, 1, 2, 2, 1]
    end

    @testset "End-to-end hessian coloring and recovery" begin
        I, J, rinfo = hessian_color_preprocess(Set([(1, 2)]), 2)
        num_colors = rinfo.num_colors
        N = length(rinfo.color)
        R = seed_matrix(rinfo)
        prepare_seed_matrix!(R, rinfo)
        @test I == [1, 2, 2]
        @test J == [1, 2, 1]
        @test R == [1.0 0.0; 0.0 1.0]
        hess = [3.4 2.1; 2.1 1.3]
        matmat = hess * R
        V = zeros(3)
        recover_from_matmat!(V, matmat, rinfo, zeros(3))
        @test V == [3.4, 1.3, 2.1]
    end
end
