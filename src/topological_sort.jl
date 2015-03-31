# Modified from Graphs.jl.
# Workaround for edge_index not being O(1) on SimpleGraph.
# edge_index was only needed to test for cycles, so
# this implementation skips that check.

# Graphs.jl is licensed under the MIT License:
#
# Copyright (c) 2012: John Myles White and other contributors.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

type DepthFirst <: Graphs.AbstractGraphVisitAlgorithm
end

type TopologicalSortVisitor <: Graphs.AbstractGraphVisitor
    vertices::Vector{Int}
    parents::Vector{Int}

    function TopologicalSortVisitor(n::Int)
        vs = Array(Int, 0)
        sizehint!(vs, n)
        new(vs, zeros(Int,n))
    end
end

function depth_first_visit_impl!{V,E}(
    graph::AbstractGraph{V,E},      # the graph
    vertex_stack,                   # an (initialized) stack of vertex
    index_stack,                    # stack with out edge indices
    vertexcolormap::Vector{Int},    # an (initialized) color-map to indicate status of vertices
    visitor::TopologicalSortVisitor)  # the visitor

    while !isempty(vertex_stack)
        u = pop!(vertex_stack)
        out_idx = pop!(index_stack)
        uegs = out_edges(u,graph)
        found_new_vertex = false

        while out_idx <= length(uegs) && !found_new_vertex
            v_edge = uegs[out_idx]
            out_idx += 1
            v = v_edge.target
            v_color = vertexcolormap[vertex_index(v, graph)]

            if v_color == 0
                found_new_vertex = true
                vertexcolormap[vertex_index(v, graph)] = 1
                if !discover_vertex!(visitor, v)
                    return
                end
                push!(vertex_stack, u)
                push!(index_stack, out_idx)
                visitor.parents[v] = u

                open_vertex!(visitor, v)
                vegs = out_edges(v, graph)
                push!(vertex_stack, v)
                push!(index_stack, 1)
            end
        end

        if !found_new_vertex
            close_vertex!(visitor, u)
            vertexcolormap[vertex_index(u, graph)] = 2
        end
    end
end

function traverse_graph{V}(
    graph::SimpleGraph,
    alg::DepthFirst,
    s::V,
    visitor::AbstractGraphVisitor,
    vertexcolormap,
    vertex_stack,
    index_stack)

    @graph_requires graph incidence_list vertex_map

    vertexcolormap[vertex_index(s, graph)] = 1
    if !discover_vertex!(visitor, s)
        return
    end

    resize!(vertex_stack,1)
    vertex_stack[1] = s
    resize!(index_stack,1)
    index_stack[1] = 1

    depth_first_visit_impl!(graph, vertex_stack, index_stack, vertexcolormap, visitor)
end



function close_vertex!(visitor::TopologicalSortVisitor, v::Int)
    push!(visitor.vertices, v)
end

function reverse_topological_sort_by_dfs{V}(graph::Graphs.AbstractGraph{V}, cmap::Vector{Int}, vertex_stack::Vector{Int}=Int[], index_stack::Vector{Int}=Int[])
    @graph_requires graph vertex_list incidence_list vertex_map

    @assert length(cmap) == num_vertices(graph)
    fill!(cmap,0)
    visitor = TopologicalSortVisitor(num_vertices(graph))

    for s in vertices(graph)
        if cmap[vertex_index(s, graph)] == 0
            traverse_graph(graph, DepthFirst(), s, visitor, cmap, vertex_stack, index_stack)
        end
    end

    visitor.vertices,visitor.parents
end
