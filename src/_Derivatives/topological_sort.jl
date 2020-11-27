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

mutable struct TopologicalSortVisitor
    vertices::Vector{Int}
    parents::Vector{Int}

    function TopologicalSortVisitor(n::Int)
        vs = Int[]
        sizehint!(vs, n)
        return new(vs, zeros(Int, n))
    end
end

function depth_first_visit_impl!(
    adjlist,
    offsets,
    vertex_stack,                   # an (initialized) stack of vertex
    index_stack,                    # stack with out edge indices
    vertexcolormap::Vector{Int},    # an (initialized) color-map to indicate status of vertices
    visitor::TopologicalSortVisitor,
)
    while !isempty(vertex_stack)
        u = pop!(vertex_stack)
        out_idx = pop!(index_stack)
        #uegs = out_edges(u,graph)
        len_uegs = offsets[u+1] - offsets[u]
        found_new_vertex = false

        while out_idx <= len_uegs && !found_new_vertex
            v = adjlist[offsets[u]+out_idx-1]
            out_idx += 1
            v_color = vertexcolormap[v]

            if v_color == 0
                found_new_vertex = true
                vertexcolormap[v] = 1
                push!(vertex_stack, u)
                push!(index_stack, out_idx)
                visitor.parents[v] = u

                push!(vertex_stack, v)
                push!(index_stack, 1)
            end
        end

        if !found_new_vertex
            close_vertex!(visitor, u)
            vertexcolormap[u] = 2
        end
    end
end

function traverse_graph(
    adjlist,
    offsets,
    s,
    visitor::TopologicalSortVisitor,
    vertexcolormap,
    vertex_stack,
    index_stack,
)
    vertexcolormap[s] = 1

    resize!(vertex_stack, 1)
    vertex_stack[1] = s
    resize!(index_stack, 1)
    index_stack[1] = 1

    return depth_first_visit_impl!(
        adjlist,
        offsets,
        vertex_stack,
        index_stack,
        vertexcolormap,
        visitor,
    )
end

function close_vertex!(visitor::TopologicalSortVisitor, v::Int)
    return push!(visitor.vertices, v)
end

function reverse_topological_sort_by_dfs(
    adjlist,
    offsets,
    num_vertices,
    cmap::Vector{Int},
    vertex_stack::Vector{Int} = Int[],
    index_stack::Vector{Int} = Int[],
)
    @assert length(cmap) == num_vertices
    fill!(cmap, 0)
    visitor = TopologicalSortVisitor(num_vertices)

    for s in 1:num_vertices
        if cmap[s] == 0
            traverse_graph(
                adjlist,
                offsets,
                s,
                visitor,
                cmap,
                vertex_stack,
                index_stack,
            )
        end
    end

    return visitor.vertices, visitor.parents
end
