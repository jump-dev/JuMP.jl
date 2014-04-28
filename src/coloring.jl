using Graphs
using DataStructures


function gen_adjlist(IJ,nel)
    g = simple_graph(nel, is_directed=false)
    for (i,j) in IJ
        i == j && continue
        add_edge!(g, i, j)
    end
    return g

end

export gen_adjlist

function reverse_dict_lookup(d, v)
    for k in keys(d)
        if d[k] == v
            return k
        end
    end
    error("Not found")
end

# acyclic coloring algorithm of Gebremdehin, Tarafdar, Manne, and Pothen
# "New Acyclic and Star Coloring Algorithms with Application to Computing Hessians"
# SIAM J. Sci. Comput. 2007
function acyclic_coloring(g)
    
    num_colors = 0
    forbiddenColors = Int[]
    firstNeighbor = Array((Int,Int),0)
    firstVisitToTree = [normalize(source(e),target(e)) => (0,0) for e in edges(g)]
    color = fill(0, num_vertices(g))
    colored(i) = (color[i] != 0)
    # disjoint set forest of edges in the graph
    S = DisjointSets{(Int,Int)}(collect(keys(firstVisitToTree)))

    function prevent_cycle(v,w,x)
        er = find_root(S, normalize(w,x))
        # reverse lookup, this isn't ideal
        e = reverse_dict_lookup(S.intmap,er)
        p,q = firstVisitToTree[e]
        if p != v
            firstVisitToTree[e] = normalize(v,w)
        elseif q != w
            forbiddenColors[color[x]] = v
        end
    end

    function grow_star(v,w)
        p,q = firstNeighbor[color[w]]
        if p != v
            firstNeighbor[color[w]] = normalize(v,w)
        else
            union!(S, normalize(v,w), normalize(p,q))
        end
    end

    function merge_trees(v,w,x)
        e1 = find_root(S, normalize(v,w))
        e2 = find_root(S, normalize(w,x))
        if e1 != e2
            union!(S, normalize(v,w), normalize(w,x))
        end
    end

    for v in 1:num_vertices(g)
        for w in out_neighbors(v,g)
            colored(w) || continue
            forbiddenColors[color[w]] = v
        end
        for w in out_neighbors(v,g)
            colored(w) || continue
            for x in out_neighbors(w,g)
                colored(x) || continue
                if forbiddenColors[color[x]] != v
                    prevent_cycle(v,w,x)
                end
            end
        end

        # find feasible color
        found = false
        for k in 1:num_colors
            if forbiddenColors[k] != v
                color[v] = k
                found = true
                break
            end
        end
        if !found
            num_colors += 1
            push!(forbiddenColors, 0)
            push!(firstNeighbor, (0, 0))
            color[v] = num_colors
        end

        for w in out_neighbors(v,g)
            colored(w) || continue
            grow_star(v,w)
        end
        for w in out_neighbors(v,g)
            colored(w) || continue
            for x in out_neighbors(w,g)
                (colored(x) && x != v) || continue
                if color[x] == color[v]
                    merge_trees(v,w,x)
                end
            end
        end
    end

    return color, num_colors
end

if VERSION >= v"0.3.0-"
    function twocolorset_of_edge(e,g,color)
        i = source(e,g)
        j = target(e,g)
        return Set([color[i], color[j]])
    end
else
    function twocolorset_of_edge(e,g,color)
        i = source(e,g)
        j = target(e,g)
        return Set(color[i], color[j])
    end
end

function recovery_preprocess(g,color)
    twocoloredges = Dict{Set{Int},Vector{(Int,Int)}}()
    twocolorvertices = Dict{Set{Int},Set{Int}}()
    for e in edges(g)
        twocolor = twocolorset_of_edge(e,g,color)
        if !haskey(twocoloredges, twocolor)
            twocoloredges[twocolor] = Array((Int,Int),0)
            twocolorvertices[twocolor] = Set{Int}()
        end
        push!(twocoloredges[twocolor], (source(e,g),target(e,g)))
        push!(twocolorvertices[twocolor], source(e,g))
        push!(twocolorvertices[twocolor], target(e,g))
    end

    twocolorgraphs = SimpleGraph[]
    # map from vertices in two-color subgraphs to original vertices
    vertexmap = Array(Vector{Int},0)
    revmap = zeros(Int,num_vertices(g))
    for twocolor in keys(twocoloredges)
        edgeset = twocoloredges[twocolor]
        vertexset = twocolorvertices[twocolor]
        n = length(vertexset)
        s = simple_graph(n, is_directed=false)

        vmap = Int[]
        for v in vertexset
            push!(vmap, v)
            revmap[v] = length(vmap)
        end

        for (i,j) in edgeset
            add_edge!(s, revmap[i], revmap[j])
        end

        push!(twocolorgraphs, s)
        push!(vertexmap, vmap)
    end

    # list the vertices in postorder
    postorder = [reverse_topological_sort_by_dfs(s) for s in twocolorgraphs]
    # identify each vertex's parent in the tree
    parents = Array(Vector{Int},0)
    for i in 1:length(twocolorgraphs)
        s = twocolorgraphs[i]
        parent = zeros(num_vertices(s))
        s = twocolorgraphs[i]
        seen = falses(num_vertices(s))
        for k in 1:num_vertices(s)
            v = postorder[i][k]
            seen[v] = true
            # find the neighbor that we haven't seen
            notseen = 0
            numseen = 0
            for w in out_neighbors(v,s)
                if seen[w]
                    numseen += 1
                else
                    notseen = w
                end
            end
            if numseen == length(out_neighbors(v,s)) - 1
                parent[v] = notseen
            else
                (numseen == length(out_neighbors(v,s))) || error("Error processing tree, invalid ordering")
            end
        end
        push!(parents, parent)
    end


    return (twocolorgraphs, vertexmap, postorder, parents)

end

function indirect_recover(hessian_matmat!, nnz, twocolorgraphs, vertexmap, postorder, parents, stored_values, color, num_colors, x, inputvals, fromcanonical, tocanonical, R, dualvec, dualout, V; structure=false)
    N = length(color)
    
    # generate vectors for hessian-vec product
    #R = zeros(N,num_colors)
    fill!(R,0.0)
    for i in 1:N
        R[i,color[i]] = 1
    end

    if !structure
        hessian_matmat!(R,x, dualvec, dualout, inputvals, fromcanonical, tocanonical)
    end
    
    # now, recover
    if structure
        I = zeros(Int, nnz+N)
        J = zeros(Int, nnz+N)
    else
        I = Int[]
        J = Int[]
        @assert length(V) == nnz+N
    end
    
    # diagonal entries
    k = 0
    if structure
        for i in 1:N
            k += 1
            I[k] = i
            J[k] = i
        end
    else
        for i in 1:N
            k += 1
            V[k] = R[i,color[i]]
        end
    end

    for t in 1:length(twocolorgraphs)
        s = twocolorgraphs[t]
        vmap = vertexmap[t]
        order = postorder[t]
        parent = parents[t]
        stored_values[1:num_vertices(s)] = 0.0

        if structure
            for z in 1:num_vertices(s)
                v = order[z]
                p = parent[v]
                (p == 0) && continue
                i = vmap[v]
                j = vmap[p]
                i,j = normalize(i,j)
                k += 1
                I[k] = i
                J[k] = j
            end
        else
            for z in 1:num_vertices(s)
                v = order[z]
                p = parent[v]
                (p == 0) && continue
                
                i = vmap[v]
                j = vmap[p]

                value = R[i,color[j]] - stored_values[v]
                stored_values[p] += value

                k += 1
                V[k] = value
            end
        end
    end

    
    

    @assert k == nnz + N

    if structure
        return I,J
    else
        return V
    end

end

export acyclic_coloring, indirect_recover

function gen_hessian_sparse_color_parametric(s::SymbolicOutput, num_total_vars)
    I,J = compute_hessian_sparsity_IJ(s)
    # remove duplicates
    M = sparse(I,J,ones(length(I)))
    I,J = findn(M)
    if length(I) == 0
        # expression is actually linear, return dummy function
        return I,J, (x,output_values,ex) -> nothing
    end


    hessian_matmat! = gen_hessian_matmat_parametric(s)
    
    g = gen_adjlist(zip(I,J), length(s.mapfromcanonical))
    
    color, num_colors = acyclic_coloring(g)

    @assert length(color) == num_vertices(g)

    R = Array(Float64,num_vertices(g),num_colors)
    
    (twocolorgraphs, vertexmap, postorder, parents) = recovery_preprocess(g, color)

    stored_values = Array(Float64,num_vertices(g))
    dualvec = Array(Dual{Float64}, num_total_vars)
    dualout = Array(Dual{Float64}, num_total_vars)

    I,J = indirect_recover(hessian_matmat!, num_edges(g), twocolorgraphs, vertexmap, postorder, parents, stored_values, color, num_colors, nothing, s.inputvals, s.mapfromcanonical, s.maptocanonical, R, dualvec, dualout, Float64[]; structure=true)
    
    function eval_h(x,output_values, ex::SymbolicOutput)
        indirect_recover(hessian_matmat!, num_edges(g), twocolorgraphs, vertexmap, postorder, parents, stored_values, color, num_colors, x, ex.inputvals, ex.mapfromcanonical, ex.maptocanonical, R, dualvec, dualout, output_values)
    end

    return I,J, eval_h


end

export gen_hessian_sparse_color_parametric

function to_H(s::SymbolicOutput, I, J, V, n)
    I2 = similar(I)
    J2 = similar(J)
    for k in 1:length(I)
        I2[k],J2[k] = normalize(s.mapfromcanonical[I[k]],s.mapfromcanonical[J[k]])
    end
    return sparse(I2,J2,V, n, n)

end

export to_H


# Topological sort using DFS -- copied from Graphs.jl until 0.3 is released...

type MyTopologicalSortVisitor <: AbstractGraphVisitor
    vertices::Vector{Int}

    function MyTopologicalSortVisitor(n::Int)
        vs = Array(Int, 0)
        sizehint(vs, n)
        new(vs)
    end
end

import Graphs: examine_neighbor!, close_vertex!


function examine_neighbor!{V}(visitor::MyTopologicalSortVisitor, u::V, v::V, vcolor::Int, ecolor::Int)
    if vcolor == 1 && ecolor == 0
        throw(ArgumentError("The input graph contains at least one loop."))
    end
end

function close_vertex!{V}(visitor::MyTopologicalSortVisitor, v::V)
    push!(visitor.vertices, v)
end

function reverse_topological_sort_by_dfs{V}(graph::AbstractGraph{V})

    cmap = zeros(Int, num_vertices(graph))
    visitor = MyTopologicalSortVisitor(num_vertices(graph))

    for s in vertices(graph)
        if cmap[vertex_index(s, graph)] == 0
            traverse_graph(graph, DepthFirst(), s, visitor, vertexcolormap=cmap)
        end
    end

    visitor.vertices
end


