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

function twocolorset_of_edge(e,g,color)
    i = source(e,g)
    j = target(e,g)
    return Set([color[i], color[j]]) # TODO: 0.2 compat
end

function generate_2color_subgraphs(g,color)
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

    return twocolorgraphs, vertexmap

end

type DFSRecoveryVisitor <: AbstractGraphVisitor
    color::Vector{Int} # color of each original vertex
    R::Matrix{Float64} # compressed Hessian product
    I::Vector{Int} # row indices
    J::Vector{Int} # col indices
    V#::AbstractVector{Float64} # output values
    k::Int # nnz counter
    structure::Bool # structure or values
    # these fields change with each subgraph
    g::SimpleGraph # two-color subgraph
    parent::Vector{Int}
    recovered_values::Vector{Float64} # recovered values of edges on this subgraph
    vmap::Vector{Int} # map from two-color vertex indices to original indices

    function DFSRecoveryVisitor(color, R, I, J, V, k, structure)
        vis = new(color, R, I, J, V, k, structure)
        @assert pointer(V) == pointer(vis.V) # julia issue #6617
        vis.parent = zeros(length(color))
        if structure
            vis.recovered_values = zeros(length(I))
        else
            vis.recovered_values = zeros(length(V))
        end
        return vis
    end
end

function initialize!(vis::DFSRecoveryVisitor,newG,vmap)
    vis.g = newG
    vis.vmap = vmap
    E = num_edges(vis.g)
    V = num_vertices(vis.g)
    vis.parent[1:V] = 0
    vis.recovered_values[1:E] = 0
end

import Graphs: close_vertex!, examine_neighbor!

# called when all of v's neighbors have been explored, so this is postorder
function close_vertex!(vis::DFSRecoveryVisitor, v)
    p = vis.parent[vertex_index(v)]
    if p == 0 # root of the tree
        return
    end
    #println("Closing vertex $(vis.vmap[v]), parent $(vis.vmap[p])")
    @assert p > 0
    # we're recovering the value corresponding to the edge that points to the parent
    s = 0.0
    paridx = 0
    for e in out_edges(v,vis.g)
        w = target(e)
        if w != p
            s += vis.recovered_values[edge_index(e)]
        else
            paridx = edge_index(e)
        end
    end
    i = vis.vmap[vertex_index(v)]
    j = vis.vmap[p]
    
    vis.k += 1
    
    if vis.structure
        i,j = normalize(i,j)
        vis.I[vis.k] = i
        vis.J[vis.k] = j
    else
        val = vis.R[i,vis.color[j]] - s
        vis.recovered_values[paridx] = val
        vis.V[vis.k] = val
    end

end

function examine_neighbor!{V}(
    vis::DFSRecoveryVisitor,
    u::V,
    v::V,
    vcolor::Int,
    ecolor::Int)
    
    #println("Exploring edge $(vis.vmap[u]) -> $(vis.vmap[v]): $vcolor $ecolor")
    if vcolor == 1 && ecolor == 0
        error("Invalid input, found cycle in two-color subtree!")
    end
    if ecolor == 1
        return # already saw this edge
    end
    # set u as parent of v
    vis.parent[vertex_index(v)] = vertex_index(u)
end

function indirect_recover(hessian_matmat!, nnz, twocolorgraphs, vertexmap, color, num_colors, x, inputvals, fromcanonical, tocanonical, R, V; structure=false)
    N = length(color)

    # generate vectors for hessian-vec product
    #R = zeros(N,num_colors)
    fill!(R,0.0)
    for i in 1:N
        R[i,color[i]] = 1
    end

    if !structure
        hessian_matmat!(R,x, inputvals, fromcanonical, tocanonical)
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
    for i in 1:N
        k += 1
        if structure
            I[k] = i
            J[k] = i
        else
            V[k] = R[i,color[i]]
        end
    end
    
    vis = DFSRecoveryVisitor(color, R, I, J, V, k, structure)

    vcmap = zeros(Int,N)
    ecmap = zeros(Int,nnz)

    # post-order depth-first search on connected components of each two-color subgraph
    for (s,vmap) in zip(twocolorgraphs, vertexmap)
        #println("G: $s")
        #println(edges(s))
        #println("Vmap: $vmap")
        
        vcmap[1:num_vertices(s)] = 0
        ecmap[1:num_edges(s)] = 0
        initialize!(vis, s, vmap)
        for v in vertices(s)
            if vcmap[vertex_index(v,s)] == 0
                # new connected component
                traverse_graph(s, DepthFirst(), v, vis, vertexcolormap=vcmap, edgecolormap=ecmap)
            end
        end
    end

    @assert vis.k == nnz + N

    if structure
        return I,J
    else
        return V
    end

end

export acyclic_coloring, indirect_recover

function gen_hessian_sparse_color_parametric(s::SymbolicOutput)
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

    R = Array(Float64,num_vertices(g),num_colors)
    
    twocolorgraphs, vertexmap = generate_2color_subgraphs(g, color)

    I,J = indirect_recover(hessian_matmat!, num_edges(g), twocolorgraphs, vertexmap, color, num_colors, nothing, s.inputvals, s.mapfromcanonical, s.maptocanonical, R, Float64[]; structure=true)
    
    function eval_h(x,output_values, ex::SymbolicOutput)
        indirect_recover(hessian_matmat!, num_edges(g), twocolorgraphs, vertexmap, color, num_colors, x, ex.inputvals, ex.mapfromcanonical, ex.maptocanonical, R, output_values)
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





