using Graphs
import DataStructures

include("topological_sort.jl")

function gen_adjlist(IJ,nel)
    edges = Edge{Int}[]
    sizehint!(edges,round(Int,length(IJ)/2))
    for (i,j) in IJ
        i == j && continue
        push!(edges,Edge(length(edges)+1,i,j))
    end
    g = graph([1:nel;], edges, is_directed=false)
    return g

end

export gen_adjlist

# workaround for slow tuples
immutable MyPair{T}
    first::T
    second::T
end

# workaround for julia issue #10208
Base.hash(x::MyPair{Int},h::UInt) = hash(x.first,hash(x.second,h))

# convert to lower triangular indices, using Pairs
# TODO: replace normalize in hessian.jl
#normalize(i,j) = (j > i) ? (j,i) : (i,j)
normalize_p(p::MyPair) = (p.second > p.first) ? MyPair(p.second,p.first) : MyPair(p.first,p.second)
normalize_p(i,j) = normalize_p(MyPair(i,j))

macro colored(i)
    esc(:((color[$i] != 0)))
end

function prevent_cycle(eg,eg2,S,firstVisitToTree,forbiddenColors,color)
    v = source(eg)
    w = target(eg)
    x = target(eg2)

    er = DataStructures.find_root(S, edge_index(eg2))
    @inbounds first = firstVisitToTree[er]
    p = source(first) # but this depends on the order?
    q = target(first)
    @inbounds if p != v
        firstVisitToTree[er] = eg
    elseif q != w
        forbiddenColors[color[x]] = v
    end
    nothing
end

function grow_star(eg,firstNeighbor,color,S)
    v = source(eg)
    w = target(eg)
    @inbounds e = firstNeighbor[color[w]]
    p = source(e)
    q = target(e)
    @inbounds if p != v
        firstNeighbor[color[w]] = eg
    else
        union!(S, edge_index(eg), edge_index(e))
    end
    nothing
end

function merge_trees(eg,eg1,S)
    e1 = DataStructures.find_root(S, edge_index(eg))
    e2 = DataStructures.find_root(S, edge_index(eg1))
    if e1 != e2
        union!(S, edge_index(eg), edge_index(eg1))
    end
    nothing
end

# acyclic coloring algorithm of Gebremdehin, Tarafdar, Manne, and Pothen
# "New Acyclic and Star Coloring Algorithms with Application to Computing Hessians"
# SIAM J. Sci. Comput. 2007
function acyclic_coloring(g)
    @assert !is_directed(g)
    # Caution: graph is undirected so edges in both directions have
    # the same index but source/target flipped.
    
    num_colors = 0
    forbiddenColors = Int[]
    firstNeighbor = Array(Edge{Int},0)
    firstVisitToTree = fill(Edge{Int}(0,0,0),num_edges(g))
    color = fill(0, num_vertices(g))
    # disjoint set forest of edges in the graph
    S = DataStructures.IntDisjointSets(num_edges(g))

    @inbounds for v in 1:num_vertices(g)
        for eg in out_edges(v,g)
            w = target(eg)
            @colored(w) || continue
            forbiddenColors[color[w]] = v
        end
        for eg in out_edges(v,g)
            w = target(eg)
            @colored(w) || continue
            for eg2 in out_edges(w,g)
                x = target(eg2)
                @colored(x) || continue
                if forbiddenColors[color[x]] != v
                    prevent_cycle(eg,eg2,S,firstVisitToTree,forbiddenColors,color)
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
            push!(firstNeighbor, Edge{Int}(0,0,0))
            color[v] = num_colors
        end

        for eg in out_edges(v,g)
            w = target(eg)
            @colored(w) || continue
            grow_star(eg,firstNeighbor,color,S)
        end
        for eg in out_edges(v,g)
            w = target(eg)
            @colored(w) || continue
            for eg2 in out_edges(w,g)
                x = target(eg2)
                (@colored(x) && x != v) || continue
                if color[x] == color[v]
                    merge_trees(eg,eg2,S)
                end
            end
        end
    end

    return color, num_colors
end

immutable RecoveryInfo
    vertexmap::Vector{Vector{Int}}
    postorder::Vector{Vector{Int}}
    parents::Vector{Vector{Int}}
    color::Vector{Int}
end

function recovery_preprocess(g,color,num_colors)
    # represent two-color subgraph as:
    # list of vertices (with map to global indices)
    # adjacency list in a single vector (with list of offsets)


    # linear index of pair of colors
    twocolorindex = zeros(Int32,num_colors, num_colors)
    seen_twocolors = 0
    g_edges = edges(g)
    # count of edges in each subgraph
    edge_count = Array(Int,0)
    for k in 1:length(g_edges)
        e = g_edges[k]
        u = source(e,g)
        v = target(e,g)
        i = min(color[u],color[v])
        j = max(color[u],color[v])
        if twocolorindex[i,j] == 0
            seen_twocolors += 1
            twocolorindex[i,j] = seen_twocolors
            push!(edge_count,0)
        end
        idx = twocolorindex[i,j]
        edge_count[idx] += 1
    end
    # edges sorted by twocolor subgraph
    sorted_edges = Array(Vector{MyPair{Int}},seen_twocolors)
    for idx in 1:seen_twocolors
        sorted_edges[idx] = Array(MyPair{Int},0)
        sizehint!(sorted_edges[idx],edge_count[idx])
    end

    for i in 1:length(g_edges)
        e = g_edges[i]
        u = source(e,g)
        v = target(e,g)
        i = min(color[u],color[v])
        j = max(color[u],color[v])
        idx = twocolorindex[i,j]
        push!(sorted_edges[idx], MyPair(u,v))
    end

    # list of unique vertices in each twocolor subgraph
    vertexmap = Array(Vector{Int},seen_twocolors)

    postorder = Array(Vector{Int},seen_twocolors)
    parents = Array(Vector{Int},seen_twocolors)

    # temporary lookup map from global index to subgraph index
    revmap = zeros(Int,num_vertices(g))

    adjcount = zeros(Int,num_vertices(g))

    cmap = zeros(Int,0) # shared storage for DFS
    vertex_stack = Int[]
    index_stack = Int[]

    for idx in 1:seen_twocolors
        my_edges = sorted_edges[idx]
        vlist = Int[]

        # build up the vertex list and adjacency count
        for k in 1:length(my_edges)
            e = my_edges[k]
            u = e.first
            v = e.second
            # seen these vertices yet?
            if revmap[u] == 0
                push!(vlist,u)
                revmap[u] = length(vlist)
                adjcount[u] = 0
            end
            if revmap[v] == 0
                push!(vlist,v)
                revmap[v] = length(vlist)
                adjcount[v] = 0
            end
            adjcount[u] += 1
            adjcount[v] += 1
        end

        # set up offsets for adjlist
        offset = Array(Int,length(vlist)+1)
        offset[1] = 1
        for k in 1:length(vlist)
            offset[k+1] = offset[k] + adjcount[vlist[k]]
            adjcount[vlist[k]] = 0
        end
        # adjlist for node u in twocolor idx starts at
        # vec[offset[u]]
        # u has global index vlist[u]
        vec = Array(Int,offset[length(vlist)+1]-1)

        # now fill in
        for k in 1:length(my_edges)
            e = my_edges[k]
            u = e.first
            v = e.second

            u_rev = revmap[u] # indices in the subgraph
            v_rev = revmap[v]
            vec[offset[u_rev]+adjcount[u]] = v_rev
            vec[offset[v_rev]+adjcount[v]] = u_rev

            adjcount[u] += 1
            adjcount[v] += 1
        end

        resize!(cmap,length(vlist))
        order, parent = reverse_topological_sort_by_dfs(vec,offset,length(vlist),cmap)

        for k in 1:length(vlist)
            # clear for reuse
            revmap[vlist[k]] = 0
        end

        postorder[idx] = order
        parents[idx] = parent
        vertexmap[idx] = vlist
    end

    return RecoveryInfo(vertexmap, postorder, parents, color)

end

function indirect_recover_structure(nnz, rinfo::RecoveryInfo)
    N = length(rinfo.color)
    
    I = zeros(Int, nnz+N)
    J = zeros(Int, nnz+N)
    
    # diagonal entries
    k = 0
    for i in 1:N
        k += 1
        I[k] = i
        J[k] = i
    end

    for t in 1:length(rinfo.postorder)
        vmap = rinfo.vertexmap[t]
        order = rinfo.postorder[t]
        parent = rinfo.parents[t]

        for z in 1:length(order)
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

    end

    @assert k == nnz + N

    return I,J
end

function indirect_recover(hessian_matmat!, nnz, rinfo::RecoveryInfo, stored_values, x, inputvals, fromcanonical, R, dualvec, dualout, V)
    N = length(rinfo.color)
    
    #R = zeros(N,num_colors)
    fill!(R,0.0)
    for i in 1:N
        R[i,rinfo.color[i]] = 1
    end

    hessian_matmat!(R,x, dualvec, dualout, inputvals, fromcanonical)
    
    # now, recover
    @assert length(V) == nnz+N
    
    # diagonal entries
    k = 0

    for i in 1:N
        k += 1
        V[k] = R[i,rinfo.color[i]]
    end

    for t in 1:length(rinfo.vertexmap)
        vmap = rinfo.vertexmap[t]
        order = rinfo.postorder[t]
        parent = rinfo.parents[t]
        stored_values[1:length(order)] = 0.0

        @inbounds for z in 1:length(order)
            v = order[z]
            p = parent[v]
            (p == 0) && continue
            
            i = vmap[v]
            j = vmap[p]

            value = R[i,rinfo.color[j]] - stored_values[v]
            stored_values[p] += value

            k += 1
            V[k] = value
        end
    end

    @assert k == nnz + N

    return V

end

export acyclic_coloring, indirect_recover

gen_hessian_sparse_color_parametric(s::SymbolicOutput, num_total_vars) =
    gen_hessian_sparse_color_parametric(s,num_total_vars,gen_hessian_matmat_parametric(s),compute_hessian_sparsity_IJ_parametric(s))

function gen_hessian_sparse_color_parametric(s::SymbolicOutput, num_total_vars, hessian_matmat!, hessian_IJ, dualvec=Array(Dual{Float64}, num_total_vars), dualout=Array(Dual{Float64}, num_total_vars))
    I,J = hessian_IJ(s)
    # remove duplicates
    M = sparse(I,J,ones(length(I)))
    I,J = findn(M)
    if length(I) == 0
        # expression is actually linear, return dummy function
        return I,J, (x,output_values,ex) -> nothing
    end


    g = gen_adjlist(zip(I,J), length(s.mapfromcanonical))

    color, num_colors = acyclic_coloring(g)

    @assert length(color) == num_vertices(g)

    R = Array(Float64,num_vertices(g),num_colors)
    
    rinfo = recovery_preprocess(g, color, num_colors)

    stored_values = Array(Float64,num_vertices(g))

    I,J = indirect_recover_structure(num_edges(g), rinfo)

    nnz = num_edges(g)
    
    function eval_h(x,output_values, ex::SymbolicOutput)
        indirect_recover(hessian_matmat!, nnz, rinfo, stored_values, x, ex.inputvals, ex.mapfromcanonical, R, dualvec, dualout, output_values)
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
