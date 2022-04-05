module Coloring

import DataStructures

include("topological_sort.jl")

# workaround for slow tuples
struct MyPair{T}
    first::T
    second::T
end

# workaround for julia issue #10208
Base.hash(x::MyPair{Int}, h::UInt) = hash(x.first, hash(x.second, h))

# indexed sparse set of integers
mutable struct IndexedSet
    nzidx::Vector{Int}
    empty::BitArray{1}
    nnz::Int
end

IndexedSet(n::Integer) = IndexedSet(zeros(Int, n), trues(n), 0)

function Base.push!(v::IndexedSet, i::Integer)
    if v.empty[i]  # new index
        v.nzidx[v.nnz+=1] = i
        v.empty[i] = false
    end
    return nothing
end

function Base.empty!(v::IndexedSet)
    nzidx = v.nzidx
    empty = v.empty
    for i in 1:v.nnz
        empty[nzidx[i]] = true
    end
    return v.nnz = 0
end

Base.length(v::IndexedSet) = length(v.nzidx)
function Base.resize!(v::IndexedSet, n::Integer)
    if n > length(v)
        @assert v.nnz == 0 # only resize empty vector
        resize!(v.nzidx, n)
        resize!(v.empty, n)
        fill!(v.empty, true)
    end
end

Base.collect(v::IndexedSet) = v.nzidx[1:v.nnz]
function Base.union!(v::IndexedSet, s)
    for x in s
        push!(v, x)
    end
    return nothing
end

# compact storage for an undirected graph
# neighbors of vertex i start at adjlist[offsets[i]]
struct UndirectedGraph
    adjlist::Vector{Int}
    edgeindex::Vector{Int} # corresponding edge number, indexed as adjlist
    offsets::Vector{Int}
    edges::Vector{MyPair{Int}}
end
num_vertices(g::UndirectedGraph) = length(g.offsets) - 1
num_edges(g::UndirectedGraph) = length(g.edges)
num_neighbors(i::Int, g::UndirectedGraph) = g.offsets[i+1] - g.offsets[i]
start_neighbors(i::Int, g::UndirectedGraph) = g.offsets[i]

function gen_adjlist(I, J, nel)
    adjcount = zeros(Int, nel)
    n_edges = 0
    for k in 1:length(I)
        i = I[k]
        j = J[k]
        i == j && continue
        n_edges += 1
        adjcount[i] += 1
        adjcount[j] += 1
    end
    offsets = Array{Int}(undef, nel + 1)
    offsets[1] = 1
    for k in 1:nel
        offsets[k+1] = offsets[k] + adjcount[k]
    end
    fill!(adjcount, 0)

    edges = Array{MyPair{Int}}(undef, n_edges)
    adjlist = Array{Int}(undef, offsets[nel+1] - 1)
    edgeindex = Array{Int}(undef, length(adjlist))
    edge_count = 0

    for k in 1:length(I)
        i = I[k]
        j = J[k]
        i == j && continue
        edge_count += 1
        adjlist[offsets[i]+adjcount[i]] = j
        edgeindex[offsets[i]+adjcount[i]] = edge_count
        adjcount[i] += 1

        adjlist[offsets[j]+adjcount[j]] = i
        edgeindex[offsets[j]+adjcount[j]] = edge_count
        adjcount[j] += 1

        edges[edge_count] = MyPair(i, j)
    end
    @assert edge_count == n_edges

    return UndirectedGraph(adjlist, edgeindex, offsets, edges)
end

struct Edge
    index::Int
    source::Int
    target::Int
end

# convert to lower triangular indices, using Pairs
normalize(i, j) = (j > i) ? (j, i) : (i, j)
function normalize_p(p::MyPair)
    return (p.second > p.first) ? MyPair(p.second, p.first) :
           MyPair(p.first, p.second)
end
normalize_p(i, j) = normalize_p(MyPair(i, j))

macro colored(i)
    return esc(:((color[$i] != 0)))
end

function prevent_cycle(
    v,
    w,
    x,
    e_idx1,
    e_idx2,
    S,
    firstVisitToTree,
    forbiddenColors,
    color,
)
    er = DataStructures.find_root!(S, e_idx2)
    @inbounds first = firstVisitToTree[er]
    p = first.source # but this depends on the order?
    q = first.target
    @inbounds if p != v
        firstVisitToTree[er] = Edge(e_idx1, v, w)
    elseif q != w
        forbiddenColors[color[x]] = v
    end
    return nothing
end

function grow_star(v, w, e_idx, firstNeighbor, color, S)
    @inbounds e = firstNeighbor[color[w]]
    p = e.source
    q = e.target
    @inbounds if p != v
        firstNeighbor[color[w]] = Edge(e_idx, v, w)
    else
        union!(S, e_idx, e.index)
    end
    return nothing
end

function merge_trees(eg, eg1, S)
    e1 = DataStructures.find_root!(S, eg)
    e2 = DataStructures.find_root!(S, eg1)
    if e1 != e2
        union!(S, eg, eg1)
    end
    return nothing
end

# acyclic coloring algorithm of Gebremdehin, Tarafdar, Manne, and Pothen
# "New Acyclic and Star Coloring Algorithms with Application to Computing Hessians"
# SIAM J. Sci. Comput. 2007
function acyclic_coloring(g::UndirectedGraph)
    if num_edges(g) == 0
        return fill(1, num_vertices(g)), 1
    end
    num_colors = 0
    forbiddenColors = Int[]
    firstNeighbor = Edge[]
    firstVisitToTree = fill(Edge(0, 0, 0), num_edges(g))
    color = fill(0, num_vertices(g))
    # disjoint set forest of edges in the graph
    S = DataStructures.IntDisjointSets(num_edges(g))

    @inbounds for v in 1:num_vertices(g)
        n_neighbor = num_neighbors(v, g)
        start_neighbor = start_neighbors(v, g)
        for k in 0:(n_neighbor-1)
            w = g.adjlist[start_neighbor+k]
            @colored(w) || continue
            forbiddenColors[color[w]] = v
        end
        for k in 0:(n_neighbor-1)
            w = g.adjlist[start_neighbor+k]
            e_idx = g.edgeindex[start_neighbor+k]
            @colored(w) || continue
            n_neighbor_w = num_neighbors(w, g)
            start_neighbor_w = start_neighbors(w, g)
            for k2 in 0:(n_neighbor_w-1)
                x = g.adjlist[start_neighbor_w+k2]
                e2_idx = g.edgeindex[start_neighbor_w+k2]
                @colored(x) || continue
                if forbiddenColors[color[x]] != v
                    prevent_cycle(
                        v,
                        w,
                        x,
                        e_idx,
                        e2_idx,
                        S,
                        firstVisitToTree,
                        forbiddenColors,
                        color,
                    )
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
            push!(firstNeighbor, Edge(0, 0, 0))
            color[v] = num_colors
        end

        for k in 0:(n_neighbor-1)
            w = g.adjlist[start_neighbor+k]
            e_idx = g.edgeindex[start_neighbor+k]
            @colored(w) || continue
            grow_star(v, w, e_idx, firstNeighbor, color, S)
        end
        for k in 0:(n_neighbor-1)
            w = g.adjlist[start_neighbor+k]
            e_idx = g.edgeindex[start_neighbor+k]
            @colored(w) || continue
            n_neighbor_w = num_neighbors(w, g)
            start_neighbor_w = start_neighbors(w, g)
            for k2 in 0:(n_neighbor_w-1)
                x = g.adjlist[start_neighbor_w+k2]
                e2_idx = g.edgeindex[start_neighbor_w+k2]
                (@colored(x) && x != v) || continue
                if color[x] == color[v]
                    merge_trees(e_idx, e2_idx, S)
                end
            end
        end
    end

    return color, num_colors
end

struct RecoveryInfo
    vertexmap::Vector{Vector{Int}}
    postorder::Vector{Vector{Int}}
    parents::Vector{Vector{Int}}
    color::Vector{Int}
    num_colors::Int
    nnz::Int # number of off-diagonal nonzeros
    local_indices::Vector{Int} # map back to global indices
end

function RecoveryInfo()
    return RecoveryInfo(
        Vector{Vector{Int}}(undef, 0),
        Vector{Vector{Int}}(undef, 0),
        Vector{Vector{Int}}(undef, 0),
        Vector{Int}(undef, 0),
        0,
        0,
        Vector{Int}(undef, 0),
    )
end

function recovery_preprocess(
    g::UndirectedGraph,
    color,
    num_colors,
    local_indices,
)
    # represent two-color subgraph as:
    # list of vertices (with map to global indices)
    # adjacency list in a single vector (with list of offsets)

    # linear index of pair of colors
    twocolorindex = zeros(Int32, num_colors, num_colors)
    seen_twocolors = 0
    # count of edges in each subgraph
    edge_count = Int[]
    for k in 1:length(g.edges)
        e = g.edges[k]
        u = e.first
        v = e.second
        i = min(color[u], color[v])
        j = max(color[u], color[v])
        if twocolorindex[i, j] == 0
            seen_twocolors += 1
            twocolorindex[i, j] = seen_twocolors
            push!(edge_count, 0)
        end
        idx = twocolorindex[i, j]
        edge_count[idx] += 1
    end
    # edges sorted by twocolor subgraph
    sorted_edges = Array{Vector{MyPair{Int}}}(undef, seen_twocolors)
    for idx in 1:seen_twocolors
        sorted_edges[idx] = MyPair{Int}[]
        sizehint!(sorted_edges[idx], edge_count[idx])
    end

    for i in 1:length(g.edges)
        e = g.edges[i]
        u = e.first
        v = e.second
        i = min(color[u], color[v])
        j = max(color[u], color[v])
        idx = twocolorindex[i, j]
        push!(sorted_edges[idx], MyPair(u, v))
    end

    # list of unique vertices in each twocolor subgraph
    vertexmap = Array{Vector{Int}}(undef, seen_twocolors)

    postorder = Array{Vector{Int}}(undef, seen_twocolors)
    parents = Array{Vector{Int}}(undef, seen_twocolors)

    # temporary lookup map from global index to subgraph index
    revmap = zeros(Int, num_vertices(g))

    adjcount = zeros(Int, num_vertices(g))

    cmap = zeros(Int, 0) # shared storage for DFS
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
                push!(vlist, u)
                revmap[u] = length(vlist)
                adjcount[u] = 0
            end
            if revmap[v] == 0
                push!(vlist, v)
                revmap[v] = length(vlist)
                adjcount[v] = 0
            end
            adjcount[u] += 1
            adjcount[v] += 1
        end

        # set up offsets for adjlist
        offset = Array{Int}(undef, length(vlist) + 1)
        offset[1] = 1
        for k in 1:length(vlist)
            offset[k+1] = offset[k] + adjcount[vlist[k]]
            adjcount[vlist[k]] = 0
        end
        # adjlist for node u in twocolor idx starts at
        # vec[offset[u]]
        # u has global index vlist[u]
        vec = Array{Int}(undef, offset[length(vlist)+1] - 1)

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

        resize!(cmap, length(vlist))
        order, parent =
            reverse_topological_sort_by_dfs(vec, offset, length(vlist), cmap)

        for k in 1:length(vlist)
            # clear for reuse
            revmap[vlist[k]] = 0
        end

        postorder[idx] = order
        parents[idx] = parent
        vertexmap[idx] = vlist
    end

    return RecoveryInfo(
        vertexmap,
        postorder,
        parents,
        color,
        num_colors,
        num_edges(g),
        local_indices,
    )
end

function indirect_recover_structure(rinfo::RecoveryInfo)
    N = length(rinfo.color)
    nnz = rinfo.nnz

    I = zeros(Int, nnz + N)
    J = zeros(Int, nnz + N)

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
            i, j = normalize(i, j)
            k += 1
            I[k] = i
            J[k] = j
        end
    end

    @assert k == nnz + N

    return I, J
end

# edgelist is nonzeros in hessian, *including* nonzeros on the diagonal
function hessian_color_preprocess(
    edgelist,
    num_total_var,
    seen_idx = IndexedSet(0),
)
    resize!(seen_idx, num_total_var)
    I = Int[]
    J = Int[]
    for (i, j) in edgelist
        push!(seen_idx, i)
        push!(seen_idx, j)
        push!(I, i)
        push!(J, j)
    end
    local_indices = sort!(seen_idx.nzidx[1:seen_idx.nnz])
    empty!(seen_idx)

    global_to_local_idx = seen_idx.nzidx # steal for storage
    for k in 1:length(local_indices)
        global_to_local_idx[local_indices[k]] = k
    end
    # only do the coloring on the local indices
    for k in 1:length(I)
        I[k] = global_to_local_idx[I[k]]
        J[k] = global_to_local_idx[J[k]]
    end

    g = gen_adjlist(I, J, length(local_indices))

    color, num_colors = acyclic_coloring(g)

    @assert length(color) == num_vertices(g)

    rinfo = recovery_preprocess(g, color, num_colors, local_indices)

    I, J = indirect_recover_structure(rinfo)
    #nnz = num_edges(g)
    # convert back to global indices
    for k in 1:length(I)
        I[k] = local_indices[I[k]]
        J[k] = local_indices[J[k]]
    end

    return I, J, rinfo
end

# allocate a seed matrix
function seed_matrix(rinfo::RecoveryInfo)
    return Array{Float64}(undef, length(rinfo.local_indices), rinfo.num_colors)
end

function prepare_seed_matrix!(R, rinfo::RecoveryInfo)
    N = length(rinfo.color) # number of local variables
    @assert N == length(rinfo.local_indices)
    @assert size(R, 1) == N
    @assert size(R, 2) == rinfo.num_colors
    fill!(R, 0.0)
    for i in 1:N
        R[i, rinfo.color[i]] = 1
    end
    return
end

# recover the hessian values from the hessian-matrix solution
# stored_values is a temporary vector of length >= length(rinfo.local_indices)
function recover_from_matmat!(V, R, rinfo::RecoveryInfo, stored_values)
    N = length(rinfo.color) # number of local variables
    nnz = rinfo.nnz

    @assert length(V) == nnz + N
    @assert length(stored_values) >= length(rinfo.local_indices)

    # diagonal entries
    k = 0

    for i in 1:N
        k += 1
        V[k] = R[i, rinfo.color[i]]
    end

    for t in 1:length(rinfo.vertexmap)
        vmap = rinfo.vertexmap[t]
        order = rinfo.postorder[t]
        parent = rinfo.parents[t]

        for z in 1:length(order)
            @inbounds stored_values[z] = 0.0
        end

        @inbounds for z in 1:length(order)
            v = order[z]
            p = parent[v]
            (p == 0) && continue

            i = vmap[v]
            j = vmap[p]

            value = R[i, rinfo.color[j]] - stored_values[v]
            stored_values[p] += value

            k += 1
            V[k] = value
        end
    end

    @assert k == nnz + N

    return
end

end
