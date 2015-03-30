using Graphs
import DataStructures

#Graphs.vertex_index(v::Int,g::Graphs.GenericGraph{Int64,Graphs.Edge{Int64},Array{Int64,1},Array{Graphs.Edge{Int64},1},Array{Array{Graphs.Edge{Int64},1},1}}) = v

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

function twocolorset_of_edge(e,g,color)
    i = source(e,g)
    j = target(e,g)
    return normalize_p(color[i],color[j])
end

immutable RecoveryInfo
    twocolorgraphs::Vector{SimpleGraph}
    vertexmap::Vector{Vector{Int}}
    postorder::Vector{Vector{Int}}
    parents::Vector{Vector{Int}}
    color::Vector{Int}
end

function recovery_preprocess(g,color; verify_acyclic::Bool=false)
    twocoloredges = Dict{MyPair{Int},Vector{MyPair{Int}}}()
    twocolorvertices = Dict{MyPair{Int},Set{Int}}()
    for e in edges(g)
        twocolor = twocolorset_of_edge(e,g,color)
        if !haskey(twocoloredges, twocolor)
            twocoloredges[twocolor] = Array(MyPair{Int},0)
            twocolorvertices[twocolor] = Set{Int}()
        end
        push!(twocoloredges[twocolor], normalize_p(source(e,g),target(e,g)))
        push!(twocolorvertices[twocolor], source(e,g))
        push!(twocolorvertices[twocolor], target(e,g))
    end

    twocolorgraphs = SimpleGraph[]
    # map from vertices in two-color subgraphs to original vertices
    vertexmap = Array(Vector{Int},0)
    postorder = Array(Vector{Int},0)
    revmap = zeros(Int,num_vertices(g))
    cmap = zeros(Int,0)
    sizehint!(vertexmap, length(twocoloredges))
    sizehint!(postorder, length(twocoloredges))
    for twocolor in keys(twocoloredges)
        edgeset = twocoloredges[twocolor]
        vertexset = twocolorvertices[twocolor]
        n = length(vertexset)
        s = simple_graph(n, is_directed=false)

        vmap = Int[]
        sizehint!(vmap, length(vertexset))
        for v in vertexset
            push!(vmap, v)
            revmap[v] = length(vmap)
        end

        for e in edgeset
            add_edge!(s, revmap[e.first], revmap[e.second])
        end
        if verify_acyclic
            @assert !test_cyclic_by_dfs(s)
        end

        push!(twocolorgraphs, s)
        push!(vertexmap, vmap)
        # list the vertices in postorder
        resize!(cmap,num_vertices(s))
        push!(postorder, reverse_topological_sort_by_dfs(s,cmap))
    end

    # identify each vertex's parent in the tree
    parent = zeros(0)
    seen = falses(0)
    parents = Array(Vector{Int},0)
    sizehint!(parents,length(twocolorgraphs))
    for i in 1:length(twocolorgraphs)
        s = twocolorgraphs[i]::SimpleGraph
        resize!(parent,num_vertices(s))
        fill!(parent,0)
        resize!(seen,num_vertices(s))
        fill!(seen,false)
        for k in 1:num_vertices(s)
            v = postorder[i][k]
            seen[v] = true
            # find the neighbor that we haven't seen
            notseen = 0
            numseen = 0
            for eg in out_edges(v,s)
                w = target(eg)
                if seen[w]
                    numseen += 1
                else
                    notseen = w
                end
            end
            if numseen == length(out_edges(v,s)) - 1
                parent[v] = notseen
            else
                (numseen == length(out_edges(v,s))) || error("Error processing tree, invalid ordering")
            end
        end
        push!(parents, parent)
    end


    return RecoveryInfo(twocolorgraphs, vertexmap, postorder, parents, color)

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

    for t in 1:length(rinfo.twocolorgraphs)
        s = rinfo.twocolorgraphs[t]::SimpleGraph
        vmap = rinfo.vertexmap[t]
        order = rinfo.postorder[t]
        parent = rinfo.parents[t]

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

    for t in 1:length(rinfo.twocolorgraphs)
        s = rinfo.twocolorgraphs[t]::SimpleGraph
        vmap = rinfo.vertexmap[t]
        order = rinfo.postorder[t]
        parent = rinfo.parents[t]
        stored_values[1:num_vertices(s)] = 0.0

        @inbounds for z in 1:num_vertices(s)
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
    
    rinfo = recovery_preprocess(g, color)

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


# Modified from Graphs.jl to avoid call to reverse()

type TopologicalSortVisitor{V} <: Graphs.AbstractGraphVisitor
    vertices::Vector{V}

    function TopologicalSortVisitor(n::Int)
        vs = Array(Int, 0)
        sizehint!(vs, n)
        new(vs)
    end
end


function Graphs.examine_neighbor!{V}(visitor::TopologicalSortVisitor{V}, u::V, v::V, vcolor::Int, ecolor::Int)
    if vcolor == 1 && ecolor == 0
        throw(ArgumentError("The input graph contains at least one loop."))
    end
end

function Graphs.close_vertex!{V}(visitor::TopologicalSortVisitor{V}, v::V)
    push!(visitor.vertices, v)
end

function reverse_topological_sort_by_dfs{V}(graph::Graphs.AbstractGraph{V}, cmap::Vector{Int})
    @graph_requires graph vertex_list incidence_list vertex_map

    @assert length(cmap) == num_vertices(graph)
    fill!(cmap,0)
    visitor = TopologicalSortVisitor{V}(num_vertices(graph))

    for s in vertices(graph)
        if cmap[vertex_index(s, graph)] == 0
            Graphs.traverse_graph(graph, DepthFirst(), s, visitor, vertexcolormap=cmap)
        end
    end

    visitor.vertices
end
