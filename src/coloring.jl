using Graphs
import DataStructures


function gen_adjlist(IJ,nel)
    g = simple_graph(nel, is_directed=false)
    for (i,j) in IJ
        i == j && continue
        add_edge!(g, i, j)
    end
    return g

end

export gen_adjlist

# acyclic coloring algorithm of Gebremdehin, Tarafdar, Manne, and Pothen
# "New Acyclic and Star Coloring Algorithms with Application to Computing Hessians"
# SIAM J. Sci. Comput. 2007
function acyclic_coloring(g)
    @assert !is_directed(g)
    
    num_colors = 0
    forbiddenColors = Int[]
    firstNeighbor = Array((Int,Int),0)
    firstVisitToTree = [normalize(source(e),target(e)) => (0,0) for e in edges(g)]
    color = fill(0, num_vertices(g))
    colored(i) = (color[i] != 0)
    # disjoint set forest of edges in the graph
    S = DataStructures.DisjointSets{(Int,Int)}(collect(keys(firstVisitToTree)))
    reverse_intmap = Array((Int,Int),length(S.intmap))
    for (k,v) in S.intmap
        reverse_intmap[v] = k
    end

    function prevent_cycle(v,w,x)
        er = DataStructures.find_root(S, normalize(w,x))
        e = reverse_intmap[er]
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
        e1 = DataStructures.find_root(S, normalize(v,w))
        e2 = DataStructures.find_root(S, normalize(w,x))
        if e1 != e2
            union!(S, normalize(v,w), normalize(w,x))
        end
    end

    for v in 1:num_vertices(g)
        for eg in out_edges(v,g)
            w = target(eg)
            colored(w) || continue
            forbiddenColors[color[w]] = v
        end
        for eg in out_edges(v,g)
            w = target(eg)
            colored(w) || continue
            for eg2 in out_edges(w,g)
                x = target(eg2)
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

        for eg in out_edges(v,g)
            w = target(eg)
            colored(w) || continue
            grow_star(v,w)
        end
        for eg in out_edges(v,g)
            w = target(eg)
            colored(w) || continue
            for eg2 in out_edges(w,g)
                x = target(eg2)
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
    return Set([color[i], color[j]])
end

immutable RecoveryInfo
    twocolorgraphs::Vector{SimpleGraph}
    vertexmap::Vector{Vector{Int}}
    postorder::Vector{Vector{Int}}
    parents::Vector{Vector{Int}}
    color::Vector{Int}
end

function recovery_preprocess(g,color; verify_acyclic::Bool=false)
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
        if verify_acyclic
            @assert !test_cyclic_by_dfs(s)
        end

        push!(twocolorgraphs, s)
        push!(vertexmap, vmap)
    end

    # list the vertices in postorder
    postorder = [reverse(topological_sort_by_dfs(s)) for s in twocolorgraphs]
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
        s = rinfo.twocolorgraphs[t]
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
        s = rinfo.twocolorgraphs[t]
        vmap = rinfo.vertexmap[t]
        order = rinfo.postorder[t]
        parent = rinfo.parents[t]
        stored_values[1:num_vertices(s)] = 0.0

        for z in 1:num_vertices(s)
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

function gen_hessian_sparse_color_parametric(s::SymbolicOutput, num_total_vars, hessian_matmat!, hessian_IJ)
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
    # vectors for hessian-vec product
    dualvec = Array(Dual{Float64}, num_total_vars)
    dualout = Array(Dual{Float64}, num_total_vars)

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
