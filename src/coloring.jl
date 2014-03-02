using Graphs
using DataStructures


function gen_adjlist(IJ)
    N = maximum([i for (i,j) in IJ])
    g = simple_adjlist(N, is_directed=false)
    for (i,j) in IJ
        i == j && continue
        add_edge!(g, i, j)
    end
    return g

end

export gen_adjlist

normalize(i,j) = (j > i) ? (j,i) : (i,j)
normalize(e) = normalize(e...)

# acyclic coloring algorithm of Gebremdehin, Tarafdar, Manne, and Pothen
# "New Acyclic and Star Coloring Algorithms with Application to Computing Hessians"
# SIAM J. Sci. Comput. 2007
function acyclic_coloring(IJ)
    IJ = collect(zip(IJ...))
    IJ_nodiag = collect(filter(ij -> ij[1] != ij[2], IJ))
    g = gen_adjlist(IJ)
    num_colors = 0
    forbiddenColors = Int[]
    firstNeighbor = Array((Int,Int),0)
    firstVisitToTree = [(i,j) => (0,0) for (i,j) in IJ_nodiag]
    color = fill(0, num_vertices(g))
    colored(i) = (color[i] != 0)
    # disjoint set forest of edges in the graph
    S = DisjointSets{(Int,Int)}(IJ_nodiag)

    function prevent_cycle(v,w,x)
        e = find_root(S, normalize(w,x))
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

    return color, S, num_colors, length(IJ_nodiag)+num_vertices(g)
end


type TwoColorTreeNode
    id::Int
    neighbors::Set{Int}
end

TwoColorTreeNode(id) = TwoColorTreeNode(id, Set{Int}())

type TwoColorTree
    nodes::Dict{Int,TwoColorTreeNode}
    leaves::Set{Int}
    color1::Int
    color2::Int
end

TwoColorTree(color1, color2) = TwoColorTree(Dict{Int,TwoColorTreeNode}(), Set{Int}(), color1, color2)

function new_node!(t::TwoColorTree, id, color)
    n = TwoColorTreeNode(id)
    @assert !haskey(t.nodes, id)
    @assert color == t.color1 || color == t.color2
    t.nodes[id] = n
end

import Graphs.add_edge!

function add_edge!(t::TwoColorTree, i, j)
    n1 = t.nodes[i]
    n2 = t.nodes[j]
    push!(n1.neighbors,j)
    push!(n2.neighbors,i)
end

function reverse_dict_lookup(d, v)
    for k in keys(d)
        if d[k] == v
            return k
        end
    end
    error("Not found")
end



function indirect_recover(hessian_matmat!, color, S, num_colors, x, inputvals, fromcanonical, tocanonical, V; structure=false)
    nnz = length(S)
    N = length(color)
    
    trees = Dict{Int,TwoColorTree}()
    
    # S has list of disjoint edges
    for e in keys(S.intmap)
        p = find_root(S, e)
        if !haskey(trees,p)
            # reverse lookup, this isn't ideal
            i,j = reverse_dict_lookup(S.intmap,p)
            color1 = color[i]
            color2 = color[j]
            @assert color1 != color2
            @assert color1 != 0 && color2 != 0
            trees[p] = TwoColorTree(color1, color2)
        end
        t = trees[p]
        i,j = e
        if !haskey(t.nodes, i)
            new_node!(t, i, color[i])
        end
        if !haskey(t.nodes, j)
            new_node!(t, j, color[j])
        end
        @assert color[i] != color[j]
        add_edge!(t, i, j)
    end

    # for each tree, find the leaves
    for t in values(trees)
        for n in values(t.nodes)
            @assert length(n.neighbors) != 0
            if length(n.neighbors) == 1 
                push!(t.leaves,n.id)
            end
        end
    end

    # generate vectors for hessian-vec product
    R = zeros(N,num_colors)
    for i in 1:N
        R[i,color[i]] = 1
    end

    if !structure
        hessian_matmat!(R,x, inputvals, fromcanonical, tocanonical)
    end
    
    # now, recover
    stored_values = zeros(N)
    if structure
        I = zeros(Int, nnz+N)
        J = zeros(Int, nnz+N)
    else
        @assert length(V) == nnz+N
    end
    
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

    for t in values(trees)
        for i in keys(t.nodes)
            stored_values[i] = 0
        end
        while !isempty(t.leaves)
            i = pop!(t.leaves)
            ni = t.nodes[i]
            @assert length(ni.neighbors) == 1
            j = first(ni.neighbors)

            k += 1
            if structure
                I[k] = i
                J[k] = j
            else
                V[k] = R[i,color[j]] - stored_values[i]
                stored_values[j] += V[k]
            end
            
            nj = t.nodes[j]
            delete!(nj.neighbors, i)
            delete!(t.nodes, i) # not necessary?
            if length(nj.neighbors) == 1
                # new leaf
                push!(t.leaves, j)
            elseif length(nj.neighbors) == 0
                delete!(t.leaves,j)
                delete!(t.nodes,j)
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
    color, S, num_colors = acyclic_coloring((I,J))
    
    I,J = indirect_recover(hessian_matmat!, color, S, num_colors, nothing, s.inputvals, s.mapfromcanonical, s.maptocanonical, nothing, structure=true)
    
    function eval_h(x,output_values, ex::SymbolicOutput)
        indirect_recover(hessian_matmat!, color, S, num_colors, x, ex.inputvals, ex.mapfromcanonical, ex.maptocanonical, output_values)
        # for now, we're inefficient
        #copy!(output_matrix.nzval, Hmat.nzval)
        #copy!(output_matrix.colptr, Hmat.colptr)
        #copy!(output_matrix.rowval, Hmat.rowval)
    end

    return I,J, eval_h


end

export gen_hessian_sparse_color_parametric

function to_H(s::SymbolicOutput, I, J, V, n)
    I2 = similar(I)
    J2 = similar(J)
    for k in 1:length(I)
        I2[k] = s.mapfromcanonical[I[k]]
        J2[k] = s.mapfromcanonical[J[k]]
    end
    return sparse(I2,J2,V, n, n)

end

export to_H





