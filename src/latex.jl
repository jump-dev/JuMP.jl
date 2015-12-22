#=
GenericRangeConstraint{SymbolicOutput} 
  terms: SymbolicOutput 
    tree: ExprNode 
      ex: Expr 
        head: Symbol call
        args: Array(Any,(3,))
          1: Symbol -
          2: ExprNode 
          3: Int64 5
        typ: Any
      parents: Array(None,(0,)) None[]
      value: Nothing nothing
      deriv: Nothing nothing
    inputnames: (Symbol,) (:x,)
    inputvals: (JuMPDict##4319{Variable},) (x[i], for all i in {1..5} free,)
    indexlist: Array(Int64,(5,)) [1,2,3,4,5]
    mapfromcanonical: Array(Int64,(5,)) [1,2,3,4,5]
    maptocanonical: Dict{Int64,Int64} len 5
      5: Int64 5
      4: Int64 4
      2: Int64 2
      3: Int64 3
      1: Int64 1
    hashval: Uint64 12263162890803976630
  lb: Float64 -Inf
  ub: Float64 0.0
=#

function print_nl(nlcon::GenericRangeConstraint{ReverseDiffSparse.SymbolicOutput})
    tree = nlcon.terms.tree

    recurse(tree::Any) = string(tree)

    function recurse(tree::ReverseDiffSparse.ExprNode)
        #println("RECURSE RDS.EN")
        #dump(tree)
        ex = tree.ex
        if ex.head == :call
            if ex.args[1] == :(-) || ex.args[1] == :(+)
                return string(recurse(ex.args[2]), " $(ex.args[1]) ", recurse(ex.args[3]))
            end
        elseif ex.head == :curly
            if ex.args[1] == :sum
                # sum{}
                # 1 :sum
                # 2 :ref (ExprNode)
                # 3 :range1
                # 4 :range2...
                num_idx = length(ex.args) - 2

                # First part: what goes in the sum
                in_sum_str = recurse(ex.args[2])
                if ex.args[2].ex.head != :ref
                    in_sum_str = "("*in_sum_str*")"
                end

                # Second part: sum over
                sums = {}
                for term = 1:num_idx
                    sum_over = ex.args[2+term]::Expr
                    sum_idx       = sum_over.args[1]
                    sum_idx_start = sum_over.args[2].args[1]
                    sum_idx_end   = sum_over.args[2].args[2]
                    push!(sums, string("\\sum_{", string(sum_idx), "=", sum_idx_start,
                                        "}^{", sum_idx_end, "}"))
                end
                return join(reverse(sums)," ") * " " * in_sum_str
            end
        elseif ex.head == :ref
            sum_var  = ex.args[1]
            sum_idxs = ex.args[2:end]
            in_sum_str = string(sum_var) * "_{" *
                            join(map(string,sum_idxs),",") * "}"
            return in_sum_str
        end
    end

    println(recurse(tree))
end