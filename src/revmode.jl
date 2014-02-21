
# generate code that can be compiled to compute the gradient

# load up derivative rules
remove_xp(x::Expr) = Expr(x.head, [remove_xp(ex) for ex in x.args]...)
remove_xp(s::Symbol) = (s == :xp) ? 1 : s
remove_xp(x) = x

replace_x(x::Expr, xvar) = Expr(x.head, [replace_x(ex,xvar) for ex in x.args]...)
replace_x(s::Symbol, xvar) = (s == :x) ? xvar : s
replace_x(x, xvar) = x

const rules = Dict()
for (funsym, exp) in Calculus.derivative_rules
    # remove xp from expression -- Calculus uses it for the chain rule
    exp = remove_xp(exp)
    # return a function that returns an expression with the given expression replacing the symbol x
    rules[funsym] = arg -> replace_x(exp,arg) 
end

#isinput(x) = isa(x, Placeholder)

function quoteTree(x::Expr, datalist::Dict, iterstack)
    if isexpr(x, :call)
        quoted = quot(x.args[1]) # this leaves the function names as symbols, instead of resolving them
        code = :(Expr(:call,$quoted))
        for y in x.args[2:end]
            push!(code.args,quoteTree(y, datalist, iterstack))
        end
        return Expr(:tuple,code,nothing)
    elseif isexpr(x, :curly)
        @assert x.args[1] == :sum # special sum syntax
        code = :(Expr(:curly,:sum))
        for ex in x.args[3:end]
            var,set = ex.args[1:2]
            push!(iterstack, :($var = first($set)))
        end
        push!(code.args,quoteTree(x.args[2],datalist,iterstack))
        # store iteration sets
        for ex in x.args[3:end]
            @assert isexpr(ex,:(=)) || isexpr(ex,:in)
            itrset = symbol(string("iterset",string(gensym())))
            datalist[itrset] = ex.args[2]
            ex.args[2] = itrset
            push!(code.args, quot(ex))
            pop!(iterstack)
        end
        return Expr(:tuple,code,nothing)
    else
        if !isexpr(x, :ref)
            error("Unrecognized expression $x")
        end
        
        # for symbolic expressions, leave them in the tree but collect the values separately.
        # if an array reference, just store the array object and not each element
        # NOTE: this assumes the values of the array will not change!
        var = x.args[1]
        datalist[var] = var

        return Expr(:tuple, quot(x),Expr(:block, iterstack..., :(isa($x,ReverseDiffSparse.Placeholder))))
    end
end

function quoteTree(x::Symbol, datalist, iterstack)
    datalist[x] = x
    return Expr(:tuple,quot(x),Expr(:block, iterstack..., :(isa($x,ReverseDiffSparse.Placeholder))))
end

quoteTree(x, datalist, iterstack) = Expr(:tuple,x,:(isa($x,ReverseDiffSparse.Placeholder)))

type SymbolicOutput
    tree
    inputnames
    inputvals
end

macro processNLExpr(x)
    datalist = Dict()
    iterstack = {}
    tree = esc(quoteTree(x, datalist, iterstack))
    inputnames = Expr(:tuple)
    inputvals = Expr(:tuple)
    for (k,v) in datalist
        push!(inputnames.args, quot(k))
        push!(inputvals.args, esc(v))
    end
    return :(SymbolicOutput(genExprGraph(inferInput($tree)), $inputnames, $inputvals))
end

export @processNLExpr

# turn each node in the expression tree into an ExprNode
# this expression is kth argument in parent expression
genExprGraph(x::(Expr,Any)) = genExprGraph(x, nothing, nothing)

function genExprGraph(t::(Expr,Any), parent, k)
    x,input = t
    if !input
        return x # collapse expressions that don't depend on the input
    end
    parentarr = parent === nothing ? [] : [(parent,k)]
    if isexpr(x, :call)
        thisnode = ExprNode(x, parentarr, nothing, nothing)
        for i in 2:length(x.args)
            x.args[i] = genExprGraph(x.args[i], thisnode, i)
        end
        return thisnode
    elseif isexpr(x, :curly)
        @assert x.args[1] == :sum
        thisnode = ExprNode(x, parentarr, nothing, nothing)
        x.args[2] = genExprGraph(x.args[2], thisnode, nothing)
        return thisnode
    else
        @assert isexpr(x, :ref)
        return ExprNode(x, parentarr, nothing, nothing)
    end
end

genExprGraph{T<:Number}(x::(T,Any), parent, k) = x[1]
genExprGraph(t, parent, k) = t[2] ? ExprNode(t[1], [(parent,k)], nothing, nothing) : t[1]

function inferInput(t::(Expr,Any))
    x,input = t
    if input != nothing
        return t
    end
    if isexpr(x, :call)
        inp = false
        for i in 2:length(x.args)
            x.args[i] = inferInput(x.args[i])
            inp |= x.args[i][2]
        end
        return (x,inp)
    elseif isexpr(x, :curly)
        x.args[2] = inferInput(x.args[2])
        return (x, true) # don't optimize out sum{}, since we need to evaluate it manually
    else
        error("Unexpected")
    end
end

inferInput{T<:Number}(x::(T,Any)) = (x[1],false)
inferInput(x) = x


function forwardpass(x::ExprNode, expr_out)
    @assert isexpr(expr_out, :block)

    # compute the value of each expression node by DFS
    # returns the variable name containing the result

    x.value = gensym() # generate a variable to represent the value of this node
    if isexpr(x.ex, :call)
        values = {}
        for i in 2:length(x.ex.args)
            push!(values, forwardpass(x.ex.args[i], expr_out))
        end
        fcall = Expr(:call, x.ex.args[1], values...)
        push!(expr_out.args, :( $(x.value) = $fcall ))
        return x.value
    elseif isexpr(x.ex, :curly)
        @assert x.ex.args[1] == :sum
        # compute value of this node, need to use a loop
        push!(expr_out.args, :( $(x.value) = zero(T) ))
        code = quote end
        valexpr = forwardpass(x.ex.args[2], code)
        code = :( $code; $(x.value) += $valexpr )
        for level in length(x.ex.args):-1:3
            code = Expr(:for, x.ex.args[level],code)
        end
        push!(expr_out.args, code)
        return x.value
    elseif isa(x.ex, Expr) || isa(x.ex, Symbol)
        # some other symbolic value, to evaluate at runtime
        push!(expr_out.args, :( $(x.value) = forwardvalue($(x.ex), __placevalues, __placeindex_in) ))
        return x.value
    else
        error("Shouldn't get here")
    end
end

forwardpass(x, expr_out) = :(forwardvalue($x, __placevalues, __placeindex_in))

forwardvalue(x::Placeholder, placevalues, placeindex_in) = placevalues[placeindex_in[getindex(x)]]
forwardvalue(x, placevalues, placeindex_in) = x

function revpass(x::ExprNode, expr_out)
    @assert isexpr(expr_out, :block)
    # compute the partial drivative wrt. each expression down the graph
    x.deriv = gensym()
    oneval = :( one(T) )
    zeroval = :( zero(T) )
    if length(x.parents) == 0
        push!(expr_out.args, :( $(x.deriv) = $oneval ) )
    else
        push!(expr_out.args, :( $(x.deriv) = $zeroval ) )
    end
    # if not all parents' derivatives are computed, come back later
    for (p, k) in x.parents
        if p.deriv === nothing
            return
        end
    end
    
    # x.deriv = sum_{p in parents} p.deriv*(\partial f_p / \partial x)
    for (p, k) in x.parents
        f = p.ex.args[1]
        if f == :(+) || (isexpr(p.ex,:curly) && f == :sum)
            push!(expr_out.args, :( $(x.deriv) += $(p.deriv) ))
        elseif f == :(-)
            if length(p.ex.args) > 2 && k == 2
                push!(expr_out.args, :( $(x.deriv) += $(p.deriv) ))
            else
                push!(expr_out.args, :( $(x.deriv) -= $(p.deriv) ))
            end
        elseif f == :(*)
            prd = gensym()
            push!(expr_out.args, :( $prd = $oneval ))
            for i in 2:length(p.ex.args)
                if i == k
                    continue
                else
                    push!(expr_out.args, :( $prd *= $(getvalue(p.ex.args[i])) ) )
                end
            end
            push!(expr_out.args, :( $(x.deriv) += $(p.deriv)*$prd ) )
        elseif f == :(^)
            if k == 2 # base
                exponent = getvalue(p.ex.args[3])
                push!(expr_out.args, 
                  :( $(x.deriv) += $(p.deriv)*$exponent*$(x.value)^($exponent-1) ))
            else
                @assert k == 3
                base = getvalue(p.ex.args[2])
                push!(expr_out.args,
                  :( $(x.deriv) += $(p.deriv)*$base^($(x.value))*log($base) ))
            end
        else
            # try one of the derivative rules
            haskey(rules, f) || error("Unrecognized function $f")
            k == 2 || error("Only one-argument version currently supported")
            fcall = rules[f](x.value)
            push!(expr_out.args, :( $(x.deriv) += $(p.deriv)*$fcall ) )
        end
    end

    # recurse through children
    if isexpr(x.ex,:call)
        for i in 2:length(x.ex.args)
            revpass(x.ex.args[i], expr_out)
        end
    elseif isexpr(x.ex,:curly)
        @assert x.ex.args[1] == :sum
        if !isa(x.ex.args[2],ExprNode) # expression inside sum doesn't depend on input
            return
        end
        # need to do a mini forward pass here for each child, because we reuse the nodes
        cleargraph(x.ex.args[2])
        code = quote end
        forwardpass(x.ex.args[2], code)
        revpass(x.ex.args[2], code)
        for level in length(x.ex.args):-1:3
            code = Expr(:for, x.ex.args[level],code)
        end
        push!(expr_out.args, code)
    else
        push!(expr_out.args, :( saverevvalue($(x.ex), $(x.deriv), __output, __placeindex_out) ))
    end
        
end

revpass(x, expr_out) = nothing

saverevvalue(x::Placeholder, val, output, placeindex_out) = (output[placeindex_out[getindex(x)]] += val)
saverevvalue(x, val, output, placeindex_out) = nothing

function genfgrad(x::SymbolicOutput)
    out = Expr(:block)
    # load data into local scope
    for i in 1:length(x.inputnames)
        push!( out.args, :( $(x.inputnames[i]) = $(x.inputvals[i]) ))
    end
    fval = forwardpass(x.tree, out)
    push!( out.args, :( beginreverse = nothing ) ) # for debugging, indicate start of reverse pass instructions
    revpass(x.tree,out)
    fname = gensym()
    # placeindex_in[i] specifies the index of the value of the ith
    # placeholder in placevalues.
    # placeindex_out[i] specifies the the index in which to output
    # the partial derivative wrt the ith placeholder
    fexpr = quote
        function $(fname){T}(__placevalues::Vector{T},__placeindex_in, __output, __placeindex_out)
            $out
            return $fval
        end
    end

    return fexpr

end

type IdentityArray
end

import Base.getindex
getindex(::IdentityArray,i) = i

function genfgrad_simple(x::SymbolicOutput)
    fexpr = genfgrad(x)
    f = eval(fexpr)
    return (xvals, out) -> (fill!(out, 0.0); f(xvals, IdentityArray(), out, IdentityArray()))
end

export genfgrad, genfgrad_simple
