
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
    rules[eval(funsym)] = arg -> replace_x(exp,arg) 
end

curlyexpr(x::Expr) = isexpr(x.args[2],:parameters) ? x.args[3] : x.args[2]
function gencurlyloop(x::Expr, code; escape=false)
    idxstart = 3
    if isexpr(x.args[2], :parameters)
        cond = escape ? esc(x.args[2].args[1]) : x.args[2].args[1]
        code = :(
        if $cond
            $code
        end)
        idxstart = 4
    end

    for level in length(x.args):-1:idxstart
        if escape
            code = Expr(:for, esc(x.args[level]),code)
        else
            code = Expr(:for, x.args[level],code)
        end
    end
    return code
end



function quoteTree(x::Expr, datalist::Dict, iterstack)
    if isexpr(x, :call)
        # evaluate, don't quote the function name
        code = :(Expr(:call,$(x.args[1])))
        for y in x.args[2:end]
            push!(code.args,quoteTree(y, datalist, iterstack))
        end
        return code
    elseif isexpr(x, :curly)
        @assert x.args[1] == :sum || x.args[1] == :prod # special sum syntax
        code = :(Expr(:curly,$(quot(x.args[1]))))
        idxstart = 3
        if isexpr(x.args[2], :parameters)
            idxstart = 4
            if length(x.args[2].args) != 1
                error("No commas after semicolon allowed in $(x.args[1]) expression, use && for multiple conditions")
            end
        end
        for ex in x.args[idxstart:end]
            # iteration variables
            push!(iterstack, ex.args[1])
        end
        if idxstart == 4
            push!(code.args,:(Expr(:parameters,$(quoteTree(x.args[2].args[1],datalist,iterstack)))))
        end
        # body of expression
        push!(code.args,quoteTree(x.args[idxstart-1],datalist,iterstack))
        # store iteration sets
        for ex in x.args[end:-1:idxstart]
            pop!(iterstack)
            @assert isexpr(ex,:(=)) || isexpr(ex,:in)
            #itrset = symbol(string("iterset",string(gensym())))
            rhs = quoteTree(ex.args[2], datalist, iterstack)
            push!(code.args, :(Expr(:(=), $(quot(ex.args[1])), $rhs)))
        end
        return code
    elseif isexpr(x, :ref)
        # for symbolic expressions, leave them in the tree but collect the values separately.
        # if an array reference, just store the array object and not each element
        # NOTE: this assumes the values of the array will not change!
        var = x.args[1]
        datalist[var] = var
        for idxvar in x.args[2:end]
            isa(idxvar, Symbol) || continue
            # if this isn't a local index variable, we need to save it
            if !(idxvar in iterstack)
                datalist[idxvar] = idxvar
            end
        end

        return quot(x)
    else
        if !(isexpr(x, :(:)) || isexpr(x, :comparison) || isexpr(x, :&&) || isexpr(x, :||))
            error("Unrecognized expression $x")
        end
        code = :(Expr($(quot(x.head))))
        for y in x.args[1:end]
            push!(code.args,quoteTree(y, datalist, iterstack))
        end
        return code
    end
end

function quoteTree(x::Symbol, datalist, iterstack)
    if !(x in iterstack)
        datalist[x] = x
    end
    return quot(x)
end

quoteTree(x, datalist, iterstack) = x

addToVarList!(l,x::Placeholder) = push!(l,getindex(x))
addToVarList!(l,x) = nothing

function genVarList(x::Expr, arrname)
    if x.head == :curly
        code = genVarList(curlyexpr(x),arrname)
        return gencurlyloop(x,code, escape=true)
    elseif x.head == :call
        return Expr(:block,[genVarList(y,arrname) for y in x.args[2:end]]...)
    else
        return :(addToVarList!($arrname,$(esc(x))))
    end
end

genVarList(x::Number,arrname) = nothing
genVarList(x::Symbol,arrname) = :(addToVarList!($arrname,$(esc(x))))

type SymbolicOutput
    tree
    inputnames
    inputvals
    indexlist # indices of placeholders as they appear in the expression
              # useful when multiple expressions have the same structure
    mapfromcanonical
    maptocanonical
    hashval # for identifying expressions with identical trees
end

macro processNLExpr(x)
    indexlist = :(idxlist = Int[]; $(genVarList(x, :idxlist)); idxlist)
    datalist = Dict()
    tree = esc(quoteTree(x, datalist, {}))
    inputnames = Expr(:tuple)
    inputvals = Expr(:tuple)
    for (k,v) in datalist
        push!(inputnames.args, quot(k))
        push!(inputvals.args, esc(v))
    end
    return :(SymbolicOutput(genExprGraph($tree), $inputnames, $inputvals,$indexlist, $(hash(x))))
end

function SymbolicOutput(tree, inputnames, inputvals, indexlist, hashval)
    # compute canonical indices and maps
    unq = unique(indexlist)
    nidx = length(unq)
    # canonical is 1:nidx
    maptocanonical = Dict{Int,Int}()
    for k in 1:nidx
        maptocanonical[unq[k]] = k
    end
    return SymbolicOutput(tree, inputnames, inputvals, indexlist, unq, maptocanonical, hashval)
end

export @processNLExpr

# turn each node in the expression tree into an ExprNode
# this expression is kth argument in parent expression
genExprGraph(x) = genExprGraph(x, nothing, nothing)

function genExprGraph(x::Expr, parent, k)
    parentarr = parent === nothing ? [] : [(parent,k)]
    if isexpr(x, :call)
        thisnode = ExprNode(x, parentarr, nothing, nothing)
        for i in 2:length(x.args)
            x.args[i] = genExprGraph(x.args[i], thisnode, i)
        end
        return thisnode
    elseif isexpr(x, :curly)
        @assert (x.args[1] == :sum) || (x.args[1] == :prod)
        thisnode = ExprNode(x, parentarr, nothing, nothing)
        if isexpr(x.args[2], :parameters) # filter conditions
            x.args[3] = genExprGraph(x.args[3], thisnode, nothing)
        else
            x.args[2] = genExprGraph(x.args[2], thisnode, nothing)
        end
        return thisnode
    else
        @assert isexpr(x, :ref)
        return ExprNode(x, parentarr, nothing, nothing)
    end
end

genExprGraph{T<:Number}(x::T, parent, k) = x
genExprGraph(x, parent, k) = ExprNode(x, [(parent,k)], nothing, nothing)

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
        oper = x.ex.args[1]
        @assert oper == :sum || oper == :prod
        # compute value of this node, need to use a loop
        if oper == :sum
            push!(expr_out.args, :( $(x.value) = zero(T) ))
        else # :prod
            push!(expr_out.args, :( $(x.value) = one(T) ))
        end
        code = quote end

        valexpr = forwardpass(curlyexpr(x.ex), code)

        if oper == :sum
            code = :( $code; $(x.value) += $valexpr )
        else # :prod
            code = :( $code; $(x.value) *= $valexpr )
        end

        push!(expr_out.args, gencurlyloop(x.ex, code))
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
forwardvalue(x, placevalues, placeindex_in) = float(x)

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
        if f == (+) || (isexpr(p.ex,:curly) && f == :sum)
            push!(expr_out.args, :( $(x.deriv) += $(p.deriv) ))
        elseif isexpr(p.ex,:curly) && f == :prod
            # potentially numerically unstable if x.value ~= 0
            push!(expr_out.args, :( $(x.deriv) += $(p.deriv)*$(p.value)/$(x.value) ))
        elseif f == (-)
            if length(p.ex.args) > 2 && k == 2
                push!(expr_out.args, :( $(x.deriv) += $(p.deriv) ))
            else
                push!(expr_out.args, :( $(x.deriv) -= $(p.deriv) ))
            end
        elseif f == (*)
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
        elseif f == (^)
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
        elseif f == (/)
            if k == 2 # numerator
                denom = getvalue(p.ex.args[3])
                push!(expr_out.args,
                  :( $(x.deriv) += $(p.deriv)/$denom ))
            else
                @assert k == 3 # denominator
                numer = getvalue(p.ex.args[2])
                push!(expr_out.args,
                  :( $(x.deriv) += -1*$(p.deriv)*$numer*($(x.value))^(-2) ))
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
        @assert x.ex.args[1] == :sum || x.ex.args[1] == :prod
        exprbody = curlyexpr(x.ex)
        if !isa(exprbody,ExprNode) # expression inside sum doesn't depend on input
            return
        end
        # need to do a mini forward pass here for each child, because we reuse the nodes
        cleargraph(exprbody)
        code = quote end
        forwardpass(exprbody, code)
        revpass(exprbody, code)
 
        push!(expr_out.args, gencurlyloop(x.ex,code))
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
        function $(fname){T}(__placevalues::Vector{T}, __placeindex_in, __output, __placeindex_out)
            $out
            return $fval
        end
    end

    return fexpr

end

# gradient evaluation parametric on "inputvals"
function genfgrad_parametric(x::SymbolicOutput)
    out = Expr(:block)
    # load data into local scope
    for i in 1:length(x.inputnames)
        push!( out.args, :( $(x.inputnames[i]) = __inputvals[$i] ))
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
        function $(fname){T}(__placevalues::Vector{T}, __placeindex_in, __output, __placeindex_out, __inputvals)
            $out
            return $fval
        end
    end

    return eval(fexpr)

end

# expression evaluation parametric on "inputvals"
# forward pass only
function genfval_parametric(x::SymbolicOutput)
    out = Expr(:block)
    # load data into local scope
    for i in 1:length(x.inputnames)
        push!( out.args, :( $(x.inputnames[i]) = __inputvals[$i] ))
    end
    fval = forwardpass(x.tree, out)
    fname = gensym()
    fexpr = quote
        function $(fname){T}(__placevalues::Vector{T}, __placeindex_in, __inputvals)
            $out
            return $fval
        end
    end

    return eval(fexpr)
end

function genfval_simple(x::SymbolicOutput)
    fexpr = genfval_parametric(x)
    f = eval(fexpr)
    return xvals -> f(xvals, IdentityArray(), x.inputvals)
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

export genfgrad, genfgrad_simple, genfgrad_parametric, genfval_parametric, genfval_simple
