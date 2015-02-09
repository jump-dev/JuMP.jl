
# generate code that can be compiled to compute the gradient

replace_x(x::Expr, xvar) = Expr(x.head, [replace_x(ex,xvar) for ex in x.args]...)
replace_x(s::Symbol, xvar) = (s == :x) ? xvar : s
replace_x(x, xvar) = x

const rules = Dict()
for (funsym, exp) in Calculus.symbolic_derivatives_1arg()
    # return a function that returns an expression with the given expression replacing the symbol x
    rules[funsym] = arg -> replace_x(exp,arg)
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



function quoteTree(x::Expr, datalist::Dict, iterstack, prefix, justrename::Bool=false)
    if isexpr(x, :call)
        # quote the function name
        code = justrename ? Expr(:call, x.args[1]) : :(Expr(:call,$(quot(x.args[1]))))
        for y in x.args[2:end]
            push!(code.args,quoteTree(y, datalist, iterstack, prefix, justrename))
        end
        return code
    elseif isexpr(x, :curly)
        @assert issum(x.args[1]) || isprod(x.args[1]) # special sum/prod syntax
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
            it = ex.args[1]
            if isa(it, Symbol) && initerstack(it, iterstack)
                error("$(it) cannot be reused as an iteration variable")
            elseif isexpr(it, :tuple) && any(x->initerstack(x,iterstack),it.args)
                error("$(it) contains duplicate iteration variables")
            end

            push!(iterstack, it)
        end
        if idxstart == 4
            push!(code.args,:(Expr(:parameters,$(quoteTree(x.args[2].args[1],datalist,iterstack, prefix, justrename)))))
        end
        # body of expression
        push!(code.args,quoteTree(x.args[idxstart-1],datalist,iterstack,prefix, justrename))
        # store iteration sets
        for ex in x.args[end:-1:idxstart]
            pop!(iterstack)
            @assert isexpr(ex,:(=)) || isexpr(ex,:in)
            #itrset = symbol(string("iterset",string(gensym())))
            rhs = quoteTree(ex.args[2], datalist, iterstack, prefix, justrename)
            push!(code.args, :(Expr(:(=), $(quot(ex.args[1])), $rhs)))
        end
        return code
    elseif isexpr(x, :ref)
        # for symbolic expressions, leave them in the tree but collect the values separately.
        # if an array reference, just store the array object and not each element
        # NOTE: this assumes the values of the array will not change!
        xold = x
        x = copy(x)
        for i in 1:length(x.args)
            x.args[i] = quoteTree(x.args[i], datalist, iterstack, prefix, true)
        end
        if length(iterstack) == 0
            return :(ifelse($xold, $(quot(x))))
        else
            return quot(x)
        end
    elseif isexpr(x, :quote)
        return x
    else
        if !(x.head in (:(:), :comparison, :&&, :||, :(.)))
            error("Unrecognized expression $x")
        end
        code = justrename ? Expr(x.head) : :(Expr($(quot(x.head))))
        for i in 1:length(x.args)
            if isexpr(x, :comparison) && iseven(i)
                push!(code.args, quot(x.args[i])) # don't rename comparison operators
            else
                push!(code.args,quoteTree(x.args[i], datalist, iterstack, prefix, justrename))
            end
        end
        return code
    end
end

# check if x is an iteration variable.
# note that multiple iteration variables may be
# present per iteration set, e.g., "(i,j) in S"
function initerstack(x::Symbol, iterstack)
    for v in iterstack
        if isa(v,Symbol)
            x == v && return true
        elseif isexpr(v,:tuple)
            (x in v.args) && return true
        else
            error("Unrecognized iteration variables $x")
        end
    end
    return false
end


function quoteTree(x::Symbol, datalist, iterstack, prefix, justrename)
    newsym = symbol("$prefix$x")
    if initerstack(x,iterstack)
        return justrename ? x : quot(x)
    else
        datalist[newsym] = x
    end
    quoted = justrename ? newsym : quot(newsym)
    #quoted = quot(newsym)
    if length(iterstack) == 0 && !justrename
        return :(ifelse($x, $quoted))
    else
        return quoted
    end
end

quoteTree(x, datalist, iterstack, prefix, justrename) = x

addToVarList!(l,inputvals,inputnames,hashsave,x::Placeholder) = push!(l,getplaceindex(x))
addToVarList!(l,inputvals,inputnames,hashsave,x) = nothing

function genVarList(x::Expr, arrname, inputvals, inputnames, hashsave)
    if x.head == :curly
        code = genVarList(curlyexpr(x),arrname, inputvals, inputnames, hashsave)
        return gencurlyloop(x,code, escape=true)
    elseif x.head == :call
        return Expr(:block,[genVarList(y,arrname, inputvals, inputnames, hashsave) for y in x.args[2:end]]...)
    elseif x.head in (:comparison, :&&, :||)
        return Expr(:block,[genVarList(y,arrname, inputvals, inputnames, hashsave) for y in x.args]...)
    else
        return :(addToVarList!($arrname,$inputvals,$inputnames,$hashsave,$(esc(x))))
    end
end

genVarList(x::Number,arrname,inputvals,inputnames,hashsave) = nothing
genVarList(x::Symbol,arrname,inputvals,inputnames,hashsave) = :(addToVarList!($arrname,$inputvals,$inputnames,$hashsave,$(esc(x))))

type SymbolicOutput
    tree
    inputnames
    inputvals
    indexlist # indices of placeholders as they appear in the expression
              # useful when multiple expressions have the same structure
    mapfromcanonical::Vector{Int}
    maptocanonical
    hashval # for identifying expressions with identical trees
    symbolprefix
end

function addToVarList!(l,inputvals,inputnames,hashsave,x::SymbolicOutput)
    append!(l,x.indexlist)
    push!(inputvals, x.inputvals...)
    push!(inputnames, x.inputnames...)
    hashsave[1] = hash(x.hashval, hashsave[1])
end

remove_prefix(x::Expr, prefix::String) = Expr(x.head, [remove_prefix(ex,prefix) for ex in x.args]...)
remove_prefix(s::Symbol, prefix::String) = beginswith(string(s),prefix) ? symbol(string(s)[length(prefix)+1:end]) : s
remove_prefix(x, prefix::String) = x
base_expression(s::SymbolicOutput) = remove_prefix(s.tree, s.symbolprefix)
export base_expression

# Type stable version of:
# isa(x,SymbolicOutput) ? x.tree : s
# We extend Base.ifelse as a hack for hygiene issues.
# See JuliaOpt/JuMP.jl#387.
Base.ifelse(x::SymbolicOutput,s::Union(Symbol,Expr)) = x.tree
Base.ifelse(::Any,s::Union(Symbol,Expr)) = s

macro processNLExpr(x)
    indexlist = genVarList(x, :idxlist, :inputvals, :inputnames, :hashsave)
    datalist = Dict()
    # in the corner case that x is just a symbol, promote it to an expr
    if isa(x, Symbol)
        x = Expr(:call,:+,x,0)
    end
    symbolprefix = string(gensym(),"##")
    tree = esc(quoteTree(x, datalist, Any[], symbolprefix))
    inputnames = Expr(:ref,:Symbol)
    inputvals = Expr(:ref,:Any)
    for (k,v) in datalist
        push!(inputnames.args, quot(k))
        push!(inputvals.args, esc(v))
    end
    #@show tree
    #@show indexlist
    return quote
        inputnames = $inputnames
        inputvals = $inputvals
        idxlist = Int[]
        hashsave = [$(hash(x))]
        $indexlist
        SymbolicOutput($tree, tuple(inputnames...), tuple(inputvals...), idxlist, hashsave[1], $symbolprefix)
    end
end

function SymbolicOutput(tree, inputnames, inputvals, indexlist, hashval, symbolprefix)
    # compute canonical indices and maps
    unq = unique(indexlist)
    nidx = length(unq)
    # canonical is 1:nidx
    maptocanonical = Dict{Int,Int}()
    for k in 1:nidx
        maptocanonical[unq[k]] = k
    end
    return SymbolicOutput(tree, inputnames, inputvals, indexlist, unq, maptocanonical, hashval, symbolprefix)
end

export @processNLExpr

# turn each node in the expression tree into an ExprNode
# this expression is kth argument in parent expression
genExprGraph(x) = genExprGraph(x, nothing, nothing)

function genExprGraph(x::Expr, parent, k)
    parentarr = parent === nothing ? [] : [(parent,k)]
    x = copy(x) # don't mutate original
    if isexpr(x, :call)
        thisnode = ExprNode(x, parentarr)
        for i in 2:length(x.args)
            x.args[i] = genExprGraph(x.args[i], thisnode, i)
        end
        return thisnode
    elseif isexpr(x, :comparison)
        thisnode = ExprNode(x, parentarr)
        for i in 1:2:length(x.args) # comparison symbols are in the middle
            x.args[i] = genExprGraph(x.args[i], thisnode, i)
        end
        return thisnode
    elseif isexpr(x, :&&) || isexpr(x, :||)
        thisnode = ExprNode(x, parentarr)
        for i in 1:length(x.args)
            x.args[i] = genExprGraph(x.args[i], thisnode, i)
        end
        return thisnode
    elseif isexpr(x, :curly)
        @assert issum(x.args[1]) || isprod(x.args[1])
        thisnode = ExprNode(x, parentarr)
        if isexpr(x.args[2], :parameters) # filter conditions
            x.args[3] = genExprGraph(x.args[3], thisnode, nothing)
        else
            x.args[2] = genExprGraph(x.args[2], thisnode, nothing)
        end
        return thisnode
    else
        @assert isexpr(x, :ref)
        return ExprNode(x, parentarr)
    end
end

genExprGraph{T<:Number}(x::T, parent, k) = x
genExprGraph(x, parent, k) = ExprNode(x, [(parent,k)], nothing, nothing)

# accumulator for prod{} terms
# we need to handle this specially so that we can compute
# derivatives when one of the terms is zero.
# (If two are zero, all derivatives are zero.)
immutable ProductAccumulator{T}
    allbutone::T # product of all but the smallest value
    smallest::T
end

ProductAccumulator{T}(::Type{T}) = ProductAccumulator(one(T),one(T))

function add_term{T}(p::ProductAccumulator{T},value::T)
    if abs(value) < abs(p.smallest)
        return ProductAccumulator(p.allbutone*p.smallest,value)
    else
        return ProductAccumulator(p.allbutone*value,p.smallest)
    end
end

asymbol(s::Symbol) = symbol(string(s,"accum"))


function forwardpass(x::ExprNode, expr_out)
    @assert isexpr(expr_out, :block)

    # compute the value of each expression node by DFS
    # returns the variable name containing the result

    x.value = gensym() # generate a variable to represent the value of this node
    if isexpr(x.ex, :call)
        values = Any[]
        for i in 2:length(x.ex.args)
            push!(values, forwardpass(x.ex.args[i], expr_out))
        end
        if x.ex.args[1] == :(^) # Use NaNMath.pow instead of ^
            fcall = Expr(:call, :pow, values...)
        else
            fcall = Expr(:call, x.ex.args[1], values...)
        end
        push!(expr_out.args, :( $(x.value) = $fcall ))
        return x.value
    elseif isexpr(x.ex, :curly)
        oper = x.ex.args[1]
        @assert issum(oper) || isprod(oper)
        # compute value of this node, need to use a loop
        if issum(oper)
            push!(expr_out.args, :( $(x.value) = zero(__T) ))
        else # :prod
            push!(expr_out.args, :( $(x.value) = one(__T) ))
            push!(expr_out.args, :( $(asymbol(x.value)) = ProductAccumulator(__T) ))
        end
        code = quote end

        valexpr = forwardpass(curlyexpr(x.ex), code)

        if issum(oper)
            code = :( $code; $(x.value) += $valexpr )
        else # :prod
            code = :( $code; $(x.value) *= $valexpr;
                      $(asymbol(x.value)) = add_term($(asymbol(x.value)),$valexpr))
        end

        push!(expr_out.args, gencurlyloop(x.ex, code))
        return x.value
    elseif isexpr(x.ex, :comparison)
        values = Any[]
        for i in 1:length(x.ex.args)
            if iseven(i)
                push!(values, x.ex.args[i]) # comparison operator
            else
                push!(values, forwardpass(x.ex.args[i], expr_out))
            end
        end
        fcall = Expr(:comparison, values...)
        push!(expr_out.args, :( $(x.value) = $fcall ))
        return x.value
    elseif isexpr(x.ex, :&&) || isexpr(x.ex, :||)
        values = Any[]
        for i in 1:length(x.ex.args)
            push!(values, forwardpass(x.ex.args[i], expr_out))
        end
        fcall = Expr(x.ex.head, values...)
        push!(expr_out.args, :( $(x.value) = $fcall ))
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

forwardvalue(x::Placeholder, placevalues, placeindex_in) = placevalues[placeindex_in[getplaceindex(x)]]
forwardvalue(x::SymbolicOutput, placevalues, placeindex_in) = error("Unexpected embedded expression $(x.tree).")
forwardvalue(x, placevalues, placeindex_in) = x

function revpass(x::ExprNode, expr_out)
    @assert isexpr(expr_out, :block)
    # compute the partial drivative wrt. each expression down the graph
    x.deriv = gensym()
    oneval = :( one(__T) )
    zeroval = :( zero(__T) )
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
        if f == :(+) || (isexpr(p.ex,:curly) && issum(f))
            push!(expr_out.args, :( $(x.deriv) += $(p.deriv) ))
        elseif isexpr(p.ex,:curly) && isprod(f)
            prodcode = quote
                if $(asymbol(p.value)).smallest == $(x.value)
                    $(x.deriv) += $(p.deriv)*$(asymbol(p.value)).allbutone
                else
                    $(x.deriv) += $(p.deriv)*$(p.value)/$(x.value)
                end
            end
            push!(expr_out.args, prodcode)
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
                  :( $(x.deriv) += $(p.deriv)*$exponent*pow($(x.value),$exponent-1) ))
            else
                @assert k == 3
                base = getvalue(p.ex.args[2])
                push!(expr_out.args,
                  :( $(x.deriv) += $(p.deriv)*pow($base,$(x.value))*log($base) ))
            end
        elseif f == :(/)
            if k == 2 # numerator
                denom = getvalue(p.ex.args[3])
                push!(expr_out.args,
                  :( $(x.deriv) += $(p.deriv)/$denom ))
            else
                @assert k == 3 # denominator
                numer = getvalue(p.ex.args[2])
                push!(expr_out.args,
                  :( $(x.deriv) += -1*$(p.deriv)*$numer*pow($(x.value),-2) ))
            end
        elseif f == :ifelse
            if k == 3 # true case
                push!(expr_out.args, :($(x.deriv) += ifelse($(getvalue(p.ex.args[2])),$(p.deriv), 0)))
            elseif k == 4
                push!(expr_out.args, :($(x.deriv) += ifelse($(getvalue(p.ex.args[2])),0,$(p.deriv))))
            else
                return # don't do any more work for the condition
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
        @assert issum(x.ex.args[1]) || isprod(x.ex.args[1])
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

saverevvalue(x::Placeholder, val, output, placeindex_out) = (output[placeindex_out[getplaceindex(x)]] += val)
saverevvalue(x, val, output, placeindex_out) = nothing

# gradient evaluation parametric on "inputvals"
function genfgrad_parametric(x::SymbolicOutput)
    out = Expr(:block)
    exgraph = genExprGraph(x.tree)
    fval = forwardpass(exgraph, out)
    push!( out.args, :( beginreverse = nothing ) ) # for debugging, indicate start of reverse pass instructions
    revpass(exgraph,out)
    # placeindex_in[i] specifies the index of the value of the ith
    # placeholder in placevalues.
    # placeindex_out[i] specifies the the index in which to output
    # the partial derivative wrt the ith placeholder
    fexpr = quote
        function _FGRAD_{__T}(__placevalues::Vector{__T}, __placeindex_in, __output, __placeindex_out)
            $out
            return $fval
        end
    end
    # add arguments for inputnames -- local data
    for i in 1:length(x.inputnames)
        push!(fexpr.args[2].args[1].args,x.inputnames[i])
    end

    return eval(:( local _FGRAD_; $fexpr; _FGRAD_))

end

# expression evaluation parametric on "inputvals"
# forward pass only
function genfval_parametric(x::SymbolicOutput)
    out = Expr(:block)
    fval = forwardpass(genExprGraph(x.tree), out)
    fexpr = quote
        function _FVAL_{__T}(__placevalues::Vector{__T}, __placeindex_in)
            $out
            return $fval
        end
    end
    # add arguments for inputnames -- local data
    for i in 1:length(x.inputnames)
        push!(fexpr.args[2].args[1].args,x.inputnames[i])
    end

    return eval(:( local _FVAL_; $fexpr; _FVAL_))
end

function genfval_simple(x::SymbolicOutput)
    f = genfval_parametric(x)
    return xvals -> f(xvals, IdentityArray(), x.inputvals...)
end

type IdentityArray
end

Base.getindex(::IdentityArray,i) = i

function genfgrad_simple(x::SymbolicOutput)
    fexpr = genfgrad_parametric(x)
    f = eval(fexpr)
    return (xvals, out) -> (fill!(out, 0.0); f(xvals, IdentityArray(), out, IdentityArray(),x.inputvals...))
end

export genfgrad, genfgrad_simple, genfgrad_parametric, genfval_parametric, genfval_simple
