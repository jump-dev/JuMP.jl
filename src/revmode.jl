
# generate code that can be compiled to compute the gradient

replace_x(x, xvar) = replace_param(x, :x, xvar)
replace_param(x::Expr, param::Symbol, xvar) = Expr(x.head, [replace_param(ex,param,xvar) for ex in x.args]...)
replace_param(s::Symbol, param::Symbol, xvar) = (s == param) ? xvar : s
replace_param(x, param::Symbol, xvar) = x

const rules = Dict()
for (funsym, exp) in Calculus.symbolic_derivatives_1arg()
    # return a function that returns an expression with the given expression replacing the symbol x
    rules[funsym] = arg -> replace_x(exp,arg)
end
# special case for nonsmooth function, we return a subgradient
rules[:abs] = arg -> replace_x(:(x > 0 ? one(x) : -one(x)),arg)

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


# `a in b` is a comparison after JuliaLang/julia#13078
const in_is_compare = VERSION >= v"0.5.0-dev+901"

function quoteTree(x::Expr, datalist::Dict, iterstack, prefix, parameters = [], justrename::Bool=false)
    if isexpr(x, :call)
        # quote the function name
        code = justrename ? Expr(:call, x.args[1]) : :(Expr(:call,$(quot(x.args[1]))))
        for y in x.args[2:end]
            push!(code.args,quoteTree(y, datalist, iterstack, prefix, parameters, justrename))
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
            push!(code.args,:(Expr(:parameters,$(quoteTree(x.args[2].args[1],datalist,iterstack, prefix, parameters, justrename)))))
        end
        # body of expression
        push!(code.args,quoteTree(x.args[idxstart-1],datalist,iterstack,prefix,parameters,justrename))
        # store iteration sets
        for ex in x.args[end:-1:idxstart]
            ex = ex::Expr
            pop!(iterstack)
            arg1, arg2 = if isexpr(ex, :(=)) || (!in_is_compare &&
                                                 isexpr(ex, :in))
                ex.args[1], ex.args[2]
            elseif in_is_compare && isexpr(ex, :comparison)
                @assert length(ex.args) == 3 && ex.args[2] === :in
                ex.args[1], ex.args[3]
            end
            # itrset = symbol(string("iterset", string(gensym())))
            rhs = quoteTree(arg2, datalist, iterstack, prefix,
                            parameters, justrename)
            push!(code.args, :(Expr(:(=), $(quot(arg1)), $rhs)))
        end
        return code
    elseif isexpr(x, :ref)
        # for symbolic expressions, leave them in the tree but collect the values separately.
        # if an array reference, just store the array object and not each element
        # NOTE: this assumes the values of the array will not change!
        xold = x
        x = copy(x)
        for i in 1:length(x.args)
            x.args[i] = quoteTree(x.args[i], datalist, iterstack, prefix, parameters, true)
        end
        if justrename # x[y[i]]
            return x
        elseif isa(xold.args[1],Expr) # v[i].x[j]
            return quot(x)
        else
            return Expr(:call,:replaceif, esc(xold.args[1]),quot(x),[quot(s) for s in x.args[2:end]]...)
        end
        #return quot(x)
    elseif isexpr(x, :quote)
        return x
    else
        if !(x.head in (:(:), :comparison, :&&, :||, :(.), :tuple))
            error("Unrecognized expression $x")
        end
        code = justrename ? Expr(x.head) : :(Expr($(quot(x.head))))
        for i in 1:length(x.args)
            if isexpr(x, :comparison) && iseven(i)
                push!(code.args, quot(x.args[i])) # don't rename comparison operators
            elseif isexpr(x.args[i], :quote) && !justrename # for dot syntax
                push!(code.args, quot(x.args[i]))
            else
                push!(code.args,quoteTree(x.args[i], datalist, iterstack, prefix, parameters, justrename))
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


function quoteTree(x::Symbol, datalist, iterstack, prefix, parameters, justrename)
    newsym = symbol("$prefix$x")
    if initerstack(x,iterstack)
        return justrename ? x : quot(x)
    else
        datalist[newsym] = x
    end
    if justrename
        return newsym
    elseif x in parameters
        return quot(newsym)
    else
        return Expr(:call,:replaceif,esc(x),quot(newsym))
    end
end

quoteTree(x, datalist, iterstack, prefix, parameters, justrename) = x



type SymbolicOutput
    tree
    inputnames
    inputvals
    indexlist # indices of placeholders as they appear in the expression
              # useful when multiple expressions have the same structure
    mapfromcanonical::Vector{Int}
    maptocanonical::Dict{Int,Int}
    hashval # for identifying expressions with identical trees
end

expression_data(s::SymbolicOutput) = s.inputvals

remove_prefix(x::Expr, prefix::AbstractString) = Expr(x.head, [remove_prefix(ex,prefix) for ex in x.args]...)
function remove_prefix(s::Symbol, prefix::AbstractString)
    str = string(s)
    if startswith(str,prefix)
        return symbol(str[search(str,"##",6)[2]+1:end])
    else
        return s
    end
end
remove_prefix(x, prefix::AbstractString) = x
base_expression(s::SymbolicOutput) = remove_prefix(s.tree, "__R")
export base_expression



macro processNLExpr(x)
    datalist = Dict()
    # in the corner case that x is just a symbol, promote it to an expr
    if isa(x, Symbol)
        x = Expr(:call,:+,x,0)
    end
    symbolprefix = string("__R",gensym(),"##")
    tree = quoteTree(x, datalist, Any[], symbolprefix)
    inputnames = Expr(:ref,:Symbol)
    inputvals = Expr(:ref,:Any)
    for (k,v) in datalist
        push!(inputnames.args, quot(k))
        push!(inputvals.args, esc(v))
    end
    #@show tree
    subexprcode = quote end
    for v in values(datalist)
        push!(subexprcode.args, :(hashval = hashif($(esc(v)),hashval)))
        push!(subexprcode.args, :(appendnames_and_vals_if!(inputnames,inputvals,$(esc(v)))))
    end
    return quote
        inputnames = $inputnames
        inputvals = $inputvals
        hashval = $(hash(symbolprefix,hash(x)))
        $subexprcode
        SymbolicOutput($tree, tuple(inputnames...), tuple(inputvals...), hashval)
    end
end

function SymbolicOutput(tree, inputnames, inputvals, hashval)
    return SymbolicOutput(tree, inputnames, inputvals, nothing, Int[], Dict{Int,Int}(), hashval)
end



type ParametricExpression{N}
    tree
    inputnames
    inputvals
    parameters::NTuple{N,Symbol}
    hashval
end

immutable ParametricExpressionWithParams{N}
    p::ParametricExpression{N}
    params::NTuple{N,Any}
end

expression_data(p::ParametricExpression) = p.inputvals
expression_data(p::ParametricExpressionWithParams) = p.p.inputvals

replace_param(x::ParametricExpressionWithParams, param::Symbol, xvar) = ParametricExpressionWithParams(x.p, map(s-> replace_param(s,param,xvar), x.params))

function Base.getindex{N}(p::ParametricExpression{N},params...)
    length(params) == N || error("Incorrect number of parameters to expression $p")
    return ParametricExpressionWithParams(p,params)
end

macro parametricExpr(args...)
    params = args[1:end-1]
    x = args[end]
    datalist = Dict()
    # in the corner case that x is just a symbol, promote it to an expr
    if isa(x, Symbol)
        x = Expr(:call,:+,x,0)
    end
    symbolprefix = string("__R",gensym(),"##")
    tree = quoteTree(x, datalist, Any[], symbolprefix, params)
    inputnames = Expr(:ref,:Symbol)
    inputvals = Expr(:ref,:Any)
    for (k,v) in datalist
        v in params && continue
        push!(inputnames.args, quot(k))
        push!(inputvals.args, esc(v))
    end
    subexprcode = quote end
    for v in values(datalist)
        v in params && continue
        push!(subexprcode.args, :(hashval = hashif($(esc(v)),hashval)))
        push!(subexprcode.args, :(appendnames_and_vals_if!(inputnames,inputvals,$(esc(v)))))
    end
    paramtup = Expr(:tuple)
    for p in params
        @assert isa(p,Symbol)
        push!(paramtup.args, quot(symbol(string("$symbolprefix$p"))))
    end
    return quote
        inputnames = $inputnames
        inputvals = $inputvals
        hashval = $(hash(symbolprefix,hash(x)))
        $subexprcode
        ParametricExpression($tree, tuple(inputnames...), tuple(inputvals...), $paramtup, hashval)
    end
end

base_expression(s::ParametricExpression) = remove_prefix(s.tree, "__R")

export @parametricExpr, @processNLExpr

hashif(x::ParametricExpression,h) = hash(x,h)
hashif(x,h) = h
function appendnames_and_vals_if!(inputnames,inputvals,x::ParametricExpression)
    for i in 1:length(x.inputnames)
        name = x.inputnames[i]::Symbol
        idx = findfirst(inputnames, name)
        if idx != 0
            if inputvals[idx] === x.inputvals[i]
                continue
            else
                sname = remove_prefix(name, "__R")
                error("Value of $sname changed value unexpectedly")
            end
        else
            push!(inputnames, name)
            push!(inputvals, x.inputvals[i])
        end
    end
end
appendnames_and_vals_if!(inputnames,inputvals,x) = nothing

# Keep the parametric expression in the tree
function replaceif{N}(x::ParametricExpression{N},s::(@compat Union{Symbol,Expr}), args...)
    length(args) == N || error("Incorrect number of parameters for expression $x")
    return x[args...]
end
replaceif(::Any,s::(@compat Union{Symbol,Expr}), args...) = s

# turn each node in the expression tree into an ExprNode
# this expression is kth argument in parent expression
genExprGraph(x) = genExprGraph(x, nothing, nothing, true, Dict())

function genExprGraph(x::Expr, parent, k, linear_so_far::Bool, expr_map)
    parentarr = parent === nothing ? [] : [(parent,k)]
    xold = x
    x = copy(x) # don't mutate original
    if haskey(expr_map, x)
        ex = expr_map[x]
        push!(ex.parents, (parent,k))
        return ex
    end
    if isexpr(x, :call)
        thisnode = ExprNode(x, parentarr, linear_so_far)
        #TODO: also handle subtraction, but need to keep track of coefficients
        if !(x.args[1] == :(+) || x.args[1] == :ifelse)
            linear_so_far = false
        end
        for i in 2:length(x.args)
            x.args[i] = genExprGraph(x.args[i], thisnode, i, linear_so_far, expr_map)
        end
        expr_map[xold] = thisnode
        return thisnode
    elseif isexpr(x, :comparison)
        thisnode = ExprNode(x, parentarr, linear_so_far)
        for i in 1:2:length(x.args) # comparison symbols are in the middle
            x.args[i] = genExprGraph(x.args[i], thisnode, i, false, Dict())
        end
        return thisnode
    elseif isexpr(x, :&&) || isexpr(x, :||)
        thisnode = ExprNode(x, parentarr, linear_so_far)
        for i in 1:length(x.args)
            x.args[i] = genExprGraph(x.args[i], thisnode, i, false, Dict())
        end
        return thisnode
    elseif isexpr(x, :curly)
        @assert issum(x.args[1]) || isprod(x.args[1])
        thisnode = ExprNode(x, parentarr, linear_so_far)
        if isprod(x.args[1])
            linear_so_far = false
        else
            @assert issum(x.args[1])
        end
        # for now, throw away expr_map when scope changes
        if isexpr(x.args[2], :parameters) # filter conditions
            x.args[3] = genExprGraph(x.args[3], thisnode, nothing, linear_so_far, Dict())
        else
            x.args[2] = genExprGraph(x.args[2], thisnode, nothing, linear_so_far, Dict())
        end
        return thisnode
    else
        @assert isexpr(x, :ref) || isexpr(x, :.)
        thisnode = ExprNode(x, parentarr, linear_so_far)
        expr_map[xold] = thisnode
        return thisnode
    end
end

genExprGraph{T<:Number}(x::T, parent, k, linear_so_far, expr_map) = x
function genExprGraph(x, parent, k, linear_so_far, expr_map)
    if (parent === nothing)
        return ExprNode(x,[], linear_so_far)
    end
    if haskey(expr_map, x)
        ex = expr_map[x]
        push!(ex.parents, (parent,k))
        return expr_map[x]
    else
        thisnode = ExprNode(x, [(parent,k)], linear_so_far)
        expr_map[x] = thisnode
        return thisnode
    end
end

# tradeoff between compilation time and derivative evaluation speed
const SPLAT_THRESHOLD = 50

function genExprGraph(x::ParametricExpressionWithParams, parent, k, linear_so_far, expr_map)
    tree = x.p.tree
    depth = expr_complexity(tree)
    if depth < SPLAT_THRESHOLD # just splat
        for i in 1:length(x.params)
            tree = replace_param(tree, x.p.parameters[i], x.params[i])
        end
        genExprGraph(tree, parent, k, linear_so_far, expr_map)
    else
        invoke(genExprGraph, (Any,Any,Any,Any,Any), x, parent, k, linear_so_far,expr_map)
    end
end

# number of nodes in expression graph
expr_complexity(x::Expr) = sum([expr_complexity(ex) for ex in x.args])
expr_complexity(x::ParametricExpressionWithParams) = expr_complexity(x.p.tree)
expr_complexity(x) = 1

addToVarList!(l,x::Placeholder) = push!(l,getplaceindex(x))
addToVarList!(l,x) = nothing
function genVarList(x::Expr, arrname)
    if x.head == :curly
        code = genVarList(curlyexpr(x),arrname)
        return gencurlyloop(x,code, escape=false)
    elseif x.head == :call
        return Expr(:block,[genVarList(y,arrname) for y in x.args[2:end]]...)
    elseif x.head in (:comparison, :&&, :||)
        return Expr(:block,[genVarList(y,arrname) for y in x.args]...)
    else
        return :(addToVarList!($arrname,$x))
    end
end
function genVarList{N}(x::ParametricExpressionWithParams{N},arrname)
    # for now, just splat it here
    tree = x.p.tree
    for i in 1:N
        tree = replace_param(tree, x.p.parameters[i], x.params[i])
    end
    return genVarList(tree,arrname)
end
genVarList(x::Number,arrname) = nothing
genVarList(x::Symbol,arrname) = :(addToVarList!($arrname,$x))

function process_indexmap(s::SymbolicOutput, indexlist)
    # compute canonical indices and maps
    unq = unique(indexlist)
    nidx = length(unq)
    # canonical is 1:nidx
    maptocanonical = Dict{Int,Int}()
    for k in 1:nidx
        maptocanonical[unq[k]] = k
    end
    s.mapfromcanonical = unq
    s.maptocanonical = maptocanonical
    return
end

function prepare_indexmap(s::SymbolicOutput,indexset)
    length(s.maptocanonical) > 0 && return
    genindexlist_parametric(s)(indexset,s.inputvals...)
    process_indexmap(s,indexset)
end

function prepare_indexmap(s::SymbolicOutput)
    # TODO: deprecate this version
    indexset = Int[]
    genindexlist_parametric(s)(indexset,s.inputvals...)
    process_indexmap(s,indexset)
end

export prepare_indexmap

function genindexlist_parametric(x::SymbolicOutput)
    out = Expr(:block)
    fexpr = quote
        function _IDXLIST_(indexlist)
            empty!(indexlist)
            $(genVarList(x.tree, :indexlist))
            return indexlist
        end
    end
    # add arguments for inputnames -- local data
    for i in 1:length(x.inputnames)
        push!(fexpr.args[2].args[1].args,x.inputnames[i])
    end

    return eval(:( let; local _IDXLIST_; $fexpr; end;))
end

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

const exprval_cache = Dict()
const exprrev_cache = Dict()

# skip_linear_sums: Don't compute values of sums which contribute to the expression linearly.
#              This is used for fusing reverse and forward modes:
#              http://link.springer.com/chapter/10.1007/978-1-4613-0075-5_35
#              Always false unless you understand what's going on.
function forwardpass(x::ExprNode, expr_out, skip_linear_sums::Bool)
    @assert isexpr(expr_out, :block)

    # compute the value of each expression node by DFS
    # returns the variable name containing the result

    if isa(x.value,Symbol) # already computed
        return x.value
    end

    x.value = gensym() # generate a variable to represent the value of this node
    if isexpr(x.ex, :call)
        values = Any[]
        for i in 2:length(x.ex.args)
            push!(values, forwardpass(x.ex.args[i], expr_out, skip_linear_sums))
        end
        if x.ex.args[1] == :(^) # Use NaNMath.pow instead of ^
            if x.ex.args[3] == 2 # special processing for x^2
                fcall = Expr(:call, :*, values[1], values[1])
            else
                fcall = Expr(:call, :pow, values...)
            end
        else
            fcall = Expr(:call, x.ex.args[1], values...)
        end
        # x.value = GenSym(abs(rand(Int))) # TODO: use on 0.4
        push!(expr_out.args, :( $(x.value) = $fcall ))
        return x.value
    elseif isexpr(x.ex, :curly)
        oper = x.ex.args[1]
        @assert issum(oper) || isprod(oper)
        # compute value of this node, need to use a loop
        if issum(oper)
            push!(expr_out.args, :( $(x.value) = zero(__T) ))
            x.linear_so_far && skip_linear_sums && return x.value # early exit
        else # :prod
            push!(expr_out.args, :( $(x.value) = one(__T) ))
            push!(expr_out.args, :( $(asymbol(x.value)) = ProductAccumulator(__T) ))
        end
        code = quote end

        valexpr = forwardpass(curlyexpr(x.ex), code, skip_linear_sums)

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
                push!(values, forwardpass(x.ex.args[i], expr_out, skip_linear_sums))
            end
        end
        fcall = Expr(:comparison, values...)
        push!(expr_out.args, :( $(x.value) = $fcall ))
        return x.value
    elseif isexpr(x.ex, :&&) || isexpr(x.ex, :||)
        values = Any[]
        for i in 1:length(x.ex.args)
            push!(values, forwardpass(x.ex.args[i], expr_out, skip_linear_sums))
        end
        fcall = Expr(x.ex.head, values...)
        push!(expr_out.args, :( $(x.value) = $fcall ))
        return x.value
    elseif isa(x.ex, Expr) || isa(x.ex, Symbol)
        # some other symbolic value, to evaluate at runtime
        # x.value = GenSym(abs(rand(Int))) # TODO: use on 0.4
        push!(expr_out.args, :( $(x.value) = forwardvalue($(x.ex), __placevalues, __placeindex_in) ))
        return x.value
    elseif isa(x.ex,ParametricExpressionWithParams)
        if !haskey(exprval_cache, x.ex.p.hashval)
            exprval_cache[x.ex.p.hashval] = genexprval(x.ex.p)
        end
        exprval = exprval_cache[x.ex.p.hashval]
        fcall = :($(exprval)(__placevalues,__placeindex_in))
        push!(fcall.args,x.ex.params...)
        push!(fcall.args,x.ex.p.inputnames...)
        push!(expr_out.args, :( $(x.value) = $fcall::__T ))
        return x.value
    else
        error("Shouldn't get here")
    end
end

forwardpass(x, expr_out, skip_linear_sums) = :(forwardvalue($x, __placevalues, __placeindex_in))

forwardvalue(x::Placeholder, placevalues, placeindex_in) = placevalues[placeindex_in[getplaceindex(x)]]
forwardvalue(x::SymbolicOutput, placevalues, placeindex_in) = error("Unexpected embedded expression $(x.tree).")
forwardvalue(x, placevalues, placeindex_in) = x

function revpass(x::ExprNode, expr_out; rootval= :(one(__T)), linear_sums=false )
    @assert isexpr(expr_out, :block)

    # if not all parents' derivatives are computed, come back later
    for (p, k) in x.parents
        if p.deriv === nothing
            return
        end
    end

    # compute the partial drivative wrt. each expression down the graph
    if isa(x.deriv,Symbol)
        return # already done
    end
    x.deriv = gensym()
    zeroval = :( zero(__T) )
    if length(x.parents) == 0
        push!(expr_out.args, :( $(x.deriv) = $rootval ) )
    else
        push!(expr_out.args, :( $(x.deriv) = $zeroval ) )
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
            if length(p.ex.args) == 3 # only two multiplicands, special case
                other = (k == 2) ? 3 : 2
                push!(expr_out.args, :( $(x.deriv) += $(p.deriv)*$(getvalue(p.ex.args[other])) ))
            else
                prd = gensym()
                push!(expr_out.args, :( $prd = one(__T) ))
                for i in 2:length(p.ex.args)
                    if i == k
                        continue
                    else
                        push!(expr_out.args, :( $prd *= $(getvalue(p.ex.args[i])) ) )
                    end
                end
                push!(expr_out.args, :( $(x.deriv) += $(p.deriv)*$prd ) )
            end
        elseif f == :(^)
            if k == 2 # base
                if p.ex.args[3] == 2 # special processing for x^2
                    push!(expr_out.args,
                      :( $(x.deriv) += $(p.deriv)*2*$(x.value) ) )
                else
                    exponent = getvalue(p.ex.args[3])
                    push!(expr_out.args,
                      :( $(x.deriv) += $(p.deriv)*$exponent*pow($(x.value),$exponent-1) ))
                end
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
            revpass(x.ex.args[i], expr_out, linear_sums=linear_sums)
        end
    elseif isexpr(x.ex,:curly)
        @assert issum(x.ex.args[1]) || isprod(x.ex.args[1])
        exprbody = curlyexpr(x.ex)
        if !isa(exprbody,ExprNode) && !x.linear_so_far # expression inside sum doesn't depend on input and haven't computed it yet
            return
        end
        # need to do a mini forward pass here for each child, because we reuse the nodes
        cleargraph(exprbody)
        code = quote end
        sval = forwardpass(exprbody, code, true)
        # did we already compute this?
        # if not, fuse forward+reverse here
        if linear_sums != false && x.linear_so_far && issum(x.ex.args[1])
            push!(code.args, :($linear_sums += $sval))
        end
        revpass(exprbody, code, linear_sums=linear_sums)
 
        push!(expr_out.args, gencurlyloop(x.ex,code))
    elseif isa(x.ex,ParametricExpressionWithParams)
        if !haskey(exprrev_cache, x.ex.p.hashval)
            exprrev_cache[x.ex.p.hashval] = genexprrev(x.ex.p)
        end
        exprrev = exprrev_cache[x.ex.p.hashval]
        fcall = :($(exprrev)($(x.deriv),__placevalues,__placeindex_in,__output,__placeindex_out))
        push!(fcall.args,x.ex.params...)
        push!(fcall.args,x.ex.p.inputnames...)
        if VERSION > v"0.4.0"
            push!(expr_out.args, :( $fcall::Void ))
        else
            push!(expr_out.args, :( $fcall::Nothing ))
        end
    else
        push!(expr_out.args, :( saverevvalue($(x.ex), $(x.deriv), __output, __placeindex_out) ))
    end
        
end

revpass(x, expr_out; linear_sums=false) = nothing


saverevvalue(x::Placeholder, val, output, placeindex_out) = (output[placeindex_out[getplaceindex(x)]] += val)
saverevvalue(x::ParametricExpressionWithParams, val, output, placeindex_out) = error("Unexpected embedded expression $(x.p.tree) (during reverse mode)")
saverevvalue(x, val, output, placeindex_out) = nothing

# gradient evaluation parametric on "inputvals"
function genfgrad_parametric(x::SymbolicOutput)
    out = Expr(:block)
    exgraph = genExprGraph(x.tree)
    fval = forwardpass(exgraph, out, true) # note we skip linear sums here
    push!( out.args, :( beginreverse = nothing ) ) # for debugging, indicate start of reverse pass instructions
    revpass(exgraph,out, linear_sums = fval)
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
    #@show fexpr

    return eval(:( let; local _FGRAD_; $fexpr; end;))

end

# expression evaluation parametric on "inputvals"
# forward pass only
function genfval_parametric(x::SymbolicOutput)
    out = Expr(:block)
    fval = forwardpass(genExprGraph(x.tree), out, false)
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

    return eval(:( let; local _FVAL_; $fexpr; end;))
end

function genexprval{N}(x::ParametricExpression{N})
    out = Expr(:block)
    fval = forwardpass(genExprGraph(x.tree), out, false)
    fexpr = quote
        function _EXPRVAL_{__T}(__placevalues::Vector{__T}, __placeindex_in)
            $out
            return convert(__T,$fval)
        end
    end
    # add arguments for parameters
    for i in 1:N
        push!(fexpr.args[2].args[1].args,x.parameters[i])
    end
    # add arguments for inputnames -- local data
    for i in 1:length(x.inputnames)
        push!(fexpr.args[2].args[1].args,x.inputnames[i])
    end
    #@show fexpr

    return eval(:( let; local _EXPRVAL_; $fexpr; end;))
end

function genexprrev{N}(x::ParametricExpression{N})
    out = Expr(:block)
    exgraph = genExprGraph(x.tree)
    fval = forwardpass(exgraph, out, true)
    revpass(exgraph, out, rootval=:(__deriv), linear_sums = fval)
    fexpr = quote
        function _EXPRREV_{__T}(__deriv::__T,__placevalues::Vector{__T}, __placeindex_in, __output, __placeindex_out)
            $out
            nothing
        end
    end
    # add arguments for parameters
    for i in 1:N
        push!(fexpr.args[2].args[1].args,x.parameters[i])
    end
    # add arguments for inputnames -- local data
    for i in 1:length(x.inputnames)
        push!(fexpr.args[2].args[1].args,x.inputnames[i])
    end
    #@show fexpr

    return eval(:( let; local _EXPRREV_; $fexpr; end;))
end

function getvalue(x::ParametricExpressionWithParams,values::Vector)
    if !haskey(exprval_cache, x.p.hashval)
        exprval_cache[x.p.hashval] = genexprval(x.p)
    end
    return exprval_cache[x.p.hashval](values, IdentityArray(), x.params..., x.p.inputvals...)
end

getvalue(x::ParametricExpression{0},values::Vector) = getvalue(x[],values)

function genfval_simple(x::SymbolicOutput,num_total_vars::Int=typemax(Int))
    if num_total_vars == typemax(Int)
        # TODO: deprecate this version
        prepare_indexmap(x,Set{Int}())
    else
        prepare_indexmap(x,IndexedSet(num_total_vars))
    end
    f = genfval_parametric(x)
    return xvals -> f(xvals, IdentityArray(), x.inputvals...)
end

type IdentityArray
end

Base.getindex(::IdentityArray,i) = i

function genfgrad_simple(x::SymbolicOutput,num_total_vars::Int=typemax(Int))
    if num_total_vars == typemax(Int)
        # TODO: deprecate this version
        prepare_indexmap(x,Set{Int}())
    else
        prepare_indexmap(x,IndexedSet(num_total_vars))
    end
    f = genfgrad_parametric(x)
    return (xvals, out) -> (fill!(out, 0.0); f(xvals, IdentityArray(), out, IdentityArray(),x.inputvals...))
end

export genfgrad, genfgrad_simple, genfgrad_parametric, genfval_parametric, genfval_simple
