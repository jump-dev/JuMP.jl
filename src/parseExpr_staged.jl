#addToExpression(ex::Number, c::Number) = ex + c

addToExpression(ex::Number, c::Number, x::Number) = ex + c*x

addToExpression(ex::Number, c::Number, x::Variable) = AffExpr([x],[c],ex)

function addToExpression(ex::Number, c::Number, x::AffExpr)
    scale!(x.coeffs, c)
    x.constant += ex
    return x
end

addToExpression(ex::Number, c::Variable, x::Variable) = QuadExpr([c],[x],[1.0],AffExpr())

function addToExpression(ex::Number, c::AffExpr, x::AffExpr)
    q = c*x
    q.aff.constant += ex
    return q
end

function addToExpression(aff::AffExpr, c::Number, x::Number)
    aff.constant += c*x
    return aff
end

function addToExpression(aff::AffExpr,c::Number,x::Variable)
    push!(aff,convert(Float64,c),x)
    return aff
end

function addToExpression(aff::AffExpr,c::Variable,x::Variable)
    q = QuadExpr([c],[x],[1.0],aff)
    return q
end

function addToExpression(aff::AffExpr,c::Number,x::AffExpr)
    append!(aff.vars, x.vars)
    append!(aff.coeffs, c*x.coeffs)
    aff.constant += c*x.constant
    return aff
end

addToExpression(aff::AffExpr,x::AffExpr,c::Variable) = addToExpression(aff,c,x)
addToExpression(aff::AffExpr,c::Variable,x::AffExpr) = QuadExpr(fill(c,length(x.vars)),
                                                                x.vars,
                                                                x.coeffs,
                                                                addToExpression(aff,x.constant,c))

addToExpression(aff::AffExpr, x::QuadExpr, c::Number) = addToExpression(aff,c,x)

addToExpression(aff::AffExpr, c::Number, x::QuadExpr) = QuadExpr(copy(x.qvars1),
                                                                 copy(x.qvars2), 
                                                                 c*x.qcoeffs, 
                                                                 addToExpression(aff,c,x.aff))

function addToExpression(ex::AffExpr, c::AffExpr, x::AffExpr)
    q = c*x
    return addToExpression(q,1.0,ex)
end

function addToExpression(quad::QuadExpr,c::Number,x::Variable)
    push!(quad.aff,convert(Float64,c),x)
    return quad
end

function addToExpression(quad::QuadExpr,c::Number,x::Number)
    (quad.aff.constant += c*x)
    return quad
end

function addToExpression(quad::QuadExpr,c::Variable,x::Variable)
    push!(quad.qvars1, x)
    push!(quad.qvars2, c)
    return quad
end

function addToExpression(quad::QuadExpr,c::Number,x::AffExpr)
    append!(quad.aff.vars, x.vars)
    append!(quad.aff.coeffs, c*x.coeffs)
    quad.aff.constant += c*x.constant
    return quad
end

addToExpression(quad::QuadExpr,x::AffExpr,c::Variable) = addToExpression(quad,c,x)
function addToExpression(quad::QuadExpr,c::Variable,x::AffExpr)
    append!(quad.qvars1, fill(c,length(x.vars)))
    append!(quad.qvars2, x.vars)
    append!(quad.qcoeffs, x.coeffs)
    addToExpression(quad.aff,x.constant,c)
    return quad
end

addToExpression(quad::QuadExpr,x::QuadExpr,c::Number) = addToExpression(quad,c,x)
function addToExpression(quad::QuadExpr,c::Number,x::QuadExpr)
    append!(quad.qvars1,x.qvars1)
    append!(quad.qvars2,x.qvars2)
    append!(quad.qcoeffs,c*x.qcoeffs)
    addToExpression(quad,c,x.aff)
    return quad
end

function addToExpression(ex::QuadExpr, c::AffExpr, x::AffExpr)
    q = c*x
    return addToExpression(ex, 1.0, q)
end

addToExpression(aff, c, x) = error("Cannot construct an affine expression with a term of type ($(typeof(c)))*($(typeof(x)))")

stagedfunction addToExpression_reorder(ex, args...)
    if !isleaftype(ex) || mapreduce(t -> !isleaftype(t), |, args)
        error("Can't process abstract types")
    end
    # how many of the multiplicands are variables?
    n_var = mapreduce(t -> (t == Variable || t == AffExpr), +, args)
    has_quad = mapreduce(t -> (t == QuadExpr), |, args)
    if n_var == 0 && !has_quad
        #println("No variables")
        return :(addToExpression(ex, 1.0, (*)(args...)))
    elseif n_var == 1 && !has_quad # linear
        #println("Linear")
        coef_expr = Expr(:call, :*)
        coef_idx = Int[]
        var_idx = 0
        for i in 1:length(args)
            if args[i] == Variable || args[i] == AffExpr
                var_idx = i
            else
                push!(coef_expr.args, :(args[$i]))
            end
        end
        return :(addToExpression(ex, $coef_expr, args[$var_idx]))
    else
        #println("Nonlinear")
        # fall back
        coef_expr = Expr(:call, :*, [:(args[$i]) for i in 1:(length(args)-1)]...)
        return :(addToExpression(ex, $coef_expr, args[$(length(args))]))
    end


end

function parseCurly(x::Expr, aff::Symbol, constantCoef)
    if !(x.args[1] == :sum || x.args[1] == :∑ || x.args[1] == :Σ) # allow either N-ARY SUMMATION or GREEK CAPITAL LETTER SIGMA
        error("Expected sum outside curly braces")
    end
    if length(x.args) < 3
        error("Need at least two arguments for sum")
    end

    # we have a filter condition
    if isexpr(x.args[2],:parameters)
        cond = x.args[2]
        if length(cond.args) != 1
            error("No commas after semicolon allowed in sum expression, use && for multiple conditions")
        end
        # generate inner loop code first and then wrap in for loops
        newaff, innercode = parseExpr(x.args[3], aff, constantCoef, aff)
        @assert aff == newaff
        code = quote
            if $(esc(cond.args[1]))
                $innercode
            end
        end
        for level in length(x.args):-1:4
            code = :(
            for $(esc(x.args[level].args[1])) in $(esc(x.args[level].args[2]))
                $code
            end)
        end
    else # no condition
        newaff, code = parseExpr(x.args[2], aff, constantCoef, aff)
        for level in length(x.args):-1:3
            code = :(
            for $(esc(x.args[level].args[1])) in $(esc(x.args[level].args[2]))
                $code
            end)
        end
        len = :len
        # precompute the number of elements to add
        # this is unncessary if we're just summing constants
        preblock = :($len += length($(esc(x.args[length(x.args)].args[2]))))
        for level in (length(x.args)-1):-1:3
            preblock = Expr(:for, esc(x.args[level]),preblock)
        end
        preblock = quote
            $len = 0
            $preblock
            if isa($aff,GenericAffExpr)
                sizehint($aff.vars,length($aff.vars)+$len)
                sizehint($aff.coeffs,length($aff.coeffs)+$len)
            elseif isa($aff,GenericQuadExpr)
                sizehint($aff.qvars1,length($aff.qvars1)+$len)
                sizehint($aff.qvars2,length($aff.qvars2)+$len)
                sizehint($aff.qcoeffs,length($aff.qcoeffs)+$len)
                sizehint($aff.aff.vars,length($aff.aff.vars)+$len)
                sizehint($aff.aff.coeffs,length($aff.aff.coeffs)+$len)
            end
        end
        code = :($preblock;$code)
    end

    
    return code
end

is_complex_expr(ex) = isa(ex,Expr) && !isexpr(ex,:ref)

# output is assigned to newaff
function parseExpr(x, aff::Symbol, coefficients::Vector, newaff::Symbol=gensym())
    #@show x
    #@show coefficients
    if !isa(x,Expr)
        # at the lowest level
        callexpr = Expr(:call,:addToExpression_reorder,aff,esc(x),coefficients...)
        return newaff, :($newaff = $callexpr)
    else
        if x.head == :call && x.args[1] == :+
            b = Expr(:block)
            aff_ = aff
            for arg in x.args[2:(end-1)]
                aff_, code = parseExpr(arg, aff_, coefficients)
                push!(b.args, code)
            end
            newaff, code = parseExpr(x.args[end], aff_, coefficients, newaff)
            push!(b.args, code)
            return newaff, b
        elseif x.head == :call && x.args[1] == :-
            if length(x.args) == 2 # unary subtraction
                return parseExpr(x.args[2], aff, vcat(coefficients, -1.0), newaff)
            else # a - b - c ...
                b = Expr(:block)
                aff_, code = parseExpr(x.args[2], aff, coefficients)
                push!(b.args, code)
                for arg in x.args[3:(end-1)]
                    aff_,code = parseExpr(arg, aff_, vcat(coefficients, -1.0))
                    push!(b.args, code)
                end
                newaff,code = parseExpr(x.args[end], aff_, vcat(coefficients, -1.0), newaff)
                push!(b.args, code)
                return newaff, b
            end
        elseif x.head == :call && x.args[1] == :*
            # we might need to recurse on multiple arguments, e.g.,
            # (x+y)*(x+y)
            n_expr = mapreduce(is_complex_expr, +, x.args)
            if n_expr == 1 # special case, only recurse on one argument and don't create temporary objects
                coefficients = Any[c for c in coefficients]
                which_idx = 0
                for i in 2:length(x.args)
                    if is_complex_expr(x.args[i])
                        which_idx = i
                    else
                        push!(coefficients, esc(x.args[i]))
                    end
                end
                return parseExpr(x.args[which_idx], aff, coefficients, newaff)
            else
                blk = Expr(:block)
                for i in 2:length(x.args)
                    if is_complex_expr(x.args[i])
                        s = gensym()
                        newaff_, parsed = parseExpr(x.args[i], s, [1.0])
                        push!(blk.args, :($s = 0.0; $parsed))
                        x.args[i] = newaff_
                    else
                        x.args[i] = esc(x.args[i])
                    end
                end
                callexpr = Expr(:call,:addToExpression_reorder,aff,vcat(coefficients,x.args[2:end])...)
                push!(blk.args, :($newaff = $callexpr))
                return newaff, blk
            end
        elseif x.head == :call && x.args[1] == :^ && is_complex_expr(x.args[2])
            x.args[3] == 2 || error("Only exponents of 2 are currently supported")
            blk = Expr(:block)
            s = gensym()
            newaff_, parsed = parseExpr(x.args[2], s, [1.0])
            push!(blk.args, :($s = 0.0; $parsed))
            push!(blk.args, :($newaff = $newaff_*$newaff_))
            return newaff, blk
        elseif x.head == :call && x.args[1] == :/
            @assert length(x.args) == 3
            numerator = x.args[2]
            denom = x.args[3]
            return parseExpr(numerator, aff, vcat(coefficients,esc(:(1/$denom))),newaff)
        elseif x.head == :curly
            return newaff, :($(parseCurly(x,aff,coefficients)); $newaff = $aff)
        else # at lowest level?
            !isexpr(x,:comparison) || error("Unexpected comparison in expression $x") 
            callexpr = Expr(:call,:addToExpression_reorder,aff,esc(x),coefficients...)
            return newaff, :($newaff = $callexpr)
        end
    end
end
