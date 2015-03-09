# flatten out expressions for output to MathProgBase format

# Lots of overlap with this code and forwardpass(). Possible to combine?
function outputpass(x::ExprNode, expr_out)
    @assert isexpr(expr_out, :block)

    # returns the variable name containing the result

    x.value = gensym() # generate a variable to represent the "value" of this node
    if isexpr(x.ex, :call)
        values = Any[]
        for i in 2:length(x.ex.args)
            push!(values, outputpass(x.ex.args[i], expr_out))
        end
        fcall = Expr(:call, :Expr, quot(:call), quot(x.ex.args[1]), values...)
        push!(expr_out.args, :( $(x.value) = $fcall ))
        return x.value
    elseif isexpr(x.ex, :curly)
        oper = x.ex.args[1]
        @assert issum(oper) || isprod(oper)
        # compute value of this node, need to use a loop
        if issum(oper)
            push!(expr_out.args, :( $(x.value) = Expr(:call,:+) ))
        else # :prod
            push!(expr_out.args, :( $(x.value) = Expr(:call,:*) ))
        end
        code = quote end

        valexpr = outputpass(curlyexpr(x.ex), code)

        code = :( $code; push!($(x.value).args, $valexpr))

        push!(expr_out.args, gencurlyloop(x.ex, code))
        # perform some basic simplifications
        defaultvalue = issum(oper) ? 0 : 1
        push!(expr_out.args, quote
            if length($(x.value).args) == 1
                $(x.value) = $defaultvalue
            elseif length($(x.value).args) == 2
                $(x.value) = $(x.value).args[2]
            end
        end)
        return x.value
    elseif isexpr(x.ex, :comparison)
        values = Any[]
        for i in 1:length(x.ex.args)
            if iseven(i)
                push!(values, quot(x.ex.args[i])) # comparison operator
            else
                push!(values, outputpass(x.ex.args[i], expr_out))
            end
        end
        fcall = Expr(:call, :Expr, quot(:comparison), values...)
        push!(expr_out.args, :( $(x.value) = $fcall ))
        return x.value
    elseif isexpr(x.ex, :&&) || isexpr(x.ex, :||)
        values = Any[]
        for i in 1:length(x.ex.args)
            push!(values, outputpass(x.ex.args[i], expr_out))
        end
        fcall = Expr(:call, :Expr, quot(x.ex.head), values...)
        push!(expr_out.args, :( $(x.value) = $fcall ))
        return x.value
    elseif isa(x.ex, Expr) || isa(x.ex, Symbol)
        # some other symbolic value, to evaluate at runtime
        push!(expr_out.args, :( $(x.value) = outputvalue($(x.ex) )))
        return x.value
    elseif isa(x.ex,ParametricExpressionWithParams)
        # just splat
        tree = x.ex.p.tree
        for i in 1:length(x.ex.params)
            tree = replace_param(tree, x.ex.p.parameters[i], x.ex.params[i])
        end
        outputpass(genExprGraph(tree), expr_out)
    else
        error("Shouldn't get here")
    end
end

outputpass(x, expr_out) = :(outputvalue($x))

outputvalue(x::Placeholder) = Expr(:ref,:x,getplaceindex(x))
outputvalue(x) = x

function to_flat_expr(x::SymbolicOutput)
    out = Expr(:block)
    # load data into local scope
    for i in 1:length(x.inputnames)
        push!( out.args, :( $(x.inputnames[i]) = $(x.inputvals[i]) ))
    end
    expr = outputpass(genExprGraph(x.tree), out)
    fexpr = quote
        let
            $out
            $expr
        end
    end
    return eval(fexpr)

end

export to_flat_expr

function genfexpr_parametric(x::SymbolicOutput)
    out = Expr(:block)
    fval = outputpass(genExprGraph(x.tree), out)
    fname = gensym()
    fexpr = quote
        function $(fname)()
            $out
            return $fval
        end
    end
    # add arguments for inputnames -- local data
    for i in 1:length(x.inputnames)
        push!(fexpr.args[2].args[1].args,x.inputnames[i])
    end

    return eval(fexpr)
end
