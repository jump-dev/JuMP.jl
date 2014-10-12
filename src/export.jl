# flatten out expressions for output to MathProgBase format

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
        return x.value
    elseif isa(x.ex, Expr) || isa(x.ex, Symbol)
        # some other symbolic value, to evaluate at runtime
        push!(expr_out.args, :( $(x.value) = outputvalue($(x.ex) )))
        return x.value
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
    expr = outputpass(x.tree, out)
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
    fval = outputpass(x.tree, out)
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
