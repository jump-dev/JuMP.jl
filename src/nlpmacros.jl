#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.




# generates code which converts an expression into a NodeData array (tape)
# parent is the index of the parent expression
# values is the name of the list of constants which appear in the expression
function parseNLExpr(m, x, tapevar, parent, values)

    if isexpr(x, :call)
        if length(x.args) == 2 # univariate
            code = :(let; end)
            block = code.args[1]
            @assert isexpr(block, :block)
            haskey(univariate_operator_to_id,x.args[1]) || error("Unrecognized function $(x.args[1]) used in nonlinear expression.")
            operatorid = univariate_operator_to_id[x.args[1]]
            push!(block.args, :(push!($tapevar, NodeData(CALLUNIVAR, $operatorid, $parent))))
            parentvar = gensym()
            push!(block.args, :($parentvar = length($tapevar)))
            push!(block.args, parseNLExpr(m, x.args[2], tapevar, parentvar, values))
            return code
        else
            code = :(let; end)
            block = code.args[1]
            @assert isexpr(block, :block)
            haskey(operator_to_id,x.args[1]) || error("Unrecognized function $(x.args[1]) used in nonlinear expression.")
            operatorid = operator_to_id[x.args[1]]
            parentvar = gensym()
            push!(block.args, :(push!($tapevar, NodeData(CALL, $operatorid, $parent))))
            push!(block.args, :($parentvar = length($tapevar)))
            for i in 1:length(x.args)-1
                push!(block.args, parseNLExpr(m, x.args[i+1], tapevar, parentvar, values))
            end
            return code
        end
    end
    if isexpr(x, :comparison)
        code = :(let; end)
        block = code.args[1]
        op = x.args[2]
        operatorid = comparison_operator_to_id[op]
        for k in 2:2:length(x.args)-1
            @assert x.args[k] == op # don't handle a <= b >= c
        end
        parentvar = gensym()
        push!(block.args, :(push!($tapevar, NodeData(COMPARISON, $operatorid, $parent))))
        push!(block.args, :($parentvar = length($tapevar)))
        for k in 1:2:length(x.args)
            push!(block.args, parseNLExpr(m, x.args[k], tapevar, parentvar, values))
        end
        return code
    end
    if isexpr(x, :&&) || isexpr(x, :||)
        code = :(let; end)
        block = code.args[1]
        op = x.head
        operatorid = logic_operator_to_id[op]
        parentvar = gensym()
        push!(block.args, :(push!($tapevar, NodeData(LOGIC, $operatorid, $parent))))
        push!(block.args, :($parentvar = length($tapevar)))
        push!(block.args, parseNLExpr(m, x.args[1], tapevar, parentvar, values))
        push!(block.args, parseNLExpr(m, x.args[2], tapevar, parentvar, values))
        return code
    end
    if isexpr(x, :curly)
        header = x.args[1]
        if length(x.args) < 3
            error("Need at least two arguments for $header")
        end
        if issum(header)
            operatorid = operator_to_id[:+]
        elseif isprod(header)
            operatorid = operator_to_id[:*]
        else
            error("Unrecognized expression $header{...}")
        end
        codeblock = :(let; end)
        block = codeblock.args[1]
        push!(block.args, :(push!($tapevar, NodeData(CALL, $operatorid, $parent))))
        parentvar = gensym()
        push!(block.args, :($parentvar = length($tapevar)))
        
        # we have a filter condition
        if isexpr(x.args[2],:parameters)
            cond = x.args[2]
            if length(cond.args) != 1
                error("No commas after semicolon allowed in $header{} expression, use && for multiple conditions")
            end
            # generate inner loop code first and then wrap in for loops
            innercode = parseNLExpr(m, x.args[3], tapevar, parentvar, values)
            code = quote
                if $(esc(cond.args[1]))
                    $innercode
                end
            end
            for level in length(x.args):-1:4
                _idxvar, idxset = parseIdxSet(x.args[level]::Expr)
                idxvar = esc(_idxvar)
                code = :(let
                    $(localvar(idxvar))
                    for $idxvar in $(esc(idxset))
                        $code
                    end
                end)
            end
            push!(block.args, code)
        else # no condition
            innercode = parseNLExpr(m, x.args[2], tapevar, parentvar, values)
            code = quote
                $innercode
            end
            for level in length(x.args):-1:3
                _idxvar, idxset = parseIdxSet(x.args[level]::Expr)
                idxvar = esc(_idxvar)
                code = :(let
                    $(localvar(idxvar))
                    for $idxvar in $(esc(idxset))
                        $code
                    end
                end)
            end
            push!(block.args, code)
        end
        return codeblock
    end
    # at the lowest level?
    return :( parseNLExpr_runtime($(esc(m)),$(esc(x)), $tapevar, $parent, $values) )

end

function parseNLExpr_runtime(m::Model, x::Number, tape, parent, values)
    push!(values, x)
    push!(tape, NodeData(VALUE, length(values), parent))
    nothing
end

# Temporary hack for deprecation of @defNLExpr syntax
const __last_model = Array(Model,1)

function parseNLExpr_runtime(m::Model, x::Variable, tape, parent, values)
    __last_model[1] = x.m
    x.m === m || error("Variable in nonlinear expression does not belong to corresponding model")
    push!(tape, NodeData(VARIABLE, x.col, parent))
    nothing
end

function parseNLExpr_runtime(m::Model, x::NonlinearExpression, tape, parent, values)
    push!(tape, NodeData(SUBEXPRESSION, x.index, parent))
    nothing
end

function parseNLExpr_runtime(m::Model, x::NonlinearParameter, tape, parent, values)
    push!(tape, NodeData(PARAMETER, x.index, parent))
    nothing
end

macro processNLExpr(m, ex)
    parsed = parseNLExpr(m, ex, :tape, -1, :values)
    quote
        tape = NodeData[]
        values = Float64[]
        $parsed
        NonlinearExprData(tape, values)
    end
end
