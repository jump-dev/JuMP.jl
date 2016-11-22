#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.


# generates code which converts an expression into a NodeData array (tape)
# parent is the index of the parent expression
# values is the name of the list of constants which appear in the expression
function parseNLExpr(m, x, tapevar, parent, values)

    if isexpr(x,:call) && length(x.args) >= 2 && (isexpr(x.args[2],:generator) || isexpr(x.args[2],:flatten))
        header = x.args[1]
        if issum(header)
            operatorid = operator_to_id[:+]
        elseif isprod(header)
            operatorid = operator_to_id[:*]
        else
            error("Unrecognized expression $header(...)")
        end
        codeblock = :(let; end)
        block = codeblock.args[1]
        push!(block.args, :(push!($tapevar, NodeData(CALL, $operatorid, $parent))))
        parentvar = gensym()
        push!(block.args, :($parentvar = length($tapevar)))


        code = parsegen(x.args[2], t -> parseNLExpr(m, t, tapevar, parentvar, values))
        push!(block.args, code)
        return codeblock
    end

    if isexpr(x, :call)
        if length(x.args) == 2 # univariate
            code = :(let; end)
            block = code.args[1]
            @assert isexpr(block, :block)
            if haskey(univariate_operator_to_id,x.args[1])
                operatorid = univariate_operator_to_id[x.args[1]]
                push!(block.args, :(push!($tapevar, NodeData(CALLUNIVAR, $operatorid, $parent))))
            else
                opname = quot(x.args[1])
                errorstring = "Unrecognized function $opname used in nonlinear expression."
                errorstring2 = "Incorrect number of arguments for $opname in nonlinear expression."
                lookupcode = quote
                    if !haskey(univariate_operator_to_id,$opname)
                        if haskey(operator_to_id,$opname)
                            error($errorstring2)
                        else
                            error($errorstring)
                        end
                    end
                    operatorid = univariate_operator_to_id[$opname]
                end
                push!(block.args, :($lookupcode; push!($tapevar, NodeData(CALLUNIVAR, operatorid, $parent))))
            end
            parentvar = gensym()
            push!(block.args, :($parentvar = length($tapevar)))
            push!(block.args, parseNLExpr(m, x.args[2], tapevar, parentvar, values))
            return code
        else
            code = :(let; end)
            block = code.args[1]
            @assert isexpr(block, :block)
            if haskey(operator_to_id,x.args[1]) # fast compile-time lookup
                operatorid = operator_to_id[x.args[1]]
                push!(block.args, :(push!($tapevar, NodeData(CALL, $operatorid, $parent))))
            elseif haskey(comparison_operator_to_id,x.args[1])
                operatorid = comparison_operator_to_id[x.args[1]]
                push!(block.args, :(push!($tapevar, NodeData(COMPARISON, $operatorid, $parent))))
            else # could be user defined
                opname = quot(x.args[1])
                errorstring = "Unrecognized function $opname used in nonlinear expression."
                errorstring2 = "Incorrect number of arguments for $opname in nonlinear expression."
                lookupcode = quote
                    if !haskey(operator_to_id,$opname)
                        if haskey(univariate_operator_to_id,$opname)
                            error($errorstring2)
                        else
                            error($errorstring)
                        end
                    end
                    operatorid = operator_to_id[$opname]
                end
                push!(block.args, :($lookupcode; push!($tapevar, NodeData(CALL, operatorid, $parent))))
            end
            parentvar = gensym()
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
        warn_curly(x)
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

function parseNLExpr_runtime(m::Model, x::Vector, tape, parent, values)
    error("Unexpected vector $x in nonlinear expression. Nonlinear expressions may contain only scalar expressions.")
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

# Construct a NonlinearExprData from a Julia expression.
# Variable objects should be spliced into the expression.
function NonlinearExprData(m::Model, ex::Expr)
    ex = spliceref(m,ex)
    nd, values = ReverseDiffSparse.expr_to_nodedata(ex)
    return NonlinearExprData(nd, values)
end
NonlinearExprData(m::Model, ex) = NonlinearExprData(m, :($ex + 0))

# recursively replace Variable(m, i) with Expr(:ref,:x,i) in ex
function spliceref(m::Model, ex::Expr)
    if ex.head == :ref # if we have x[1] already in there, something is wrong
        error("Unrecognized expression $ex.")
    end
    return Expr(ex.head,map(e -> spliceref(m,e), ex.args)...)
end
function spliceref(m::Model, v::Variable)
    v.m === m || error("Variable $v does not belong to this model")
    return Expr(:ref, :x, linearindex(v))
end
spliceref(m::Model, ex) = ex
