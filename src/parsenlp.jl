#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
        if issum(x.args[1]) || isprod(x.args[1])
            opname = x.args[1]
            errorstring = "$opname() can appear in nonlinear expressions " *
            " only if the argument is a generator statement, for example, " *
            "$opname(x[i] for i in 1:N)."
            return :(error($errorstring))
        end
        if length(x.args) == 2 # univariate
            code = :(let; end)
            block = code.args[1]
            @assert isexpr(block, :block)
            if haskey(univariate_operator_to_id,x.args[1])
                operatorid = univariate_operator_to_id[x.args[1]]
                push!(block.args, :(push!($tapevar, NodeData(CALLUNIVAR, $operatorid, $parent))))
            else
                opname = quot(x.args[1])
                errorstring = "Unrecognized function \"$(x.args[1])\" used in nonlinear expression."
                errorstring2 = "Incorrect number of arguments for \"$(x.args[1])\" in nonlinear expression."
                lookupcode = quote
                    if $(esc(m)).nlpdata === nothing
                        error($errorstring)
                    end
                    if !haskey($(esc(m)).nlpdata.user_operators.univariate_operator_to_id,$opname)
                        if haskey($(esc(m)).nlpdata.user_operators.multivariate_operator_to_id,$opname)
                            error($errorstring2)
                        else
                            error($errorstring)
                        end
                    end
                    operatorid = $(esc(m)).nlpdata.user_operators.univariate_operator_to_id[$opname] + Derivatives.USER_UNIVAR_OPERATOR_ID_START - 1
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
                errorstring = "Unrecognized function \"$(x.args[1])\" used in nonlinear expression."
                errorstring2 = "Incorrect number of arguments for \"$(x.args[1])\" in nonlinear expression."
                lookupcode = quote
                    if $(esc(m)).nlpdata === nothing
                        error($errorstring)
                    end
                    if !haskey($(esc(m)).nlpdata.user_operators.multivariate_operator_to_id,$opname)
                        if haskey($(esc(m)).nlpdata.user_operators.univariate_operator_to_id,$opname)
                            error($errorstring2)
                        else
                            error($errorstring)
                        end
                    end
                    operatorid = $(esc(m)).nlpdata.user_operators.multivariate_operator_to_id[$opname] + Derivatives.USER_OPERATOR_ID_START - 1
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
        error_curly(x)
    end
    # at the lowest level?
    return :( parseNLExpr_runtime($(esc(m)),$(esc(x)), $tapevar, $parent, $values) )

end

function parseNLExpr_runtime(m::Model, x::Number, tape, parent, values)
    push!(values, x)
    push!(tape, NodeData(VALUE, length(values), parent))
    nothing
end

function parseNLExpr_runtime(m::Model, x::VariableRef, tape, parent, values)
    x.m === m || error("Variable in nonlinear expression does not belong to corresponding model")
    push!(tape, NodeData(MOIVARIABLE, x.index.value, parent))
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

function parseNLExpr_runtime(m::Model, x::AbstractArray, tape, parent, values)
    error("Unexpected array $x in nonlinear expression. Nonlinear expressions may contain only scalar expressions.")
end

function expression_complexity(ex::Expr)
    return isempty(ex.args) ? 1 : sum(expression_complexity, ex.args)
end
expression_complexity(other) = 1

macro processNLExpr(m, ex)
    # This is an arbitrary cutoff. See issue #1355.
    if expression_complexity(ex) > 5000
        Compat.@warn "Processing a very large nonlinear expression with " *
                     "@NLexpression/@NLconstraint/@NLobjective. This may be " *
                     "very slow. Consider using setNLobjective() and " *
                     "addNLconstraint() instead of the macros or " *
                     "reformulating the expressions using sum() and prod() " *
                     "to make them more compact. The macros are designed to " *
                     "process smaller, human-readable expressions."
    end
    parsed = parseNLExpr(m, ex, :tape, -1, :values)
    quote
        tape = NodeData[]
        values = Float64[]
        $parsed
        NonlinearExprData(tape, values)
    end
end

function Derivatives.expr_to_nodedata(ex::VariableRef,nd::Vector{NodeData},values::Vector{Float64},parentid,r::Derivatives.UserOperatorRegistry)
    push!(nd, NodeData(MOIVARIABLE, ex.index.value, parentid))
    nothing
end

function Derivatives.expr_to_nodedata(ex::NonlinearExpression,nd::Vector{NodeData},values::Vector{Float64},parentid,r::Derivatives.UserOperatorRegistry)
    push!(nd, NodeData(SUBEXPRESSION, ex.index, parentid))
    nothing
end

function Derivatives.expr_to_nodedata(ex::NonlinearParameter,nd::Vector{NodeData},values::Vector{Float64},parentid,r::Derivatives.UserOperatorRegistry)
    push!(nd, NodeData(PARAMETER, ex.index, parentid))
    nothing
end

# Construct a NonlinearExprData from a Julia expression.
# VariableRef objects should be spliced into the expression.
function NonlinearExprData(m::Model, ex::Expr)
    initNLP(m)
    checkexpr(m,ex)
    nd, values = Derivatives.expr_to_nodedata(ex,m.nlpdata.user_operators)
    return NonlinearExprData(nd, values)
end
NonlinearExprData(m::Model, ex) = NonlinearExprData(m, :($ex + 0))

# Error if:
# 1) Unexpected expression
# 2) VariableRef doesn't match the model
function checkexpr(m::Model, ex::Expr)
    if ex.head == :ref # if we have x[1] already in there, something is wrong
        error("Unrecognized expression $ex. JuMP variable objects and input coefficients should be spliced directly into expressions.")
    end
    for e in ex.args
        checkexpr(m, e)
    end
    return
end
function checkexpr(m::Model, v::VariableRef)
    v.m === m || error("Variable $v does not belong to this model")
    return
end
checkexpr(m::Model, ex) = nothing
