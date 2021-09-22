#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Returns the block expression inside a :let that holds the code to be run.
# The other block (not returned) is for declaring variables in the scope of the
# let.
function _let_code_block(ex::Expr)
    @assert isexpr(ex, :let)
    return ex.args[2]
end

function _error_curly(x)
    return Base.error(
        "The curly syntax (sum{},prod{},norm2{}) is no longer supported. Expression: $x.",
    )
end

function _warn_auto_register(opname, N)
    @warn("""Function $(opname) automatically registered with $N arguments.

    Calling the function with a different number of arguments will result in an
    error.

    While you can safely ignore this warning, we recommend that you manually
    register the function as follows:
    ```Julia
    model = Model()
    register(model, :$opname, $N, $opname; autodiff = true)
    ```""")
    return
end

# generates code which converts an expression into a NodeData array (tape)
# parent is the index of the parent expression
# values is the name of the list of constants which appear in the expression
function _parse_NL_expr(m, x, tapevar, parent, values)
    if isexpr(x, :block)
        error(
            "`begin...end` blocks are not supported in nonlinear macros. The " *
            "nonlinear expression must be a single statement.",
        )
    end
    if isexpr(x, :call) &&
       length(x.args) >= 2 &&
       (isexpr(x.args[2], :generator) || isexpr(x.args[2], :flatten))
        header = x.args[1]
        if _is_sum(header)
            operatorid = operator_to_id[:+]
        elseif _is_prod(header)
            operatorid = operator_to_id[:*]
        else
            error("Unrecognized expression $header(...)")
        end
        codeblock = :(
            let
            end
        )
        block = _let_code_block(codeblock)
        push!(
            block.args,
            :(push!($tapevar, NodeData(CALL, $operatorid, $parent))),
        )
        parentvar = gensym()
        push!(block.args, :($parentvar = length($tapevar)))

        code = _MA.rewrite_generator(
            x.args[2],
            t -> _parse_NL_expr(m, t, tapevar, parentvar, values),
        )
        push!(block.args, code)
        return codeblock
    end

    if isexpr(x, :call)
        if isexpr(x.args[1], :.)
            # Functions like foo.bar cannot possibly be registered, because you
            # can register only with a symbol name.
            errorstring =
                "Unexpected function $(x.args[1]). See the " *
                "documentation on how to register a function."
            return :(error($errorstring))
        end
        if _is_sum(x.args[1]) || _is_prod(x.args[1])
            opname = x.args[1]
            errorstring =
                "$opname() can appear in nonlinear expressions " *
                " only if the argument is a generator statement, for example, " *
                "$opname(x[i] for i in 1:N)."
            return :(error($errorstring))
        end
        if length(x.args) == 2 && !isexpr(x.args[2], :...) # univariate
            code = :(
                let
                end
            )
            block = _let_code_block(code)
            @assert isexpr(block, :block)
            if haskey(univariate_operator_to_id, x.args[1])
                operatorid = univariate_operator_to_id[x.args[1]]
                push!(
                    block.args,
                    :(push!(
                        $tapevar,
                        NodeData(CALLUNIVAR, $operatorid, $parent),
                    )),
                )
            else
                opname = quot(x.args[1])
                f = x.args[1]
                errorstring = """
                Unrecognized function \"$(f)\" used in nonlinear expression.

                If the function exists, but is not within the scope of this call,
                you should register it as a user-defined function before building
                the model. For example:
                ```julia
                model = Model()
                register(model, :$(f), 1, $(f), autodiff=true)
                # ... variables and constraints ...
                ```
                """
                errorstring2 = "Incorrect number of arguments for \"$(x.args[1])\" in nonlinear expression."
                lookupcode = quote
                    if $(esc(m)).nlp_data === nothing
                        try
                            register(
                                $(esc(m)),
                                $opname,
                                1,
                                $(esc(x.args[1]));
                                autodiff = true,
                            )
                            _warn_auto_register($opname, 1)
                        catch
                            error($errorstring)
                        end
                    end
                    if !haskey(
                        $(esc(
                            m,
                        )).nlp_data.user_operators.univariate_operator_to_id,
                        $opname,
                    )
                        if haskey(
                            $(esc(
                                m,
                            )).nlp_data.user_operators.multivariate_operator_to_id,
                            $opname,
                        )
                            error($errorstring2)
                        else
                            try
                                register(
                                    $(esc(m)),
                                    $opname,
                                    1,
                                    $(esc(x.args[1]));
                                    autodiff = true,
                                )
                                _warn_auto_register($opname, 1)
                            catch
                                error($errorstring)
                            end
                        end
                    end
                    operatorid =
                        $(esc(
                            m,
                        )).nlp_data.user_operators.univariate_operator_to_id[$opname] +
                        _Derivatives.USER_UNIVAR_OPERATOR_ID_START - 1
                end
                push!(
                    block.args,
                    :(
                        $lookupcode;
                        push!(
                            $tapevar,
                            NodeData(CALLUNIVAR, operatorid, $parent),
                        )
                    ),
                )
            end
            parentvar = gensym()
            push!(block.args, :($parentvar = length($tapevar)))
            push!(
                block.args,
                _parse_NL_expr(m, x.args[2], tapevar, parentvar, values),
            )
            return code
        else
            code = :(
                let
                end
            )
            block = _let_code_block(code)
            @assert isexpr(block, :block)
            if haskey(operator_to_id, x.args[1]) # fast compile-time lookup
                operatorid = operator_to_id[x.args[1]]
                push!(
                    block.args,
                    :(push!($tapevar, NodeData(CALL, $operatorid, $parent))),
                )
            elseif haskey(comparison_operator_to_id, x.args[1])
                operatorid = comparison_operator_to_id[x.args[1]]
                push!(
                    block.args,
                    :(push!(
                        $tapevar,
                        NodeData(COMPARISON, $operatorid, $parent),
                    )),
                )
            else # could be user defined
                opname = quot(x.args[1])
                N = length(x.args) - 1
                f = x.args[1]
                errorstring = """
                Unrecognized function \"$(f)\" used in nonlinear expression.

                If the function exists, but is not within the scope of this call,
                you should register it as a user-defined function before building
                the model. For example:
                ```julia
                model = Model()
                register(model, :$(f), $(N), $(f), autodiff=true)
                # ... variables and constraints ...
                ```
                """
                errorstring2 = "Incorrect number of arguments for \"$(x.args[1])\" in nonlinear expression."
                lookupcode = quote
                    if $(esc(m)).nlp_data === nothing
                        try
                            register(
                                $(esc(m)),
                                $opname,
                                $N,
                                $(esc(x.args[1]));
                                autodiff = true,
                            )
                            _warn_auto_register($opname, $N)
                        catch
                            error($errorstring)
                        end
                    end
                    if !haskey(
                        $(esc(
                            m,
                        )).nlp_data.user_operators.multivariate_operator_to_id,
                        $opname,
                    )
                        if haskey(
                            $(esc(
                                m,
                            )).nlp_data.user_operators.univariate_operator_to_id,
                            $opname,
                        )
                            error($errorstring2)
                        else
                            try
                                register(
                                    $(esc(m)),
                                    $opname,
                                    $N,
                                    $(esc(x.args[1]));
                                    autodiff = true,
                                )
                                _warn_auto_register($opname, $N)
                            catch
                                error($errorstring)
                            end
                        end
                    end
                    operatorid =
                        $(esc(
                            m,
                        )).nlp_data.user_operators.multivariate_operator_to_id[$opname] +
                        _Derivatives.USER_OPERATOR_ID_START - 1
                end
                push!(
                    block.args,
                    :($lookupcode;
                    push!($tapevar, NodeData(CALL, operatorid, $parent))),
                )
            end
            parentvar = gensym()
            push!(block.args, :($parentvar = length($tapevar)))
            for i in 1:length(x.args)-1
                arg = x.args[i+1]
                if isexpr(arg, :...)
                    if !isa(arg.args[1], Symbol)
                        error(
                            "Unexpected expression in $x. JuMP supports " *
                            "splatting only symbols. For example, x... is " *
                            "ok, but (x + 1)..., [x; y]... and g(f(y)...) " *
                            "are not.",
                        )
                    end
                    push!(
                        block.args,
                        quote
                            for val in $(esc(arg.args[1]))
                                _parse_NL_expr_runtime(
                                    $(esc(m)),
                                    val,
                                    $tapevar,
                                    $parentvar,
                                    $values,
                                )
                            end
                        end,
                    )
                else
                    push!(
                        block.args,
                        _parse_NL_expr(m, arg, tapevar, parentvar, values),
                    )
                end
            end
            return code
        end
    end
    if isexpr(x, :comparison)
        code = :(
            let
            end
        )
        block = _let_code_block(code)
        op = x.args[2]
        operatorid = comparison_operator_to_id[op]
        for k in 2:2:length(x.args)-1
            @assert x.args[k] == op # don't handle a <= b >= c
        end
        parentvar = gensym()
        push!(
            block.args,
            :(push!($tapevar, NodeData(COMPARISON, $operatorid, $parent))),
        )
        push!(block.args, :($parentvar = length($tapevar)))
        for k in 1:2:length(x.args)
            push!(
                block.args,
                _parse_NL_expr(m, x.args[k], tapevar, parentvar, values),
            )
        end
        return code
    end
    if isexpr(x, :&&) || isexpr(x, :||)
        code = :(
            let
            end
        )
        block = _let_code_block(code)
        op = x.head
        operatorid = logic_operator_to_id[op]
        parentvar = gensym()
        push!(
            block.args,
            :(push!($tapevar, NodeData(LOGIC, $operatorid, $parent))),
        )
        push!(block.args, :($parentvar = length($tapevar)))
        push!(
            block.args,
            _parse_NL_expr(m, x.args[1], tapevar, parentvar, values),
        )
        push!(
            block.args,
            _parse_NL_expr(m, x.args[2], tapevar, parentvar, values),
        )
        return code
    end
    if isexpr(x, :curly)
        _error_curly(x)
    end
    if isexpr(x, :...)
        error("Unexpected splatting expression $x.")
    end
    # at the lowest level?
    return :(_parse_NL_expr_runtime(
        $(esc(m)),
        $(esc(x)),
        $tapevar,
        $parent,
        $values,
    ))
end

function _parse_NL_expr_runtime(m::Model, x::Real, tape, parent, values)
    push!(values, x)
    push!(tape, NodeData(VALUE, length(values), parent))
    return nothing
end

function _parse_NL_expr_runtime(m::Model, x::VariableRef, tape, parent, values)
    if owner_model(x) !== m
        error(
            "Variable in nonlinear expression does not belong to the " *
            "corresponding model",
        )
    end
    push!(tape, NodeData(MOIVARIABLE, x.index.value, parent))
    return nothing
end

function _parse_NL_expr_runtime(
    m::Model,
    x::NonlinearExpression,
    tape,
    parent,
    values,
)
    push!(tape, NodeData(SUBEXPRESSION, x.index, parent))
    return nothing
end

function _parse_NL_expr_runtime(
    m::Model,
    x::NonlinearParameter,
    tape,
    parent,
    values,
)
    push!(tape, NodeData(PARAMETER, x.index, parent))
    return nothing
end

function _parse_NL_expr_runtime(
    m::Model,
    x::AbstractArray,
    tape,
    parent,
    values,
)
    return error(
        "Unexpected array $x in nonlinear expression. Nonlinear expressions may contain only scalar expressions.",
    )
end

function _parse_NL_expr_runtime(
    m::Model,
    x::GenericQuadExpr,
    tape,
    parent,
    values,
)
    push!(tape, NodeData(CALL, operator_to_id[:+], parent))
    sum_parent = length(tape)
    _parse_NL_expr_runtime(m, x.aff, tape, sum_parent, values)
    for (xy, c) in x.terms
        push!(tape, NodeData(CALL, operator_to_id[:*], sum_parent))
        mult_parent = length(tape)
        _parse_NL_expr_runtime(m, xy.a, tape, mult_parent, values)
        _parse_NL_expr_runtime(m, xy.b, tape, mult_parent, values)
        if !isone(c)  # Optimization: no need for * node.
            _parse_NL_expr_runtime(m, c, tape, mult_parent, values)
        end
    end
    return
end

function _parse_NL_expr_runtime(
    m::Model,
    x::GenericAffExpr,
    tape,
    parent,
    values,
)
    push!(tape, NodeData(CALL, operator_to_id[:+], parent))
    sum_parent = length(tape)
    if !iszero(x.constant)
        _parse_NL_expr_runtime(m, x.constant, tape, sum_parent, values)
    end
    for (v, c) in x.terms
        if isone(c)  # Optimization: no need for * node.
            _parse_NL_expr_runtime(m, v, tape, sum_parent, values)
        else
            push!(tape, NodeData(CALL, operator_to_id[:*], sum_parent))
            mult_parent = length(tape)
            _parse_NL_expr_runtime(m, c, tape, mult_parent, values)
            _parse_NL_expr_runtime(m, v, tape, mult_parent, values)
        end
    end
    return
end

function _parse_NL_expr_runtime(m::Model, x, tape, parent, values)
    return error(
        "Unexpected object $x (of type $(typeof(x)) in nonlinear expression.",
    )
end

function _parse_NL_expr_runtime(m, x, tape, parent, values)
    return error(
        "Encountered an error parsing nonlinear expression: we don't support " *
        "models of type $(typeof(m)). In general, JuMP's nonlinear features " *
        "don't work with JuMP-extensions.",
    )
end

function _expression_complexity(ex::Expr)
    return isempty(ex.args) ? 1 : sum(_expression_complexity, ex.args)
end
_expression_complexity(other) = 1

# This is separated from the macro version to make it available for other @NL*
# macros.
function _process_NL_expr(model, ex)
    # This is an arbitrary cutoff. See issue #1355.
    if _expression_complexity(ex) > 5000
        @warn "Processing a very large nonlinear expression with " *
              "@NLexpression/@NLconstraint/@NLobjective. This may be very " *
              "slow. Consider using setNLobjective() and addNLconstraint() " *
              "instead of the macros or reformulating the expressions using " *
              "sum() and prod() to make them more compact. The macros are " *
              "designed to process smaller, human-readable expressions."
    end
    parsed = _parse_NL_expr(model, ex, :tape, -1, :values)
    return quote
        tape = NodeData[]
        values = Float64[]
        $parsed
        _NonlinearExprData(tape, values)
    end
end

macro _process_NL_expr(model, ex)
    return _process_NL_expr(model, ex)
end

function _Derivatives.expr_to_nodedata(
    ex::VariableRef,
    nd::Vector{NodeData},
    values::Vector{Float64},
    parentid,
    r::_Derivatives.UserOperatorRegistry,
)
    push!(nd, NodeData(MOIVARIABLE, ex.index.value, parentid))
    return nothing
end

function _Derivatives.expr_to_nodedata(
    ex::NonlinearExpression,
    nd::Vector{NodeData},
    values::Vector{Float64},
    parentid,
    r::_Derivatives.UserOperatorRegistry,
)
    push!(nd, NodeData(SUBEXPRESSION, ex.index, parentid))
    return nothing
end

function _Derivatives.expr_to_nodedata(
    ex::NonlinearParameter,
    nd::Vector{NodeData},
    values::Vector{Float64},
    parentid,
    r::_Derivatives.UserOperatorRegistry,
)
    push!(nd, NodeData(PARAMETER, ex.index, parentid))
    return nothing
end

function _Derivatives.expr_to_nodedata(
    ex::GenericAffExpr,
    nd::Vector{NodeData},
    values::Vector{Float64},
    parentid,
    r::_Derivatives.UserOperatorRegistry,
)
    push!(nd, NodeData(CALL, operator_to_id[:+], parentid))
    sum_parent = length(nd)
    if !iszero(ex.constant)
        _Derivatives.expr_to_nodedata(ex.constant, nd, values, sum_parent, r)
    end
    for (v, c) in ex.terms
        if isone(c)  # Optimization: no need for * node.
            _Derivatives.expr_to_nodedata(v, nd, values, sum_parent, r)
        else
            push!(nd, NodeData(CALL, operator_to_id[:*], sum_parent))
            mult_parent = length(nd)
            _Derivatives.expr_to_nodedata(c, nd, values, mult_parent, r)
            _Derivatives.expr_to_nodedata(v, nd, values, mult_parent, r)
        end
    end
    return
end

function _Derivatives.expr_to_nodedata(
    ex::GenericQuadExpr,
    nd::Vector{NodeData},
    values::Vector{Float64},
    parentid,
    r::_Derivatives.UserOperatorRegistry,
)
    push!(nd, NodeData(CALL, operator_to_id[:+], parentid))
    sum_parent = length(nd)
    _Derivatives.expr_to_nodedata(ex.aff, nd, values, sum_parent, r)
    for (xy, c) in ex.terms
        push!(nd, NodeData(CALL, operator_to_id[:*], sum_parent))
        mult_parent = length(nd)
        _Derivatives.expr_to_nodedata(xy.a, nd, values, mult_parent, r)
        _Derivatives.expr_to_nodedata(xy.b, nd, values, mult_parent, r)
        if !isone(c)  # Optimization: no need for * node.
            _Derivatives.expr_to_nodedata(c, nd, values, mult_parent, r)
        end
    end
    return
end

# Construct a _NonlinearExprData from a Julia expression.
# VariableRef objects should be spliced into the expression.
function _NonlinearExprData(m::Model, ex::Expr)
    _init_NLP(m)
    _check_expr(m, ex)
    nd, values = _Derivatives.expr_to_nodedata(ex, m.nlp_data.user_operators)
    return _NonlinearExprData(nd, values)
end
_NonlinearExprData(m::Model, ex) = _NonlinearExprData(m, :($ex + 0))

# Error if:
# 1) Unexpected expression
# 2) VariableRef doesn't match the model
function _check_expr(m::Model, ex::Expr)
    if ex.head == :ref # if we have x[1] already in there, something is wrong
        error(
            "Unrecognized expression $ex. JuMP variable objects and input coefficients should be spliced directly into expressions.",
        )
    end
    for e in ex.args
        _check_expr(m, e)
    end
    return
end
function _check_expr(m::Model, v::VariableRef)
    owner_model(v) === m || error("Variable $v does not belong to this model.")
    return
end
_check_expr(m::Model, ex) = nothing
