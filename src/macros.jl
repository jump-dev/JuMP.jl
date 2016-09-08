#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta

issum(s::Symbol) = (s == :sum) || (s == :∑) || (s == :Σ)
isprod(s::Symbol) = (s == :prod) || (s == :∏)

# for backward compatibility with 0.4
function comparison_to_call(ex)
    if isexpr(ex,:comparison) && length(ex.args) == 3
        return Expr(:call,ex.args[2],ex.args[1],ex.args[3])
    else
        return ex
    end
end

function flatten_error(ex)
    isexpr(ex, :flatten) || return
    error("The generator syntax \"$ex\" with multiple \"for\" statements is not yet supported. For JuMP, you should contract statements like \"for i in 1:N for j in 1:i\" into \"for i in 1:N, j in 1:i\"")
end


include("parseExpr_staged.jl")

###############################################################################
# buildrefsets
# Unexported. Takes as input an object representing a name, associated index
# sets, and conditions on those sets, for example
# buildrefsets(:(x[i=1:3,[:red,:blue]],k=S; i+k <= 6))
# Used internally in macros to build JuMPContainers and constraints. Returns
#       refcall:  Expr to reference a particular element, e.g. :(x[i,j,k])
#       idxvars:  Index names used in referencing, e.g.g {:i,:j,:k}
#       idxsets:  Index sets for indexing, e.g. {1:3, [:red,:blue], S}
#       idxpairs: Vector of IndexPair
#       condition: Expr containing any condition present for indexing
# Note in particular that it does not actually evaluate the condition, and so
# it returns just the cartesian product of possible indices.
function buildrefsets(expr::Expr, cname)
    c = copy(expr)
    idxvars = Any[]
    idxsets = Any[]
    idxpairs = IndexPair[]
    # Creating an indexed set of refs
    refcall = Expr(:ref, cname)
    if isexpr(c, :typed_vcat) || isexpr(c, :ref)
        shift!(c.args)
    end
    condition = :()
    if isexpr(c, :vcat) || isexpr(c, :typed_vcat)
        if isexpr(c.args[1], :parameters)
            @assert length(c.args[1].args) == 1
            condition = shift!(c.args).args[1]
        else
            condition = pop!(c.args)
        end
    end

    for s in c.args
        parse_done = false
        if isa(s, Expr)
            parse_done, idxvar, _idxset = tryParseIdxSet(s::Expr)
            if parse_done
                idxset = esc(_idxset)
                push!(idxpairs, IndexPair(idxvar, _idxset))
            end
        end
        if !parse_done # No index variable specified
            idxvar = gensym()
            idxset = esc(s)
            push!(idxpairs, IndexPair(nothing,s))
        end
        push!(idxvars, idxvar)
        push!(idxsets, idxset)
        push!(refcall.args, esc(idxvar))
    end
    return refcall, idxvars, idxsets, idxpairs, condition
end

buildrefsets(c, cname)  = (cname, Any[], Any[], IndexPair[], :())
buildrefsets(c) = buildrefsets(c, getname(c))

###############################################################################
# getloopedcode
# Unexported. Takes a bit of code and corresponding looping information and
# returns that code nested in corresponding loops, along with preceding code
# to construct an appropriate container. Input is:
#       c: symbolic representation of name and appropriate indexing sets, if
#          any. E.g. :(myvar) or :(x[i=1:3,[:red,:blue]])
#       code: inner loop code kernel to be nested in the loops
#       condition: a boolean expression to be evaluated before each kernel.
#                  If none, pass :().
#       idxvars: As defined for buildrefsets
#       idxsets: As defined for buildrefsets
#       idxpairs: As defined for buildrefsets
#       sym: A symbol or expression containing the element type of the
#            resulting container, e.g. :AffExpr or :Variable
function getloopedcode(varname, code, condition, idxvars, idxsets, idxpairs, sym; lowertri=false)

    # if we don't have indexing, just return to avoid allocating stuff
    if isempty(idxsets)
        return code
    end

    hascond = (condition != :())

    if lowertri
        @assert !hascond
        @assert length(idxvars)  == 2
        @assert length(idxpairs) == 2
        @assert !hasdependentsets(idxvars, idxsets)

        i, j = esc(idxvars[1]), esc(idxvars[2])
        expr = copy(code)
        vname = expr.args[1].args[1]
        expr.args[1] = :tmp
        code = quote
            let
                $(localvar(i))
                $(localvar(j))
                for $i in $(idxsets[1]), $j in $(idxsets[2])
                    $i <= $j || continue
                    $expr
                    $vname[$i,$j] = tmp
                    $vname[$j,$i] = tmp
                end
            end
        end
    else
        if hascond
            code = quote
                $(esc(condition)) || continue
                $code
            end
        end
        for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
            code = quote
                let
                    $(localvar(esc(idxvar)))
                    for $(esc(idxvar)) in $idxset
                        $code
                    end
                end
            end
        end
    end
    if hascond || hasdependentsets(idxvars,idxsets)
        # force a JuMPDict
        N = length(idxsets)
        mac = :($varname = JuMPDict{$(sym),$N}())
    else
        mac = gendict(varname, sym, idxsets...)
    end
    return quote
        $mac
        $code
        nothing
    end
end

localvar(x::Symbol) = _localvar(x)
localvar(x::Expr) = Expr(:block, _localvar(x)...)
_localvar(x::Symbol) = :(local $(esc(x)))
function _localvar(x::Expr)
    @assert x.head == :escape
    args = Any[]
    for t in x.args
        if isa(t, Symbol)
            push!(args, :(local $(esc(t))))
        else
            @assert isa(t, Expr)
            if t.head == :tuple
                append!(args, map(_localvar, t.args))
            else
                error("Internal error defining local variables in macros; please file an issue at https://github.com/JuliaOpt/JuMP.jl/issues/new")
            end
        end
    end
    args
end

getname(c::Symbol) = c
getname(c::Void) = ()
getname(c::AbstractString) = c
function getname(c::Expr)
    if c.head == :string
        return c
    else
        return c.args[1]
    end
end

validmodel(m::AbstractModel, name) = nothing
validmodel(m::MathProgBase.MathProgCallbackData, name) = error("Expected $name to be a JuMP model, but it is a callback object. Use of this macro is not supported within callbacks.")
validmodel(m, name) = error("Expected $name to be a JuMP model, but it has type ", typeof(m))

function assert_validmodel(m, macrocode)
    # assumes m is already escaped
    quote
        validmodel($m, $(quot(m.args[1])))
        $macrocode
    end
end

function _canonicalize_sense(sns::Symbol)
    if sns == :(==)
        return (:(==),false)
    elseif sns == :(>=) || sns == :(≥)
        return (:(>=),false)
    elseif sns == :(<=) || sns == :(≤)
        return (:(<=),false)
    elseif sns == :(.==)
        return (:(==),true)
    elseif sns == :(.>=) || sns == :(.≥)
        return (:(>=),true)
    elseif sns == :(.<=) || sns == :(.≤)
        return (:(<=),true)
    else
        error("Unrecognized sense $sns")
    end
end

function _construct_constraint!(args...)
    warn("_construct_constraint! is deprecated. Use constructconstraint! instead")
    constructconstraint(args...)
end

# two-argument constructconstraint! is used for one-sided constraints.
# Right-hand side is zero.
constructconstraint!(v::Variable, sense::Symbol) = constructconstraint(convert(AffExpr,v), sense)
function constructconstraint!(aff::AffExpr, sense::Symbol)
    offset = aff.constant
    aff.constant = 0.0
    if sense == :(<=) || sense == :≤
        return LinearConstraint(aff, -Inf, -offset)
    elseif sense == :(>=) || sense == :≥
        return LinearConstraint(aff, -offset, Inf)
    elseif sense == :(==)
        return LinearConstraint(aff, -offset, -offset)
    else
        error("Cannot handle ranged constraint")
    end
end

function constructconstraint!(aff::AffExpr, lb, ub)
    offset = aff.constant
    aff.constant = 0.0
    LinearConstraint(aff, lb-offset, ub-offset)
end

constructconstraint!(quad::QuadExpr, sense::Symbol) = QuadConstraint(quad, sense)

function constructconstraint!(normexpr::SOCExpr, sense::Symbol)
    # check that the constraint is SOC representable
    if sense == :(<=)
        SOCConstraint(normexpr)
    elseif sense == :(>=)
        SOCConstraint(-normexpr)
    else
        error("Invalid sense $sense in SOC constraint")
    end
end

constructconstraint!(x::Array, sense::Symbol) = map(c->constructconstraint!(c,sense), x)

_vectorize_like(x::Number, y::Array{AffExpr}) = fill(x, size(y))
function _vectorize_like{R<:Number}(x::Array{R}, y::Array{AffExpr})
    for i in 1:max(ndims(x),ndims(y))
        size(x,i) == size(y,i) || error("Unequal sizes for ranged constraint")
    end
    x
end

function constructconstraint!(x::Array{AffExpr}, lb, ub)
    LB = _vectorize_like(lb,x)
    UB = _vectorize_like(ub,x)
    ret = similar(x, LinearConstraint)
    map!(ret, 1:length(ret)) do i
        constructconstraint!(x[i], LB[i], UB[i])
    end
end

# three-argument constructconstraint! is used for two-sided constraints.
function constructconstraint!(aff::AffExpr, lb::Real, ub::Real)
    offset = aff.constant
    aff.constant = 0.0
    LinearConstraint(aff,lb-offset,ub-offset)
end

constructconstraint!(aff::Variable, lb::Real, ub::Real) = constructconstraint!(convert(AffExpr,v),lb,ub)

constructconstraint!(q::QuadExpr, lb, ub) = error("Two-sided quadratic constraints not supported. (Try @NLconstraint instead.)")

constraint_error(args, str) = error("In @constraint($(join(args,","))): ", str)

macro constraint(args...)
    # Pick out keyword arguments
    if isexpr(args[1],:parameters) # these come if using a semicolon
        kwargs = args[1]
        args = args[2:end]
    else
        kwargs = Expr(:parameters)
    end
    append!(kwargs.args, collect(filter(x -> isexpr(x, :kw), args))) # comma separated
    args = collect(filter(x->!isexpr(x, :kw), args))

    if length(args) < 2
        if length(kwargs.args) > 0
            constraint_error(args, "Not enough positional arguments")
        else
            constraint_error(args, "Not enough arguments")
        end
    end
    m = args[1]
    x = args[2]
    extra = args[3:end]

    m = esc(m)
    # Two formats:
    # - @constraint(m, a*x <= 5)
    # - @constraint(m, myref[a=1:5], a*x <= 5)
    length(extra) > 1 && constraint_error(args, "Too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : gensym()
    x = length(extra) == 1 ? extra[1] : x

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    variable = gensym()
    quotvarname = quot(getname(c))
    escvarname  = anonvar ? variable : esc(getname(c))

    if isa(x, Symbol)
        constraint_error(args, "Incomplete constraint specification $x. Are you missing a comparison (<=, >=, or ==)?")
    end

    (x.head == :block) &&
        constraint_error(args, "Code block passed as constraint. Perhaps you meant to use @constraints instead?")
    if VERSION < v"0.5.0-dev+3231"
        x = comparison_to_call(x)
    end

    # Strategy: build up the code for non-macro addconstraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c, variable)
    # Build the constraint
    if isexpr(x, :call)
        # Simple comparison - move everything to the LHS
        @assert length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[1])
        lhs = :($(x.args[2]) - $(x.args[3]))
        addconstr = (vectorized ? :addVectorizedConstraint : :addconstraint)
        newaff, parsecode = parseExprToplevel(lhs, :q)
        constraintcall = :($addconstr($m, constructconstraint!($newaff,$(quot(sense)))))
        for kw in kwargs.args
            @assert isexpr(kw, :kw)
            push!(constraintcall.args, esc(kw))
        end
        code = quote
            q = zero(AffExpr)
            $parsecode
            $(refcall) = $constraintcall
        end
    elseif isexpr(x, :comparison)
        # Ranged row
        (lsign,lvectorized) = _canonicalize_sense(x.args[2])
        (rsign,rvectorized) = _canonicalize_sense(x.args[4])
        if (lsign != :(<=)) || (rsign != :(<=))
            constraint_error(args, "Only ranged rows of the form lb <= expr <= ub are supported.")
        end
        ((vectorized = lvectorized) == rvectorized) || constraint_error("Signs are inconsistently vectorized")
        addconstr = (lvectorized ? :addVectorizedConstraint : :addconstraint)
        x_str = string(x)
        lb_str = string(x.args[1])
        ub_str = string(x.args[5])
        newaff, parsecode = parseExprToplevel(x.args[3],:aff)

        newlb, parselb = parseExprToplevel(x.args[1],:lb)
        newub, parseub = parseExprToplevel(x.args[5],:ub)

        constraintcall = :($addconstr($m, constructconstraint!($newaff,$newlb,$newub)))
        for kw in kwargs.args
            @assert isexpr(kw, :kw)
            push!(constraintcall.args, esc(kw))
        end
        code = quote
            aff = zero(AffExpr)
            $parsecode
            lb = 0.0
            $parselb
            ub = 0.0
            $parseub
        end
        if vectorized
            code = quote
                $code
                lbval, ubval = $newlb, $newub
            end
        else
            code = quote
                $code
                CoefType = coeftype($newaff)
                try
                    lbval = convert(CoefType, $newlb)
                catch
                    constraint_error(args, string("Expected ",$lb_str," to be a ", CoefType, "."))
                end
                try
                    ubval = convert(CoefType, $newub)
                catch
                    constraint_error(args, string("Expected ",$ub_str," to be a ", CoefType, "."))
                end
            end
        end
        code = quote
            $code
            $(refcall) = $constraintcall
        end
    else
        # Unknown
        constraint_error(args, string("Constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2\n" * "       lb <= expr <= ub"))
    end
    return assert_validmodel(m, quote
        $(getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :ConstraintRef))
        registercon($m, $quotvarname, $variable)
        $(anonvar ? variable : :($escvarname = $variable))
    end)
end

macro SDconstraint(m, x)
    m = esc(m)

    if isa(x, Symbol)
        error("in @SDConstraint: Incomplete constraint specification $x. Are you missing a comparison (<= or >=)?")
    end

    (x.head == :block) &&
        error("Code block passed as constraint.")
    if VERSION < v"0.5.0-dev+3231"
        x = comparison_to_call(x)
    end
    isexpr(x,:call) && length(x.args) == 3 || error("in @SDconstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2")
    # Build the constraint
    # Simple comparison - move everything to the LHS
    sense = x.args[1]
    if sense == :⪰
        sense = :(>=)
    elseif sense == :⪯
        sense = :(<=)
    end
    sense,_ = _canonicalize_sense(sense)
    lhs = :()
    if sense == :(>=)
        lhs = :($(x.args[2]) - $(x.args[3]))
    elseif sense == :(<=)
        lhs = :($(x.args[3]) - $(x.args[2]))
    else
        error("Invalid sense $sense in SDP constraint")
    end
    newaff, parsecode = parseExprToplevel(lhs, :q)
    assert_validmodel(m, quote
        q = zero(AffExpr)
        $parsecode
        addconstraint($m, SDConstraint($newaff))
    end)
end

macro LinearConstraint(x)
    (x.head == :block) &&
        error("Code block passed as constraint. Perhaps you meant to use @LinearConstraints instead?")

    if VERSION < v"0.5.0-dev+3231"
        x = comparison_to_call(x)
    end
    if isexpr(x, :call) && length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[1])
        # Simple comparison - move everything to the LHS
        vectorized &&
            error("in @LinearConstraint ($(string(x))): Cannot add vectorized constraints")
        lhs = :($(x.args[2]) - $(x.args[3]))
        return quote
            newaff = @Expression($(esc(lhs)))
            c = constructconstraint!(newaff,$(quot(sense)))
            isa(c, LinearConstraint) ||
                error("Constraint in @LinearConstraint is really a $(typeof(c))")
            c
        end
    elseif isexpr(x, :comparison)
        # Ranged row
        (lsense,lvectorized) = _canonicalize_sense(x.args[2])
        (rsense,rvectorized) = _canonicalize_sense(x.args[4])
        if (lsense != :<=) || (rsense != :<=)
            error("in @constraint ($(string(x))): only ranged rows of the form lb <= expr <= ub are supported.")
        end
        (lvectorized || rvectorized) &&
            error("in @LinearConstraint ($(string(x))): Cannot add vectorized constraints")
        lb = x.args[1]
        ub = x.args[5]
        return quote
            if !isa($(esc(lb)),Number)
                error(string("in @LinearConstraint (",$(string(x)),"): expected ",$(string(lb))," to be a number."))
            elseif !isa($(esc(ub)),Number)
                error(string("in @LinearConstraint (",$(string(x)),"): expected ",$(string(ub))," to be a number."))
            end
            newaff = @Expression($(esc(x.args[3])))
            offset = newaff.constant
            newaff.constant = 0.0
            isa(newaff,AffExpr) || error("Ranged quadratic constraints are not allowed")
            LinearConstraint(newaff,$(esc(lb))-offset,$(esc(ub))-offset)
        end
    else
        # Unknown
        error("in @LinearConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2\n" * "       lb <= expr <= ub")
    end
end

macro QuadConstraint(x)
    (x.head == :block) &&
        error("Code block passed as constraint. Perhaps you meant to use @QuadConstraints instead?")

    if VERSION < v"0.5.0-dev+3231"
        x = comparison_to_call(x)
    end
    if isexpr(x, :call) && length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[1])
        # Simple comparison - move everything to the LHS
        vectorized &&
            error("in @QuadConstraint ($(string(x))): Cannot add vectorized constraints")
        lhs = :($(x.args[2]) - $(x.args[3]))
        return quote
            newaff = @Expression($(esc(lhs)))
            q = constructconstraint!(newaff,$(quot(sense)))
            isa(q, QuadConstraint) || error("Constraint in @QuadConstraint is really a $(typeof(q))")
            q
        end
    elseif isexpr(x, :comparison)
        error("Ranged quadratic constraints are not allowed")
    else
        # Unknown
        error("in @QuadConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2")
    end
end

macro SOCConstraint(x)
    (x.head == :block) &&
        error("Code block passed as constraint. Perhaps you meant to use @SOCConstraints instead?")

    if VERSION < v"0.5.0-dev+3231"
        x = comparison_to_call(x)
    end
    if isexpr(x, :call) && length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[1])
        # Simple comparison - move everything to the LHS
        vectorized &&
            error("in @SOCConstraint ($(string(x))): Cannot add vectorized constraints")
        lhs = :($(x.args[2]) - $(x.args[3]))
        return quote
            newaff = @Expression($(esc(lhs)))
            q = constructconstraint!(newaff,$(quot(sense)))
            isa(q, SOCConstraint) || error("Constraint in @SOCConstraint is really a $(typeof(q))")
            q
        end
    elseif isexpr(x, :comparison)
        error("Ranged second-order cone constraints are not allowed")
    else
        # Unknown
        error("in @SOCConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2")
    end
end

for (mac,sym) in [(:LinearConstraints, Symbol("@LinearConstraint")),
                  (:QuadConstraints,   Symbol("@QuadConstraint")),
                  (:SOCConstraints,    Symbol("@SOCConstraint"))]
    @eval begin
        macro $mac(x)
            x.head == :block || error(string("Invalid syntax for ", $(quot(sym))))
            @assert x.args[1].head == :line
            code = Expr(:vect)
            for it in x.args
                if it.head == :line
                    # do nothing
                elseif it.head == :comparison # regular constraint
                    push!(code.args, Expr(:macrocall, $sym, esc(it)))
                elseif it.head == :tuple # constraint ref
                    if all([isexpr(arg,:comparison) for arg in it.args]...)
                        # the user probably had trailing commas at end of lines, e.g.
                        # @LinearConstraints(m, begin
                        #     x <= 1,
                        #     x >= 1
                        # end)
                        error(string("Invalid syntax in ", $(quot(sym)), ". Do you have commas at the end of a line specifying a constraint?"))
                    end
                    error(string($(quot(sym)), " does not currently support the two argument syntax for specifying groups of constraints in one line."))
                else
                    error("Unexpected constraint expression $it")
                end
            end
            return code
        end
    end
end

for (mac,sym) in [(:constraints,  Symbol("@constraint")),
                  (:NLconstraints,Symbol("@NLconstraint")),
                  (:SDconstraints,Symbol("@SDconstraint")),
                  (:variables,Symbol("@variable")),
                  (:expressions, Symbol("@expression")),
                  (:NLexpressions, Symbol("@NLexpression"))]
    @eval begin
        macro $mac(m, x)
            x.head == :block || error("Invalid syntax for @$mac")
            @assert x.args[1].head == :line
            code = quote end
            for it in x.args
                if isexpr(it, :line)
                    # do nothing
                elseif isexpr(it, :tuple) # line with commas
                    args = []
                    for ex in it.args
                        if isexpr(ex, :tuple) # embedded tuple
                            append!(args, ex.args)
                        else
                            push!(args, ex)
                        end
                    end
                    args_esc = []
                    for ex in args
                        if isexpr(ex, :(=))
                            push!(args_esc,Expr(:kw, ex.args[1], esc(ex.args[2])))
                        else
                            push!(args_esc, esc(ex))
                        end
                    end
                    mac = Expr(:macrocall,$(quot(sym)), esc(m), args_esc...)
                    push!(code.args, mac)
                else # stand-alone symbol or expression
                    push!(code.args,Expr(:macrocall,$(quot(sym)), esc(m), esc(it)))
                end
            end
            push!(code.args, :(nothing))
            return code
        end
    end
end

macro objective(m, args...)
    m = esc(m)
    if length(args) != 2
        # Either just an objective sene, or just an expression.
        error("in @objective: needs three arguments: model, objective sense (Max or Min) and expression.")
    end
    sense, x = args
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    newaff, parsecode = parseExprToplevel(x, :q)
    code = quote
        q = zero(AffExpr)
        $parsecode
        setobjective($m, $(esc(sense)), $newaff)
    end
    return assert_validmodel(m, code)
end

# Return a standalone, unnamed expression
# ex = @Expression(2x + 3y)
# Currently for internal use only.
macro Expression(x)
    newaff, parsecode = parseExprToplevel(x, :q)
    return quote
        q = 0.0
        $parsecode
        $newaff
    end
end


macro expression(args...)
    if length(args) == 3
        m = esc(args[1])
        c = args[2]
        x = args[3]
    elseif length(args) == 1
        m = nothing
        c = gensym()
        x = args[1]
        Base.warn_once("The one-argument version of @defExpr is deprecated. The corresponding JuMP model is now required as the first argument, and a name for the expression or collection of expressions is required as the second argument. The new syntax is @expression(<JuMP model>, <name of expression(s)>, <expression>)")
    elseif length(args) == 2
        m = nothing
        c = args[1]
        x = args[2]
        Base.warn_once("The two-argument version of @defExpr is deprecated. The corresponding JuMP model is now required as the first argument. The new syntax is @expression(<JuMP model>, <name of expression(s)>, <expression>)")
    else
        error("@expression: needs three arguments.")
    end

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    variable = gensym()
    escvarname  = anonvar ? variable : esc(getname(c))

    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c, variable)
    newaff, parsecode = parseExprToplevel(x, :q)
    code = quote
        q = 0.0
        $parsecode
    end
    if isa(c,Expr)
        code = quote
            $code
            (isa($newaff,AffExpr) || isa($newaff,Number) || isa($newaff,Variable)) || error("Collection of expressions with @expression must be linear. For quadratic expressions, use your own array.")
        end
    end
    code = quote
        $code
        $(refcall) = $newaff
    end
    code = getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :AffExpr)
    if m === nothing # deprecated usage
        return quote
            $code
            $(anonvar ? variable : :($escvarname = $variable))
        end
    else
        # don't do anything with the model, but check that it's valid anyway
        return assert_validmodel(m, quote
            $code
            $(anonvar ? variable : :($escvarname = $variable))
        end)
    end
end

function hasdependentsets(idxvars, idxsets)
    # check if any index set depends on a previous index var
    for i in 2:length(idxsets)
        for v in idxvars[1:(i-1)]
            if dependson(idxsets[i],v)
                return true
            end
        end
    end
    return false
end

dependson(ex::Expr,s::Symbol) = any(a->dependson(a,s), ex.args)
dependson(ex::Symbol,s::Symbol) = (ex == s)
dependson(ex,s::Symbol) = false
function dependson(ex1,ex2)
    @assert isa(ex2, Expr)
    @assert ex2.head == :tuple
    any(s->dependson(ex1,s), ex2.args)
end

function isdependent(idxvars,idxset,i)
    for (it,idx) in enumerate(idxvars)
        it == i && continue
        dependson(idxset, idx) && return true
    end
    return false
end

esc_nonconstant(x::Number) = x
esc_nonconstant(x) = esc(x)

const EMPTYSTRING = UTF8String("")

variable_error(args, str) = error("In @variable($(join(args,","))): ", str)

macro variable(args...)
    m = esc(args[1])

    extra = vcat(args[2:end]...)
    # separate out keyword arguments
    kwargs = filter(ex->isexpr(ex,:kw), extra)
    extra = filter(ex->!isexpr(ex,:kw), extra)

    # if there is only a single non-keyword argument, this is an anonymous
    # variable spec and the one non-kwarg is the model
    if length(kwargs) == length(args)-1
        x = gensym()
        anon_singleton = true
    else
        x = shift!(extra)
        if x in [:Cont,:Int,:Bin,:SemiCont,:SemiInt,:SDP]
            variable_error(args, "Ambiguous variable name $x detected. This conflicts with a possible variable type. Rename the variable.")
        end
        anon_singleton = false
    end

    t = :Cont
    gottype = false
    haslb = false
    hasub = false
    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    if VERSION < v"0.5.0-dev+3231"
        x = comparison_to_call(x)
    end
    explicit_comparison = false
    if isexpr(x,:comparison) # two-sided
        explicit_comparison = true
        haslb = true
        hasub = true
        if x.args[2] == :>= || x.args[2] == :≥
            # ub >= x >= lb
            x.args[4] == :>= || x.args[4] == :≥ || variable_error(args, "Invalid variable bounds")
            var = x.args[3]
            lb = esc_nonconstant(x.args[5])
            ub = esc_nonconstant(x.args[1])
        elseif x.args[2] == :<= || x.args[2] == :≤
            # lb <= x <= u
            var = x.args[3]
            (x.args[4] != :<= && x.args[4] != :≤) &&
                variable_error(args, "Expected <= operator after variable name.")
            lb = esc_nonconstant(x.args[1])
            ub = esc_nonconstant(x.args[5])
        else
            variable_error(args, "Use the form lb <= ... <= ub.")
        end
    elseif isexpr(x,:call)
        explicit_comparison = true
        if x.args[1] == :>= || x.args[1] == :≥
            # x >= lb
            var = x.args[2]
            @assert length(x.args) == 3
            lb = esc_nonconstant(x.args[3])
            haslb = true
            ub = Inf
        elseif x.args[1] == :<= || x.args[1] == :≤
            # x <= ub
            var = x.args[2]
            # NB: May also be lb <= x, which we do not support
            #     We handle this later in the macro
            @assert length(x.args) == 3
            ub = esc_nonconstant(x.args[3])
            hasub = true
            lb = -Inf
        elseif x.args[1] == :(==)
            # fixed variable
            var = x.args[2]
            @assert length(x.args) == 3
            lb = esc(x.args[3])
            haslb = true
            ub = esc(x.args[3])
            hasub = true
            gottype = true
            t = :Fixed
        else
            # Its a comparsion, but not using <= ... <=
            variable_error(args, "Unexpected syntax $(string(x)).")
        end
    else
        # No bounds provided - free variable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        var = x
        lb = -Inf
        ub = Inf
    end

    anonvar = isexpr(var, :vect) || isexpr(var, :vcat) || anon_singleton
    anonvar && explicit_comparison && error("Cannot use explicit bounds via >=, <= with an anonymous variable")
    variable = gensym()
    quotvarname = anonvar ? :(:__anon__) : quot(getname(var))
    escvarname  = anonvar ? variable     : esc(getname(var))

    # process keyword arguments
    value = NaN
    obj = nothing
    inconstraints = nothing
    coefficients = nothing
    for ex in kwargs
        kwarg = ex.args[1]
        if kwarg == :start
            value = esc(ex.args[2])
        elseif kwarg == :objective
            obj = esc(ex.args[2])
        elseif kwarg == :inconstraints
            inconstraints = esc(ex.args[2])
        elseif kwarg == :coefficients
            coefficients = esc(ex.args[2])
        elseif kwarg == :basename
            quotvarname = esc(ex.args[2])
        elseif kwarg == :lowerbound
            haslb && variable_error(args, "Cannot specify variable lowerbound twice")
            lb = esc_nonconstant(ex.args[2])
            haslb = true
        elseif kwarg == :upperbound
            hasub && variable_error(args, "Cannot specify variable upperbound twice")
            ub = esc_nonconstant(ex.args[2])
            hasub = true
        else
            variable_error(args, "Unrecognized keyword argument $kwarg")
        end
    end

    if (obj !== nothing || inconstraints !== nothing || coefficients !== nothing) &&
       (obj === nothing || inconstraints === nothing || coefficients === nothing)
        variable_error(args, "Must provide 'objective', 'inconstraints', and 'coefficients' arguments all together for column-wise modeling")
    end

    sdp = any(t -> (t == :SDP), extra)
    symmetric = (sdp || any(t -> (t == :Symmetric), extra))
    extra = filter(x -> (x != :SDP && x != :Symmetric), extra) # filter out SDP and sym tag

    # Determine variable type (if present).
    # Types: default is continuous (reals)
    if length(extra) > 0
        if t == :Fixed
            variable_error(args, "Unexpected extra arguments when declaring a fixed variable")
        end
        if extra[1] in [:Bin, :Int, :SemiCont, :SemiInt]
            gottype = true
            t = extra[1]
        end

        if t == :Bin
            if (lb != -Inf || ub != Inf) && !(lb == 0.0 && ub == 1.0)
            variable_error(args, "Bounds other than [0, 1] may not be specified for binary variables.\nThese are always taken to have a lower bound of 0 and upper bound of 1.")
            else
                lb = 0.0
                ub = 1.0
            end
        end

        !gottype && variable_error(args, "Syntax error")
    end

    # Handle the column generation functionality
    if coefficients !== nothing
        !isa(var,Symbol) &&
        variable_error(args, "Can only create one variable at a time when adding to existing constraints.")

        return assert_validmodel(m, quote
            $variable = Variable($m,$lb,$ub,$(quot(t)),$obj,$inconstraints,$coefficients,UTF8String(string($quotvarname)),$value)
            $(anonvar ? variable : :($escvarname = $variable))
        end)
    end

    if isa(var,Symbol)
        # Easy case - a single variable
        sdp && variable_error(args, "Cannot add a semidefinite scalar variable")
        code = :($variable = Variable($m,$lb,$ub,$(quot(t)),UTF8String(string($quotvarname)),$value))
        if !anonvar
            code = quote
                $code
                registervar($m, $quotvarname, $variable)
                $escvarname = $variable
            end
        end
        return assert_validmodel(m, code)
    end
    isa(var,Expr) || variable_error(args, "Expected $var to be a variable name")

    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(var, variable)
    clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? () : idxsets[i])

    code = :( $(refcall) = Variable($m, $lb, $ub, $(quot(t)), EMPTYSTRING, $value) )
    if symmetric
        # Sanity checks on SDP input stuff
        condition == :() ||
            variable_error(args, "Cannot have conditional indexing for SDP variables")
        length(idxvars) == length(idxsets) == 2 ||
            variable_error(args, "SDP variables must be 2-dimensional")
        !symmetric || (length(idxvars) == length(idxsets) == 2) ||
            variable_error(args, "Symmetric variables must be 2-dimensional")
        hasdependentsets(idxvars, idxsets) &&
            variable_error(args, "Cannot have index dependencies in symmetric variables")
        for _rng in idxsets
            isexpr(_rng, :escape) ||
                variable_error(args, "Internal error 1")
            rng = _rng.args[1] # undo escaping
            (isexpr(rng,:(:)) && rng.args[1] == 1 && length(rng.args) == 2) ||
                variable_error(args, "Index sets for SDP variables must be ranges of the form 1:N")
        end

        if !(lb == -Inf && ub == Inf)
            variable_error(args, "Semidefinite or symmetric variables cannot be provided bounds")
        end
        return assert_validmodel(m, quote
            $(esc(idxsets[1].args[1].args[2])) == $(esc(idxsets[2].args[1].args[2])) || error("Cannot construct symmetric variables with nonsquare dimensions")
            (Compat.issymmetric($lb) && Compat.issymmetric($ub)) || error("Bounds on symmetric  variables must be symmetric")
            $(getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :Variable; lowertri=symmetric))
            $(if sdp
                quote
                    push!($(m).varCones, (:SDP, first($variable).col : last($variable).col))
                end
            end)
            push!($(m).dictList, $variable)
            registervar($m, $quotvarname, $variable)
            storecontainerdata($m, $variable, $quotvarname,
                               $(Expr(:tuple,idxsets...)),
                               $idxpairs, $(quot(condition)))
            $(anonvar ? variable : :($escvarname = $variable))
        end)
    else
        coloncheckcode = Expr(:call,:coloncheck,refcall.args[2:end]...)
        code = :($coloncheckcode; $code)
        return assert_validmodel(m, quote
            $(getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :Variable))
            isa($variable, JuMPContainer) && pushmeta!($variable, :model, $m)
            push!($(m).dictList, $variable)
            registervar($m, $quotvarname, $variable)
            storecontainerdata($m, $variable, $quotvarname,
                               $(Expr(:tuple,map(clear_dependencies,1:length(idxsets))...)),
                               $idxpairs, $(quot(condition)))
            $(anonvar ? variable : :($escvarname = $variable))
        end)
    end
end

storecontainerdata(m::Model, variable, varname, idxsets, idxpairs, condition) =
    m.varData[variable] = JuMPContainerData(varname, map(collect,idxsets), idxpairs, condition)

macro constraintref(var)
    if isa(var,Symbol)
        # easy case
        return esc(:(local $var))
    else
        if !isexpr(var,:ref)
            error("Syntax error: Expected $var to be of form var[...]")
        end

        varname = var.args[1]
        idxsets = var.args[2:end]
        idxpairs = IndexPair[]

        code = quote
            $(esc(gendict(varname, :ConstraintRef, idxsets...)))
            nothing
        end
        return code
    end
end

macro NLobjective(m, sense, x)
    m = esc(m)
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    return assert_validmodel(m, quote
        initNLP($m)
        setobjectivesense($m, $(esc(sense)))
        ex = @processNLExpr($m, $(esc(x)))
        $m.nlpdata.nlobj = ex
        $m.obj = zero(QuadExpr)
        $m.internalModelLoaded = false
        nothing
    end)
end

macro NLconstraint(m, x, extra...)
    m = esc(m)
    # Two formats:
    # - @NLconstraint(m, a*x <= 5)
    # - @NLconstraint(m, myref[a=1:5], sin(x^a) <= 5)
    length(extra) > 1 && error("in @NLconstraint: too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : gensym()
    x = length(extra) == 1 ? extra[1] : x

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    variable = gensym()
    quotvarname = anonvar ? :(:__anon__) : quot(getname(c))
    escvarname  = anonvar ? variable : esc(getname(c))

    if VERSION < v"0.5.0-dev+3231"
        x = comparison_to_call(x)
    end

    # Strategy: build up the code for non-macro addconstraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c, variable)
    # Build the constraint
    if isexpr(x, :call) # one-sided constraint
        # Simple comparison - move everything to the LHS
        op = x.args[1]
        if op == :(==)
            lb = 0.0
            ub = 0.0
        elseif op == :(<=) || op == :(≤)
            lb = -Inf
            ub = 0.0
        elseif op == :(>=) || op == :(≥)
            lb = 0.0
            ub = Inf
        else
            error("in @NLconstraint ($(string(x))): expected comparison operator (<=, >=, or ==).")
        end
        lhs = :($(x.args[2]) - $(x.args[3]))
        code = quote
            c = NonlinearConstraint(@processNLExpr($m, $(esc(lhs))), $lb, $ub)
            push!($m.nlpdata.nlconstr, c)
            $(refcall) = ConstraintRef{Model,NonlinearConstraint}($m, length($m.nlpdata.nlconstr))
        end
    elseif isexpr(x, :comparison)
        # ranged row
        if (x.args[2] != :<= && x.args[2] != :≤) || (x.args[4] != :<= && x.args[4] != :≤)
            error("in @NLconstraint ($(string(x))): only ranged rows of the form lb <= expr <= ub are supported.")
        end
        lb = x.args[1]
        ub = x.args[5]
        code = quote
            if !isa($(esc(lb)),Number)
                error(string("in @NLconstraint (",$(string(x)),"): expected ",$(string(lb))," to be a number."))
            elseif !isa($(esc(ub)),Number)
                error(string("in @NLconstraint (",$(string(x)),"): expected ",$(string(ub))," to be a number."))
            end
            c = NonlinearConstraint(@processNLExpr($m, $(esc(x.args[3]))), $(esc(lb)), $(esc(ub)))
            push!($m.nlpdata.nlconstr, c)
            $(refcall) = ConstraintRef{Model,NonlinearConstraint}($m, length($m.nlpdata.nlconstr))
        end
    else
        # Unknown
        error("in @NLconstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2")
    end
    looped = getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :(ConstraintRef{Model,NonlinearConstraint}))
    return assert_validmodel(m, quote
        initNLP($m)
        $m.internalModelLoaded = false
        $looped
        registercon($m, $quotvarname, $variable)
        $(anonvar ? variable : :($escvarname = $variable))
    end)
end

macro NLexpression(args...)
    if length(args) <= 2
        s = IOBuffer()
        print(s,args[1])
        if length(args) == 2
            print(s,",")
            print(s,args[2])
        end
        msg = """
        in @NLexpression($(takebuf_string(s))): three arguments are required.
        Note that the syntax of @NLexpression has recently changed:
        The first argument should be the model to which the expression is attached.
        The second is the name of the expression (or collection of expressions).
        The third is the expression itself.
        Example:
        @NLexpression(m, my_expr, x^2/y)
        @NLexpression(m, my_expr_collection[i=1:2], sin(z[i])^2)
        Support for the old syntax (with the model omitted) will be removed in an upcoming release.
        """
        Base.warn(msg)
        m = :(__last_model[1])
        if length(args) == 2
            c = args[1]
            x = args[2]
        else
            c = gensym()
            x = args[1]
        end
    else
        @assert length(args) == 3
        m, c, x = args
        m = esc(m)
    end

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    variable = gensym()
    escvarname  = anonvar ? variable : esc(getname(c))

    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c, variable)
    code = quote
        $(refcall) = NonlinearExpression($m, @processNLExpr($m, $(esc(x))))
    end
    return assert_validmodel(m, quote
        $(getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :NonlinearExpression))
        $(anonvar ? variable : :($escvarname = $variable))
    end)
end

# syntax is @NLparameter(m, p[i=1] == 2i)
macro NLparameter(m, ex)
    m = esc(m)
    if VERSION < v"0.5.0-dev+3231"
        ex = comparison_to_call(ex)
    end
    @assert isexpr(ex, :call)
    @assert length(ex.args) == 3
    @assert ex.args[1] == :(==)
    c = ex.args[2]
    x = ex.args[3]

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    if anonvar
        error("In @NLparameter($m, $ex): Anonymous nonlinear parameter syntax is not currently supported")
    end
    variable = gensym()
    escvarname  = anonvar ? variable : esc(getname(c))

    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c, variable)
    code = quote
        $(refcall) = newparameter($m, $(esc(x)))
    end
    return assert_validmodel(m, quote
        $(getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :NonlinearParameter))
        $(anonvar ? variable : :($escvarname = $variable))
    end)
end
