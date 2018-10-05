#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta

issum(s::Symbol) = (s == :sum) || (s == :∑) || (s == :Σ)
isprod(s::Symbol) = (s == :prod) || (s == :∏)

function curly_to_generator(x)
    # we have a filter condition
    x = copy(x)
    @assert isexpr(x,:curly)
    if isexpr(x.args[2],:parameters)
        cond = x.args[2].args[1]
        body = x.args[3]
        if length(x.args) == 3 # no iteration set!
            push!(x.args,:(_ in 1))
        end
        for i in length(x.args):-1:4
            if i == length(x.args)
                body = Expr(:generator,body,Expr(:filter,cond,x.args[i]))
            else
                body = Expr(:generator,body,x.args[i])
            end
        end
    else
        cond = nothing
        body = x.args[2]
        if length(x.args) == 2 # no iteration set!
            push!(x.args,:(_ in 1))
        end
        for i in length(x.args):-1:3
            body = Expr(:generator,body,x.args[i])
        end
    end
    if isexpr(body.args[1],:generator)
        body = Expr(:flatten,body)
    end
    name = x.args[1]
    if name == :norm2
        return Expr(:call,:norm,body)
    elseif name == :norm1
        return Expr(:call,:norm,body,1)
    elseif name == :norminf
        return Expr(:call,:norm,body,Inf)
    elseif name == :norm∞
        return Expr(:call,:norm,body,Inf)
    else
        return Expr(:call,name,body)
    end
end

function warn_curly(x)
    genform = curly_to_generator(x)
    if length(genform.args) == 2
        # don't print extra parens
        genstr = "$(genform.args[1])$(genform.args[2])"
    else
        genstr = "$genform"
    end
    warn_once("The curly syntax (sum{},prod{},norm2{}) is deprecated in favor of the new generator syntax (sum(),prod(),norm()).")
    warn_once("Replace $x with $genstr.")
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
    @static if VERSION >= v"0.7-"
        # On 0.7, :(t[i;j]) is a :ref, while t[i,j;j] is a :typed_vcat.
        # In both cases :t is the first arg.
        if isexpr(c, :typed_vcat) || isexpr(c, :ref)
            popfirst!(c.args)
        end
        condition = :()
        if isexpr(c, :vcat) || isexpr(c, :typed_vcat)
            # Parameters appear as plain args at the end.
            if length(c.args) > 2
                error("Unsupported syntax $c.")
            elseif length(c.args) == 2
                condition = pop!(c.args)
            end # else no condition.
        elseif isexpr(c, :ref) || isexpr(c, :vect)
            # Parameters appear at the front.
            if isexpr(c.args[1], :parameters)
                if length(c.args[1].args) != 1
                    error("Invalid syntax: $c. Multiple semicolons are not " *
                          "supported.")
                end
                condition = popfirst!(c.args).args[1]
            end
        end
        if isexpr(c, :vcat) || isexpr(c, :typed_vcat) || isexpr(c, :ref)
            if isexpr(c.args[1], :parameters)
                @assert length(c.args[1].args) == 1
                condition = popfirst!(c.args).args[1]
            end # else no condition.
        end
    else
        if isexpr(c, :typed_vcat) || isexpr(c, :ref)
            popfirst!(c.args)
        end
        condition = :()
        if isexpr(c, :vcat) || isexpr(c, :typed_vcat)
            if isexpr(c.args[1], :parameters)
                @assert length(c.args[1].args) == 1
                condition = popfirst!(c.args).args[1]
            else
                condition = pop!(c.args)
            end
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
        tmp = gensym()
        expr.args[1] = tmp
        code = quote
            let
                $(localvar(i))
                $(localvar(j))
                for $i in $(idxsets[1]), $j in $(idxsets[2])
                    $i <= $j || continue
                    $expr
                    $vname[$i,$j] = $tmp
                    $vname[$j,$i] = $tmp
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
        $varname = $mac
        $code
        nothing
    end
end

localvar(x::Symbol) = _localvar(x)
localvar(x::Expr) = Expr(:block, _localvar(x)...)
_localvar(x::Symbol) = :(local $(esc(x)))
function _localvar(x::Expr)
    @assert x.head in (:escape,:tuple)
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

function addkwargs!(call, kwargs)
    kwsymbol = :(=) # changed by julia PR #19868
    for kw in kwargs
        @assert isexpr(kw, kwsymbol)
        push!(call.args, esc(Expr(:kw, kw.args...)))
    end
end

getname(c::Symbol) = c
getname(c::Nothing) = ()
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

"""
    macro_return(model, code, variable)

Return a block of code that 1. runs the code block `code` in a local scope and 2. returns the value of a local variable named `variable`. `variable` must reference a variable defined by `code`.
"""
function macro_return(code, variable)
    return quote
        let
            # The let block ensures that all variables created behaves like
            # local variables, see https://github.com/JuliaOpt/JuMP.jl/issues/1496
            # To make $variable accessible from outside we need to return it at
            # the end of the block
            $code
            $variable
        end
    end
end

"""
    macro_assign_and_return(code, variable, name;
                            register_fun::Union{Nothing, Function}=nothing,
                            model=nothing)

Return runs `code` in a local scope which returns the value of `variable`
and then assign `variable` to `name`.
If `register_fun` is given, `register_fun(model, name, variable)` is called.
"""
function macro_assign_and_return(code, variable, name;
                                 register_fun::Union{Nothing, Function}=nothing,
                                 model=nothing)
    macro_code = macro_return(code, variable)
    return quote
        $variable = $macro_code
        $(if register_fun !== nothing
              :($register_fun($model, $(quot(name)), $variable))
          end)
        # This assignment should be in the scope calling the macro
        $(esc(name)) = $variable
    end
end

function _construct_constraint!(args...)
    warn("_construct_constraint! is deprecated. Use constructconstraint! instead")
    constructconstraint!(args...)
end

# two-argument constructconstraint! is used for one-sided constraints.
# Right-hand side is zero.
constructconstraint!(v::Variable, sense::Symbol) = constructconstraint!(convert(AffExpr,v), sense)
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
constructconstraint!(x::AbstractArray, sense::Symbol) = constructconstraint!(collect(x), sense)

_vectorize_like(x::Number, y::AbstractArray{AffExpr}) = (ret = similar(y, typeof(x)); fill!(ret, x))
function _vectorize_like(x::AbstractArray{R}, y::AbstractArray{AffExpr}) where R<:Number
    for i in 1:max(ndims(x),ndims(y))
        _size(x,i) == _size(y,i) || error("Unequal sizes for ranged constraint")
    end
    x
end

function constructconstraint!(x::AbstractArray{AffExpr}, lb, ub)
    LB = _vectorize_like(lb,x)
    UB = _vectorize_like(ub,x)
    ret = similar(x, LinearConstraint)
    map!(ret, eachindex(ret)) do i
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

constructconstraint!(x::AbstractMatrix, ::PSDCone) = SDConstraint(x)

constraint_error(args, str) = error("In @constraint($(join(args,","))): ", str)

macro constraint(args...)
    # Pick out keyword arguments
    if isexpr(args[1],:parameters) # these come if using a semicolon
        kwargs = args[1]
        args = args[2:end]
    else
        kwargs = Expr(:parameters)
    end
    kwsymbol = :(=) # changed by julia PR #19868
    append!(kwargs.args, filter(x -> isexpr(x, kwsymbol), collect(args))) # comma separated
    args = filter(x->!isexpr(x, kwsymbol), collect(args))

    if length(args) < 2
        if length(kwargs.args) > 0
            constraint_error(args, "Not enough positional arguments. Do you have a constraint like x=1 instead of x==1?")
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

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat) || length(extra) != 1
    variable = gensym()

    if isa(x, Symbol)
        constraint_error(args, "Incomplete constraint specification $x. Are you missing a comparison (<=, >=, or ==)?")
    end

    (x.head == :block) &&
        constraint_error(args, "Code block passed as constraint. Perhaps you meant to use @constraints instead?")

    # Strategy: build up the code for non-macro addconstraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c, variable)
    # JuMP accepts constraint syntax of the form @constraint(m, foo in bar).
    # This will be rewritten to a call to constructconstraint!(foo, bar). To
    # extend JuMP to accept set-based constraints of this form, it is necessary
    # to add the corresponding methods to constructconstraint!. Note that this
    # will likely mean that bar will be some custom type, rather than e.g. a
    # Symbol, since we will likely want to dispatch on the type of the set
    # appearing in the constraint.
    if isexpr(x, :call)
        if x.args[1] == :in
            @assert length(x.args) == 3
            newaff, parsecode = parseExprToplevel(x.args[2], :q)
            constraintcall = :(addconstraint($m, constructconstraint!($newaff,$(esc(x.args[3])))))
        else
            # Simple comparison - move everything to the LHS
            @assert length(x.args) == 3
            (sense,vectorized) = _canonicalize_sense(x.args[1])
            lhs = :($(x.args[2]) - $(x.args[3]))
            addconstr = (vectorized ? :addVectorizedConstraint : :addconstraint)
            newaff, parsecode = parseExprToplevel(lhs, :q)
            constraintcall = :($addconstr($m, constructconstraint!($newaff,$(quot(sense)))))
        end
        addkwargs!(constraintcall, kwargs.args)
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
        addkwargs!(constraintcall, kwargs.args)
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
                    constraint_error($args, string("Expected ",$lb_str," to be a ", CoefType, "."))
                end
                try
                    ubval = convert(CoefType, $newub)
                catch
                    constraint_error($args, string("Expected ",$ub_str," to be a ", CoefType, "."))
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
    creation_code = getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :ConstraintRef)
    if anonvar
        macro_code = macro_return(creation_code, variable)
    else
        macro_code = macro_assign_and_return(creation_code, variable, getname(c),
                                             register_fun = registercon,
                                             model = m)
    end
end

macro SDconstraint(m, x)
    m = esc(m)

    if isa(x, Symbol)
        error("in @SDconstraint: Incomplete constraint specification $x. Are you missing a comparison (<= or >=)?")
    end

    (x.head == :block) &&
        error("Code block passed as constraint.")
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
    cref = gensym()
    code = quote
        q = zero(AffExpr)
        $parsecode
        $cref = addconstraint($m, constructconstraint!($newaff, PSDCone()))
    end
    assert_validmodel(m, macro_return(code, cref))
end

macro LinearConstraint(x)
    (x.head == :block) &&
        error("Code block passed as constraint. Perhaps you meant to use @LinearConstraints instead?")

    if isexpr(x, :call) && length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[1])
        # Simple comparison - move everything to the LHS
        vectorized &&
            error("in @LinearConstraint ($(string(x))): Cannot add vectorized constraints")
        lhs = :($(x.args[2]) - $(x.args[3]))
        c = gensym()
        code = quote
            newaff = @Expression($(esc(lhs)))
            $c = constructconstraint!(newaff,$(quot(sense)))
            isa($c, LinearConstraint) ||
                error("Constraint in @LinearConstraint is really a $(typeof($c))")
        end
        return macro_return(code, c)
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
        c = gensym()
        code = quote
            if !isa($(esc(lb)),Number)
                error(string("in @LinearConstraint (",$(string(x)),"): expected ",$(string(lb))," to be a number."))
            elseif !isa($(esc(ub)),Number)
                error(string("in @LinearConstraint (",$(string(x)),"): expected ",$(string(ub))," to be a number."))
            end
            newaff = @Expression($(esc(x.args[3])))
            offset = newaff.constant
            newaff.constant = 0.0
            isa(newaff,AffExpr) || error("Ranged quadratic constraints are not allowed")
            $c = LinearConstraint(newaff,$(esc(lb))-offset,$(esc(ub))-offset)
        end
        return macro_return(code, c)
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

    if isexpr(x, :call) && length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[1])
        # Simple comparison - move everything to the LHS
        vectorized &&
            error("in @QuadConstraint ($(string(x))): Cannot add vectorized constraints")
        lhs = :($(x.args[2]) - $(x.args[3]))
        q = gensym()
        code = quote
            newaff = @Expression($(esc(lhs)))
            $q = constructconstraint!(newaff,$(quot(sense)))
            isa($q, QuadConstraint) || error("Constraint in @QuadConstraint is really a $(typeof($q))")
        end
        return macro_return(code, q)
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

    if isexpr(x, :call) && length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[1])
        # Simple comparison - move everything to the LHS
        vectorized &&
            error("in @SOCConstraint ($(string(x))): Cannot add vectorized constraints")
        lhs = :($(x.args[2]) - $(x.args[3]))
        q = gensym()
        code = quote
            newaff = @Expression($(esc(lhs)))
            $q = constructconstraint!(newaff,$(quot(sense)))
            isa($q, SOCConstraint) || error("Constraint in @SOCConstraint is really a $(typeof($q))")
        end
        return macro_return(code, q)
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
    if VERSION ≥ v"0.7-"
        @eval begin
            macro $mac(x)
                x.head == :block || error(string("Invalid syntax for @", $(string(mac))))
                @assert x.args[1] isa LineNumberNode || x.args[1].head == :line
                lastline = x.args[1]
                code = Expr(:vect)
                for it in x.args
                    if it isa LineNumberNode
                        lastline = it
                    elseif it.head == :comparison || (it.head == :call && it.args[1] in (:<=,:≤,:>=,:≥,:(==))) # regular constraint
                        mac = esc(Expr(:macrocall, $(quot(sym)), lastline, it))
                        push!(code.args, mac)
                    elseif it.head == :tuple # constraint ref
                        if all([isexpr(arg,:comparison) for arg in it.args]...)
                            # the user probably had trailing commas at end of lines, e.g.
                            # @LinearConstraints(m, begin
                            #     x <= 1,
                            #     x >= 1
                            # end)
                            error(string("Invalid syntax in @", $(string(mac)), ". Do you have commas at the end of a line specifying a constraint?"))
                        end
                        error("@", string($(string(mac)), " does not currently support the two argument syntax for specifying groups of constraints in one line."))
                    else
                        error("Unexpected constraint expression $it")
                    end
                end
                return code
            end
        end
    else
        @eval begin
            macro $mac(x)
                x.head == :block || error(string("Invalid syntax for @", $(string(mac))))
                @assert x.args[1] isa LineNumberNode || x.args[1].head == :line
                code = Expr(:vect)
                for it in x.args
                    if it isa LineNumberNode || it.head == :line
                        # do nothing
                    elseif it.head == :comparison || (it.head == :call && it.args[1] in (:<=,:≤,:>=,:≥,:(==))) # regular constraint
                        push!(code.args, Expr(:macrocall, $sym, esc(it)))
                    elseif it.head == :tuple # constraint ref
                        if all([isexpr(arg,:comparison) for arg in it.args]...)
                            # the user probably had trailing commas at end of lines, e.g.
                            # @LinearConstraints(m, begin
                            #     x <= 1,
                            #     x >= 1
                            # end)
                            error(string("Invalid syntax in @", $(string(mac)), ". Do you have commas at the end of a line specifying a constraint?"))
                        end
                        error("@", string($(string(mac)), " does not currently support the two argument syntax for specifying groups of constraints in one line."))
                    else
                        error("Unexpected constraint expression $it")
                    end
                end
                return code
            end
        end
    end
end

add_JuMP_prefix(s::Symbol) = Expr(:., JuMP, :($(QuoteNode(s))))

for (mac,sym) in [(:constraints,  Symbol("@constraint")),
                  (:NLconstraints,Symbol("@NLconstraint")),
                  (:SDconstraints,Symbol("@SDconstraint")),
                  (:variables,Symbol("@variable")),
                  (:expressions, Symbol("@expression")),
                  (:NLexpressions, Symbol("@NLexpression"))]
    if VERSION ≥ v"0.7-"
        @eval begin
            macro $mac(m, x)
                x.head == :block || error("Invalid syntax for @",$(string(mac)))
                @assert isa(x.args[1], LineNumberNode)
                lastline = x.args[1]
                code = quote end
                for it in x.args
                    if isa(it, LineNumberNode)
                        lastline = it
                    elseif isexpr(it, :tuple) # line with commas
                        args = []
                        # Keyword arguments have to appear like:
                        # x, (start = 10, lowerbound = 5)
                        # because of the precedence of "=".
                        for ex in it.args
                            if isexpr(ex, :tuple) # embedded tuple
                                append!(args, ex.args)
                            else
                                push!(args, ex)
                            end
                        end
                        mac = esc(Expr(:macrocall, $(add_JuMP_prefix(sym)), lastline, m, args...))
                        push!(code.args, mac)
                    else # stand-alone symbol or expression
                        push!(code.args, esc(Expr(:macrocall, $(add_JuMP_prefix(sym)), lastline, m, it)))
                    end
                end
                push!(code.args, :(nothing))
                return code
            end
        end
    else
        @eval begin
            macro $mac(m, x)
                x.head == :block || error("Invalid syntax for @",$(string(mac)))
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
                            isexpr(ex, :(=))
                            push!(args_esc, esc(ex))
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
    return assert_validmodel(m, macro_return(code, newaff))
end

# Return a standalone, unnamed expression
# ex = @Expression(2x + 3y)
# Currently for internal use only.
macro Expression(x)
    newaff, parsecode = parseExprToplevel(x, :q)
    code = quote
        q = 0.0
        $parsecode
    end
    return macro_return(code, newaff)
end


macro expression(args...)
    if length(args) == 3
        m = esc(args[1])
        c = args[2]
        x = args[3]
    elseif length(args) == 2
        m = esc(args[1])
        c = gensym()
        x = args[2]
    else
        error("@expression: needs at least two arguments.")
    end

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    variable = gensym()

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
    if anonvar
        macro_code = macro_return(code, variable)
    else
        macro_code = macro_assign_and_return(code, variable, getname(c))
    end
    if m === nothing # deprecated usage
        macro_code
    else
        return assert_validmodel(m, macro_code)
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
esc_nonconstant(x::Expr) = isexpr(x,:quote) ? x : esc(x)
esc_nonconstant(x) = esc(x)

# Returns the type of what `constructvariable!` would return with these starting positional arguments.
variabletype(m::Model) = Variable
# Returns a new variable belonging to the model `m`. Additional positional arguments can be used to dispatch the call to a different method.
# The return type should only depends on the positional arguments for `variabletype` to make sense.
function constructvariable!(m::Model, _error::Function, lowerbound::Number, upperbound::Number, category::Symbol, objective::Number, inconstraints::Vector, coefficients::Vector{Float64}, basename::AbstractString, start::Number; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    Variable(m, lowerbound, upperbound, category == :Default ? :Cont : category, objective, inconstraints, coefficients, basename, start)
end

function constructvariable!(m::Model, _error::Function, lowerbound::Number, upperbound::Number, category::Symbol, basename::AbstractString, start::Number; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    Variable(m, lowerbound, upperbound, category == :Default ? :Cont : category, basename, start)
end

const EMPTYSTRING = ""

variable_error(args, str) = error("In @variable($(join(args,","))): ", str)

# @variable(m, expr, extra...; kwargs...)
# where `extra` is a list of extra positional arguments and `kwargs` is a list of keyword arguments.
#
# It creates a new variable (resp. a container of new variables) belonging to the model `m` using `constructvariable!` to create the variable (resp. each variable of the container).
# The following modifications will be made to the arguments before they are passed to `constructvariable!`:
# * The `expr` argument will not be passed but the expression will be parsed to determine the kind of container needed (if one is needed) and
#   additional information that will alter what is passed with the keywords `lowerbound`, `upperbound`, `basename` and `start`.
# * The `SDP` and `Symmetric` positional arguments in `extra` will not be passed to `constructvariable!`. Instead,
#    * the `Symmetric` argument will check that the container is symmetric and only allocate one variable for each pair of non-diagonal entries.
#    * the `SDP` argument will do the same as `Symmetric` but in addition it will specify that the variables created belongs to the SDP cone in the `varCones` field of the model.
#   Moreover, if a Cont, Int, Bin, SemiCont or SemiInt is passed in `extra`, it will alter what is passed with the keyword `category`.
# * The keyword arguments start, objective, inconstraints, coefficients, basename, lowerbound, upperbound, category may not be passed as is to
#   `constructvariable!` since they may be altered by the parsing of `expr` and we may need to pass it pointwise if it is a container since
#   `constructvariable!` is called separately for each variable of the container. Moreover it will be passed as positional argument to `constructvariable!`.
#   If `objective`, `inconstraints` and `coefficients` are not present, they won't be passed as positional argument but for the other five arguments,
#   default values are passed when they are not present.
# * A custom error function is passed as positional argument to print the full @variable call before the error message.
#
# Examples (... is the custom error function):
# * `@variable(m, x >= 0)` is equivalent to `x = constructvariable!(m, msg -> error("In @variable(m, x >= 0): ", msg), 0, Inf, :Cont, "x", NaN)
# * `@variable(m, x[1:N,1:N], Symmetric, Poly(X))` is equivalent to
#   ```
#   x = Matrix{...}(N, N)
#   for i in 1:N
#       for j in 1:N
#           x[i,j] = x[j,i] = constructvariable!(m, Poly(X), msg -> error("In @variable(m, x[1:N,1:N], Symmetric, Poly(X)): ", msg), -Inf, Inf, :Cont, "", NaN)
#       end
#   end
#   ```
macro variable(args...)
    _error(str) = variable_error(args, str)

    m = esc(args[1])

    extra = vcat(args[2:end]...)
    # separate out keyword arguments
    kwsymbol = :(=)
    kwargs = filter(ex->isexpr(ex,kwsymbol), extra)
    extra = filter(ex->!isexpr(ex,kwsymbol), extra)

    # if there is only a single non-keyword argument, this is an anonymous
    # variable spec and the one non-kwarg is the model
    if length(kwargs) == length(args)-1
        x = gensym()
        anon_singleton = true
    else
        x = popfirst!(extra)
        if x in [:Cont,:Int,:Bin,:SemiCont,:SemiInt,:SDP]
            _error("Ambiguous variable name $x detected. Use the \"category\" keyword argument to specify a category for an anonymous variable.")
        end
        anon_singleton = false
    end

    t = quot(:Default)
    gottype = false
    haslb = false
    hasub = false
    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    explicit_comparison = false
    if isexpr(x,:comparison) # two-sided
        explicit_comparison = true
        haslb = true
        hasub = true
        if x.args[2] == :>= || x.args[2] == :≥
            # ub >= x >= lb
            x.args[4] == :>= || x.args[4] == :≥ || _error("Invalid variable bounds")
            var = x.args[3]
            lb = esc_nonconstant(x.args[5])
            ub = esc_nonconstant(x.args[1])
        elseif x.args[2] == :<= || x.args[2] == :≤
            # lb <= x <= u
            var = x.args[3]
            (x.args[4] != :<= && x.args[4] != :≤) &&
                _error("Expected <= operator after variable name.")
            lb = esc_nonconstant(x.args[1])
            ub = esc_nonconstant(x.args[5])
        else
            _error("Use the form lb <= ... <= ub.")
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
            t = quot(:Fixed)
        else
            # Its a comparsion, but not using <= ... <=
            _error("Unexpected syntax $(string(x)).")
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

    if !isa(getname(var),Symbol) && !anonvar
        warn_once("Expression $(getname(var)) should not be used as a variable name. Use the \"anonymous\" syntax $(getname(var)) = @variable(m, ...) instead.")
    end

    # process keyword arguments
    value = NaN
    obj = nothing
    inconstraints = nothing
    coefficients = nothing
    extra_kwargs = []
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
            haslb && _error("Cannot specify variable lowerbound twice")
            lb = esc_nonconstant(ex.args[2])
            haslb = true
        elseif kwarg == :upperbound
            hasub && _error("Cannot specify variable upperbound twice")
            ub = esc_nonconstant(ex.args[2])
            hasub = true
        elseif kwarg == :category
            (t == quot(:Fixed)) && _error("Unexpected extra arguments when declaring a fixed variable")
            t = esc_nonconstant(ex.args[2])
            gottype = true
        else
            push!(extra_kwargs, ex)
        end
    end

    if (obj !== nothing || inconstraints !== nothing || coefficients !== nothing) &&
       (obj === nothing || inconstraints === nothing || coefficients === nothing)
        _error("Must provide 'objective', 'inconstraints', and 'coefficients' arguments all together for column-wise modeling")
    end

    sdp = any(t -> (t == :SDP), extra)
    symmetric = (sdp || any(t -> (t == :Symmetric), extra))
    extra = filter(x -> (x != :SDP && x != :Symmetric), extra) # filter out SDP and sym tag
    for ex in extra
        if ex in var_cats
            if t != quot(:Default) && t != quot(ex)
                _error("A variable cannot be both of category $cat and $category. Please specify only one category.")
            else
                t = quot(ex)
            end
        end
    end
    extra = esc.(filter(ex -> !(ex in var_cats), extra))

    # Handle the column generation functionality
    if coefficients !== nothing
        !isa(var,Symbol) &&
        _error("Can only create one variable at a time when adding to existing constraints.")

        variablecall = :( constructvariable!($m, $(extra...), $_error, $lb, $ub, $t, $obj, $inconstraints, $coefficients, string($quotvarname), $value) )
        addkwargs!(variablecall, extra_kwargs)
        code = :($variable = $variablecall)
        if anonvar
            macro_code = macro_return(code, variable)
        else
            macro_code = macro_assign_and_return(code, variable, getname(var))
        end
        return assert_validmodel(m, macro_code)
    end

    if isa(var,Symbol)
        # Easy case - a single variable
        sdp && _error("Cannot add a semidefinite scalar variable")
        variablecall = :( constructvariable!($m, $(extra...), $_error, $lb, $ub, $t, string($quotvarname), $value) )
        addkwargs!(variablecall, extra_kwargs)
        code = :($variable = $variablecall)
        code = :($variable = $variablecall)
        if anonvar
            macro_code = macro_return(code, variable)
        else
            macro_code = macro_assign_and_return(code, variable, getname(var),
                                                 register_fun = registervar,
                                                 model = m)
        end
        return assert_validmodel(m, macro_code)
    end
    isa(var,Expr) || _error("Expected $var to be a variable name")

    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(var, variable)
    clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? () : idxsets[i])

    # Code to be used to create each variable of the container.
    variablecall = :( constructvariable!($m, $(extra...), $_error, $lb, $ub, $t, EMPTYSTRING, $value) )
    addkwargs!(variablecall, extra_kwargs)
    code = :( $(refcall) = $variablecall )
    # Determine the return type of constructvariable!. This is needed to create the container holding them.
    vartype = :( variabletype($m, $(extra...)) )

    if symmetric
        # Sanity checks on SDP input stuff
        condition == :() ||
            _error("Cannot have conditional indexing for SDP variables")
        length(idxvars) == length(idxsets) == 2 ||
            _error("SDP variables must be 2-dimensional")
        !symmetric || (length(idxvars) == length(idxsets) == 2) ||
            _error("Symmetric variables must be 2-dimensional")
        hasdependentsets(idxvars, idxsets) &&
            _error("Cannot have index dependencies in symmetric variables")
        for _rng in idxsets
            isexpr(_rng, :escape) ||
                _error("Internal error 1")
            rng = _rng.args[1] # undo escaping
            if VERSION >= v"0.7-"
                (isexpr(rng,:call) && length(rng.args) == 3 && rng.args[1] == :(:) && rng.args[2] == 1) ||
                    _error("Index sets for SDP variables must be ranges of the form 1:N")
            else
                (isexpr(rng,:(:)) && rng.args[1] == 1 && length(rng.args) == 2) ||
                    _error("Index sets for SDP variables must be ranges of the form 1:N")
            end
        end

        if !(lb == -Inf && ub == Inf)
            _error("Semidefinite or symmetric variables cannot be provided bounds")
        end
        @static if VERSION >= v"0.7-"
            # 1:3 is parsed as (:call, :, 1, 3)
            dimension_check = :($(esc(idxsets[1].args[1].args[3])) ==
                                $(esc(idxsets[2].args[1].args[3])))
        else
            # 1:3 is parsed as (:, 1, 3)
            dimension_check = :($(esc(idxsets[1].args[1].args[2])) ==
                                $(esc(idxsets[2].args[1].args[2])))
        end
        creation_code = quote
            $dimension_check || error("Cannot construct symmetric variables with nonsquare dimensions")
            (issymmetric($lb) && issymmetric($ub)) || error("Bounds on symmetric  variables must be symmetric")
            $(getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, vartype; lowertri=symmetric))
            $(if sdp
                quote
                    push!($(m).varCones, (:SDP, first($variable).col : last($variable).col))
                end
            end)
            push!($(m).dictList, $variable)
            !$anonvar && registervar($m, $quotvarname, $variable)
            storecontainerdata($m, $variable, $quotvarname,
                               $(Expr(:tuple,idxsets...)),
                               $idxpairs, $(quot(condition)))
        end
    else
        coloncheckcode = Expr(:call,:coloncheck,refcall.args[2:end]...)
        code = :($coloncheckcode; $code)
        creation_code = quote
            $(getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, vartype))
            isa($variable, JuMPContainer) && pushmeta!($variable, :model, $m)
            push!($(m).dictList, $variable)
            !$anonvar && registervar($m, $quotvarname, $variable)
            storecontainerdata($m, $variable, $quotvarname,
                               $(Expr(:tuple,map(clear_dependencies,1:length(idxsets))...)),
                               $idxpairs, $(quot(condition)))
        end
    end
    if anonvar
        macro_code = macro_return(creation_code, variable)
    else
        macro_code = macro_assign_and_return(creation_code, variable, getname(var))
    end
    return assert_validmodel(m, macro_code)
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
            $(gendict(esc(varname), :ConstraintRef, esc.(idxsets)...))
            nothing
        end
        return code
    end
end

macro NLobjective(m, sense, x)
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    ex = gensym()
    code = quote
        initNLP($(esc(m)))
        setobjectivesense($(esc(m)), $(esc(sense)))
        $ex = $(processNLExpr(m, x))
        $(esc(m)).nlpdata.nlobj = $ex
        $(esc(m)).obj = zero(QuadExpr)
        $(esc(m)).internalModelLoaded = false
    end
    return assert_validmodel(esc(m), macro_return(code, ex))
end

macro NLconstraint(m, x, extra...)
    esc_m = esc(m)
    # Two formats:
    # - @NLconstraint(m, a*x <= 5)
    # - @NLconstraint(m, myref[a=1:5], sin(x^a) <= 5)
    length(extra) > 1 && error("in @NLconstraint: too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : gensym()
    x = length(extra) == 1 ? extra[1] : x

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat) || length(extra) != 1
    variable = gensym()

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
            c = NonlinearConstraint($(processNLExpr(m, lhs)), $lb, $ub)
            push!($esc_m.nlpdata.nlconstr, c)
            $(refcall) = ConstraintRef{Model,NonlinearConstraint}($esc_m, length($esc_m.nlpdata.nlconstr))
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
            c = NonlinearConstraint($(processNLExpr(m, x.args[3])), $(esc(lb)), $(esc(ub)))
            push!($esc_m.nlpdata.nlconstr, c)
            $(refcall) = ConstraintRef{Model,NonlinearConstraint}($esc_m, length($esc_m.nlpdata.nlconstr))
        end
    else
        # Unknown
        error("in @NLconstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2")
    end
    looped = getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :(ConstraintRef{Model,NonlinearConstraint}))
    creation_code = quote
        initNLP($esc_m)
        $esc_m.internalModelLoaded = false
        $looped
    end
    if anonvar
        macro_code = macro_return(creation_code, variable)
    else
        macro_code = macro_assign_and_return(creation_code, variable, getname(c),
                                             register_fun = registercon,
                                             model = esc_m)
    end
    return assert_validmodel(esc_m, macro_code)
end

macro NLexpression(args...)
    if length(args) <= 1
        error("in @NLexpression: To few arguments ($(length(args))); must pass the model and nonlinear expression as arguments.")
    elseif length(args) == 2
        m, x = args
        c = gensym()
    elseif length(args) == 3
        m, c, x = args
    else
        error("in @NLexpression: To many arguments ($(length(args))).")
    end

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    variable = gensym()

    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c, variable)
    code = quote
        $(refcall) = NonlinearExpression($(esc(m)), $(processNLExpr(m, x)))
    end
    creation_code = getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :NonlinearExpression)
    if anonvar
        macro_code = macro_return(creation_code, variable)
    else
        macro_code = macro_assign_and_return(creation_code, variable, getname(c))
    end
    return assert_validmodel(esc(m), macro_code)
end

# syntax is @NLparameter(m, p[i=1] == 2i)
macro NLparameter(m, ex)
    m = esc(m)
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

    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c, variable)
    code = quote
        $(refcall) = newparameter($m, $(esc(x)))
    end
    creation_code = getloopedcode(variable, code, condition, idxvars, idxsets, idxpairs, :NonlinearParameter)
    if anonvar
        macro_code = macro_return(creation_code, variable)
    else
        macro_code = macro_assign_and_return(creation_code, variable, getname(c))
    end
    return assert_validmodel(m, macro_code)
end
