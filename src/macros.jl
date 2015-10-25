#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta

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
function buildrefsets(expr::Expr)
    c = copy(expr)
    isexpr(c,:ref) || isexpr(c,:typed_vcat) || error("Unrecognized name in construction macro; expected $(string(c)) to be of the form name[...]")
    idxvars = Any[]
    idxsets = Any[]
    idxpairs = IndexPair[]
    # Creating an indexed set of refs
    cname = shift!(c.args)
    refcall = Expr(:ref,esc(cname))
    condition = :()
    if isexpr(c, :typed_vcat)
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

buildrefsets(c::Symbol)  = (esc(c), Any[], Any[], IndexPair[], :())
buildrefsets(c::Void) = (gensym(), Any[], Any[], IndexPair[], :())

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
function getloopedcode(c::Expr, code, condition, idxvars, idxsets, idxpairs, sym; lowertri=false)
    varname = getname(c)
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
        mac = :($(esc(varname)) = JuMPDict{$(sym),$N}())
    else
        mac = Expr(:macrocall,symbol("@gendict"),esc(varname),sym,idxpairs,idxsets...)
    end
    return quote
        $mac
        $code
        nothing
    end
end

getloopedcode(c, code, condition, idxvars, idxsets, idxpairs, sym) = code

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
getname(c::Expr) = c.args[1]

validmodel(m::Model, name) = nothing
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

constructconstraint!(q::QuadExpr, lb, ub) = error("Two-sided quadratic constraints not supported. (Try @addNLConstraint instead.)")

macro addConstraint(args...)
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
            error("in @addConstraint($(join(args,','))) with keyword arguments: ($(join(kwargs.args,','))): not enough arguments")
        else
            error("in @addConstraint($(join(args,','))): not enough arguments")
        end
    end
    m = args[1]
    x = args[2]
    extra = args[3:end]

    m = esc(m)
    # Two formats:
    # - @addConstraint(m, a*x <= 5)
    # - @addConstraint(m, myref[a=1:5], a*x <= 5)
    length(extra) > 1 && error("in @addConstraint: too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : nothing
    x = length(extra) == 1 ? extra[1] : x

    (x.head == :block) &&
        error("Code block passed as constraint. Perhaps you meant to use addConstraints instead?")
    (x.head != :comparison) &&
        error("in @addConstraint ($(string(x))): expected comparison operator (<=, >=, or ==).")

    # Strategy: build up the code for non-macro addconstraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs
    crefflag = isa(c,Expr)
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c)
    # Build the constraint
    if length(x.args) == 3
        # Simple comparison - move everything to the LHS
        (sense,vectorized) = _canonicalize_sense(x.args[2])
        lhs = :($(x.args[1]) - $(x.args[3]))
        addconstr = (vectorized ? :addVectorizedConstraint : :addConstraint)
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
    elseif length(x.args) == 5
        # Ranged row
        (lsign,lvectorized) = _canonicalize_sense(x.args[2])
        (rsign,rvectorized) = _canonicalize_sense(x.args[4])
        if (lsign != :(<=)) || (rsign != :(<=))
            error("in @addConstraint ($(string(x))): only ranged rows of the form lb <= expr <= ub are supported.")
        end
        ((vectorized = lvectorized) == rvectorized) || error("in @addConstraint ($(string(x))): signs are inconsistently vectorized")
        addconstr = (lvectorized ? :addVectorizedConstraint : :addConstraint)
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
                    error(string("in @addConstraint (",$x_str,"): expected ",$lb_str," to be a ", CoefType, "."))
                end
                try
                    ubval = convert(CoefType, $newub)
                catch
                    error(string("in @addConstraint (",$x_str,"): expected ",$ub_str," to be a ", CoefType, "."))
                end
            end
        end
        code = quote
            $code
            $(refcall) = $constraintcall
        end
    else
        # Unknown
        error("in @addConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2\n" * "       lb <= expr <= ub")
    end
    return assert_validmodel(m, getloopedcode(c, code, condition, idxvars, idxsets, idxpairs, :ConstraintRef))
end

macro addSDPConstraint(m, x)
    m = esc(m)

    (x.head == :block) &&
        error("Code block passed as constraint.")
    (x.head != :comparison) &&
        error("in @addSDPConstraint ($(string(x))): expected comparison operator (<=, or >=).")

    length(x.args) == 3 || error("in @addSDPConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n")
    # Build the constraint
    # Simple comparison - move everything to the LHS
    sense = x.args[2]
    if sense == :⪰
        sense = :(>=)
    elseif sense == :⪯
        sense = :(<=)
    end
    sense,_ = _canonicalize_sense(sense)
    lhs = :()
    if sense == :(>=)
        lhs = :($(x.args[1]) - $(x.args[3]))
    elseif sense == :(<=)
        lhs = :($(x.args[3]) - $(x.args[1]))
    else
        error("Invalid sense $sense in SDP constraint")
    end
    newaff, parsecode = parseExprToplevel(lhs, :q)
    assert_validmodel(m, quote
        q = zero(AffExpr)
        $parsecode
        c = SDPConstraint($newaff)
        push!($(m).sdpconstr, c)
        c
    end)
end

macro LinearConstraint(x)
    (x.head == :block) &&
        error("Code block passed as constraint. Perhaps you meant to use @LinearConstraints instead?")
    (x.head != :comparison) &&
        error("in @LinearConstraint ($(string(x))): expected comparison operator (<=, >=, or ==).")

    if length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[2])
        # Simple comparison - move everything to the LHS
        vectorized &&
            error("in @LinearConstraint ($(string(x))): Cannot add vectorized constraints")
        lhs = :($(x.args[1]) - $(x.args[3]))
        return quote
            newaff = @defExpr($(esc(lhs)))
            c = constructconstraint!(newaff,$(quot(sense)))
            isa(c, LinearConstraint) ||
                error("Constraint in @LinearConstraint is really a $(typeof(c))")
            c
        end
    elseif length(x.args) == 5
        # Ranged row
        (lsense,lvectorized) = _canonicalize_sense(x.args[2])
        (rsense,rvectorized) = _canonicalize_sense(x.args[4])
        if (lsense != :<=) || (rsense != :<=)
            error("in @addConstraint ($(string(x))): only ranged rows of the form lb <= expr <= ub are supported.")
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
            newaff = @defExpr($(esc(x.args[3])))
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
    (x.head != :comparison) &&
        error("in @QuadConstraint ($(string(x))): expected comparison operator (<=, >=, or ==).")

    if length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[2])
        # Simple comparison - move everything to the LHS
        vectorized &&
            error("in @QuadConstraint ($(string(x))): Cannot add vectorized constraints")
        lhs = :($(x.args[1]) - $(x.args[3]))
        return quote
            newaff = @defExpr($(esc(lhs)))
            q = constructconstraint!(newaff,$(quot(sense)))
            isa(q, QuadConstraint) || error("Constraint in @QuadConstraint is really a $(typeof(q))")
            q
        end
    elseif length(x.args) == 5
        error("Ranged quadratic constraints are not allowed")
    else
        # Unknown
        error("in @QuadConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2\n")
    end
end

macro SOCConstraint(x)
    (x.head == :block) &&
        error("Code block passed as constraint. Perhaps you meant to use @SOCConstraints instead?")
    (x.head != :comparison) &&
        error("in @SOCConstraint ($(string(x))): expected comparison operator (<=, >=, or ==).")

    if length(x.args) == 3
        (sense,vectorized) = _canonicalize_sense(x.args[2])
        # Simple comparison - move everything to the LHS
        vectorized &&
            error("in @SOCConstraint ($(string(x))): Cannot add vectorized constraints")
        lhs = :($(x.args[1]) - $(x.args[3]))
        return quote
            newaff = @defExpr($(esc(lhs)))
            q = constructconstraint!(newaff,$(quot(sense)))
            isa(q, SOCConstraint) || error("Constraint in @SOCConstraint is really a $(typeof(q))")
            q
        end
    elseif length(x.args) == 5
        error("Ranged second-order cone constraints are not allowed")
    else
        # Unknown
        error("in @SOCConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n")
    end
end

for (mac,sym) in [(:LinearConstraints, symbol("@LinearConstraint")),
                  (:QuadConstraints,   symbol("@QuadConstraint")),
                  (:SOCConstraints,    symbol("@SOCConstraint"))]
    @eval begin
        macro $mac(x)
            x.head == :block || error("Invalid syntax for @$mac")
            @assert x.args[1].head == :line
            code = Expr(:vect)
            for it in x.args
                if it.head == :line
                    # do nothing
                elseif it.head == :comparison # regular constraint
                    push!(code.args, Expr(:macrocall, $sym, esc(it)))
                elseif it.head == :tuple # constraint ref
                    error("@$mac does not currently support groups of constraints")
                else
                    error("Unexpected constraint expression $it")
                end
            end
            return code
        end
    end
end

for (mac,sym) in [(:addConstraints,  symbol("@addConstraint")),
                  (:addNLConstraints,symbol("@addNLConstraint"))]
    @eval begin
        macro $mac(m, x)
            x.head == :block || error("Invalid syntax for @$mac")
            @assert x.args[1].head == :line
            code = quote end
            for it in x.args
                if it.head == :line
                    # do nothing
                elseif it.head == :comparison # regular constraint
                    mac = Expr(:macrocall,$(quot(sym)), esc(m), esc(it))
                    code = quote
                            $code
                            $mac
                            end
                elseif it.head == :tuple # constraint ref or kwargs
                    for i in 1:length(it.args)
                        if isexpr(it.args[i], :(=))
                            it.args[i] = Expr(:kw, it.args[i].args[1], esc(it.args[i].args[2]))
                        else
                            it.args[i] = esc(it.args[i])
                        end
                    end
                    mac = Expr(:macrocall,$(quot(sym)), esc(m), it.args...)
                    code = quote
                            $code
                            $mac
                            end
                else
                    error("Unexpected constraint expression $it")
                end
            end
            return quote
                    $code
                    nothing
                end
        end
    end
end

macro setObjective(m, args...)
    m = esc(m)
    if length(args) != 2
        # Either just an objective sene, or just an expression.
        error("in @setObjective: needs three arguments: model, objective sense (Max or Min) and expression.")
    end
    sense, x = args
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    newaff, parsecode = parseExprToplevel(x, :q)
    code = quote
        q = zero(AffExpr)
        $parsecode
        setObjective($m, $(esc(sense)), $newaff)
    end
    return assert_validmodel(m, code)
end

macro defExpr(args...)
    if length(args) == 1
        c = nothing
        x = args[1]
    elseif length(args) == 2
        c = args[1]
        x = args[2]
    else
        error("in @defExpr: needs either one or two arguments.")
    end

    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c)
    newaff, parsecode = parseExprToplevel(x, :q)
    code = quote
        q = 0.0
        $parsecode
    end
    if isa(c,Expr)
        code = quote
            $code
            (isa($newaff,AffExpr) || isa($newaff,Number) || isa($newaff,Variable)) || error("Collection of expressions with @defExpr must be linear. For quadratic expressions, use your own array.")
        end
    end
    code = quote
        $code
        $(refcall) = $newaff
    end
    return getloopedcode(c, code, condition, idxvars, idxsets, idxpairs, :AffExpr)
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

const EMPTYSTRING = utf8("")

macro defVar(args...)
    length(args) <= 1 &&
        error("in @defVar: expected model as first argument, then variable information.")
    m = esc(args[1])
    x = args[2]
    extra = vcat(args[3:end]...)

    t = :Cont
    gottype = false
    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    if isexpr(x,:comparison)
        # We have some bounds
        if x.args[2] == :>= || x.args[2] == :≥
            if length(x.args) == 5
                # ub >= x >= lb
                x.args[4] == :>= || x.args[4] == :≥ || error("Invalid variable bounds")
                var = x.args[3]
                lb = esc_nonconstant(x.args[5])
                ub = esc_nonconstant(x.args[1])
            else
                # x >= lb
                var = x.args[1]
                @assert length(x.args) == 3
                lb = esc_nonconstant(x.args[3])
                ub = Inf
            end
        elseif x.args[2] == :<= || x.args[2] == :≤
            if length(x.args) == 5
                # lb <= x <= u
                var = x.args[3]
                (x.args[4] != :<= && x.args[4] != :≤) &&
                    error("in @defVar ($var): expected <= operator after variable name.")
                lb = esc_nonconstant(x.args[1])
                ub = esc_nonconstant(x.args[5])
            else
                # x <= ub
                var = x.args[1]
                # NB: May also be lb <= x, which we do not support
                #     We handle this later in the macro
                @assert length(x.args) == 3
                ub = esc_nonconstant(x.args[3])
                lb = -Inf
            end
        elseif x.args[2] == :(==)
            # fixed variable
            var = x.args[1]
            @assert length(x.args) == 3
            lb = esc(x.args[3])
            ub = esc(x.args[3])
            gottype = true
            t = :Fixed
        else
            # Its a comparsion, but not using <= ... <=
            error("in @defVar ($(string(x))): use the form lb <= ... <= ub.")
        end
    else
        # No bounds provided - free variable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        var = x
        lb = -Inf
        ub = Inf
    end

    # separate out keyword arguments
    kwargs = filter(ex->isexpr(ex,:kw), extra)
    extra = filter(ex->!isexpr(ex,:kw), extra)

    # process keyword arguments
    value = NaN
    obj = nothing
    inconstraints = nothing
    coefficients = nothing
    for ex in kwargs
        if ex.args[1] == :start
            value = esc(ex.args[2])
        elseif ex.args[1] == :objective
            obj = esc(ex.args[2])
        elseif ex.args[1] == :inconstraints
            inconstraints = esc(ex.args[2])
        elseif ex.args[1] == :coefficients
            coefficients = esc(ex.args[2])
        else
            error("in @defVar ($var): Unrecognized keyword argument $(ex.args[1])")
        end
    end

    if (obj !== nothing || inconstraints !== nothing || coefficients !== nothing) &&
       (obj === nothing || inconstraints === nothing || coefficients === nothing)
        error("in @defVar ($var): Must provide 'objective', 'inconstraints', and 'coefficients' arguments all together for column-wise modeling")
    end

    sdp = any(t -> (t == :SDP), extra)
    symmetric = (sdp || any(t -> (t == :Symmetric), extra))
    extra = filter(x -> (x != :SDP && x != :Symmetric), extra) # filter out SDP and sym tag

    # Determine variable type (if present).
    # Types: default is continuous (reals)
    if length(extra) > 0
        if t == :Fixed
            error("in @defVar ($var): unexpected extra arguments when declaring a fixed variable")
        end
        if extra[1] in [:Bin, :Int, :SemiCont, :SemiInt]
            gottype = true
            t = extra[1]
        end

        if t == :Bin
            if (lb != -Inf || ub != Inf) && !(lb == 0.0 && ub == 1.0)
            error("in @defVar ($var): bounds other than [0, 1] may not be specified for binary variables.\nThese are always taken to have a lower bound of 0 and upper bound of 1.")
            else
                lb = 0.0
                ub = 1.0
            end
        end

        !gottype && error("in @defVar ($var): syntax error")
    end

    # Handle the column generation functionality
    if coefficients !== nothing
        !isa(var,Symbol) &&
        error("in @defVar ($var): can only create one variable at a time when adding to existing constraints.")

        return assert_validmodel(m, quote
            $(esc(var)) = Variable($m,$lb,$ub,$(quot(t)),$obj,$inconstraints,$coefficients,$(utf8(string(var))),$value)
            nothing
        end)
    end

    if isa(var,Symbol)
        # Easy case - a single variable
        sdp && error("Cannot add a semidefinite scalar variable")
        return assert_validmodel(m, quote
            $(esc(var)) = Variable($m,$lb,$ub,$(quot(t)),$(utf8(string(var))),$value)
            registervar($m, $(quot(var)), $(esc(var)))
        end)
    end
    isa(var,Expr) || error("in @defVar: expected $var to be a variable name")

    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(var)
    clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? nothing : idxsets[i])

    code = :( $(refcall) = Variable($m, $lb, $ub, $(quot(t)), EMPTYSTRING, $value) )
    if symmetric
        # Sanity checks on SDP input stuff
        condition == :() ||
            error("Cannot have conditional indexing for SDP variables")
        length(idxvars) == length(idxsets) == 2 ||
            error("SDP variables must be 2-dimensional")
        !symmetric || (length(idxvars) == length(idxsets) == 2) ||
            error("Symmetric variables must be 2-dimensional")
        hasdependentsets(idxvars, idxsets) &&
            error("Cannot have index dependencies in symmetric variables")
        for _rng in idxsets
            isexpr(_rng, :escape) ||
                error("Internal error 1")
            rng = _rng.args[1] # undo escaping
            (isexpr(rng,:(:)) && rng.args[1] == 1 && length(rng.args) == 2) ||
                error("Index sets for SDP variables must be ranges of the form 1:N")
        end

        if sdp && !(lb == -Inf && ub == Inf)
            error("Semidefinite variables cannot be provided bounds")
        end

        looped = getloopedcode(var, code, condition, idxvars, idxsets, idxpairs, :Variable; lowertri=symmetric)
        varname = esc(getname(var))
        code = quote
            $(esc(idxsets[1].args[1].args[2])) == $(esc(idxsets[2].args[1].args[2])) || error("Cannot construct symmetric variables with nonsquare dimensions")
            (issym($lb) && issym($ub)) || error("Bounds on symmetric variables must be symmetric")
            $looped
            push!($(m).dictList, $varname)
        end
        if sdp
            code = :($code; push!($(m).varCones, (:SDP, first($varname).col : last($varname).col)))
        end
        return assert_validmodel(m, quote
            $code
            registervar($m, $(quot(getname(var))), $varname)
            storecontainerdata($m, $varname, $(quot(getname(var))),
                               $(Expr(:tuple,idxsets...)),
                               $idxpairs, $(quot(condition)))
            $varname
        end)
    else
        looped = getloopedcode(var, code, condition, idxvars, idxsets, idxpairs, :Variable)
        varname = esc(getname(var))
        return assert_validmodel(m, quote
            $looped
            push!($(m).dictList, $varname)
            registervar($m, $(quot(getname(var))), $varname)
            storecontainerdata($m, $varname, $(quot(getname(var))),
                               $(Expr(:tuple,map(clear_dependencies,1:length(idxsets))...)),
                               $idxpairs, $(quot(condition)))
            isa($varname, JuMPContainer) && pushmeta!($varname, :model, $m)
            $varname
        end)
    end
end

storecontainerdata(m::Model, variable, varname, idxsets, idxpairs, condition) =
    m.varData[variable] = JuMPContainerData(varname, idxsets, idxpairs, condition)

macro defConstrRef(var)
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

        mac = Expr(:macrocall,symbol("@gendict"),varname,:ConstraintRef,idxpairs, idxsets...)
        code = quote
            $(esc(mac))
            nothing
        end
        return code
    end
end

macro setNLObjective(m, sense, x)
    m = esc(m)
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    code = quote
        initNLP($m)
        setObjectiveSense($m, $(esc(sense)))
        ex = @processNLExpr($(esc(x)))
        $m.nlpdata.nlobj = ex
        $m.obj = zero(QuadExpr)
        $m.internalModelLoaded = false
        nothing
    end
    return assert_validmodel(m, code)
end

macro addNLConstraint(m, x, extra...)
    m = esc(m)
    # Two formats:
    # - @addNLConstraint(m, a*x <= 5)
    # - @addNLConstraint(m, myref[a=1:5], sin(x^a) <= 5)
    length(extra) > 1 && error("in @addNLConstraint: too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : nothing
    x = length(extra) == 1 ? extra[1] : x

    (x.head != :comparison) &&
        error("in @addNLConstraint ($(string(x))): expected comparison operator (<=, >=, or ==).")

    # Strategy: build up the code for non-macro addconstraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(c)
    # Build the constraint
    if length(x.args) == 3
        # Simple comparison - move everything to the LHS
        op = x.args[2]
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
            error("in @addNLConstraint ($(string(x))): expected comparison operator (<=, >=, or ==).")
        end
        lhs = :($(x.args[1]) - $(x.args[3]))
        code = quote
            c = NonlinearConstraint(@processNLExpr($(esc(lhs))), $lb, $ub)
            push!($m.nlpdata.nlconstr, c)
            push!($m.nlpdata.nlconstrlist, c.terms)
            $(refcall) = ConstraintRef{NonlinearConstraint}($m, length($m.nlpdata.nlconstr))
        end
    elseif length(x.args) == 5
        # ranged row
        if (x.args[2] != :<= && x.args[2] != :≤) || (x.args[4] != :<= && x.args[4] != :≤)
            error("in @addNLConstraint ($(string(x))): only ranged rows of the form lb <= expr <= ub are supported.")
        end
        lb = x.args[1]
        ub = x.args[5]
        code = quote
            if !isa($(esc(lb)),Number)
                error(string("in @addNLConstraint (",$(string(x)),"): expected ",$(string(lb))," to be a number."))
            elseif !isa($(esc(ub)),Number)
                error(string("in @addNLConstraint (",$(string(x)),"): expected ",$(string(ub))," to be a number."))
            end
            c = NonlinearConstraint(@processNLExpr($(esc(x.args[3]))), $(esc(lb)), $(esc(ub)))
            push!($m.nlpdata.nlconstr, c)
            push!($m.nlpdata.nlconstrlist, c.terms)
            $(refcall) = ConstraintRef{NonlinearConstraint}($m, length($m.nlpdata.nlconstr))
        end
    else
        # Unknown
        error("in @addNLConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2\n")
    end
    looped = getloopedcode(c, code, condition, idxvars, idxsets, idxpairs, :(ConstraintRef{NonlinearConstraint}))
    code = quote
        initNLP($m)
        $looped
        $m.internalModelLoaded = false
        nothing
    end

    return assert_validmodel(m, code)
end

macro defNLExpr(x, extra...)
    # Two formats:
    # - @defNLExpr(a*x <= 5)
    # - @defNLExpr(myref[a=1:5], sin(x^a))
    length(extra) > 1 && error("in @defNLExpr: too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : nothing
    x = length(extra) == 1 ? extra[1] : x

    refcall, idxvars, idxsets, idxpairs = buildrefsets(c)
    varname = isexpr(refcall,:ref) ? refcall.args[1] : refcall
    macrocall = Expr(:macrocall, symbol("@parametricExpr"), [esc(v) for v in idxvars]..., esc(x))
    return :($(varname) = $macrocall)
end
