#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta

if VERSION > v"0.4.0-"
    include(joinpath("v0.4","parseExpr_staged.jl"))
else
    include(joinpath("v0.3","parseExpr_0.3.jl"))
end

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
        if isa(s,Expr) && (s.head == :(=) || s.head == :in)
            idxvar = s.args[1]
            idxset = esc(s.args[2])
            push!(idxpairs, IndexPair(s.args[1],s.args[2]))
        else
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
buildrefsets(c::Nothing) = (gensym(), Any[], Any[], IndexPair[], :())

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
function getloopedcode(c::Expr, code, condition, idxvars, idxsets, idxpairs, sym)
    varname = getname(c)
    hascond = (condition != :())

    if hascond
        code = quote
            $(esc(condition)) || continue
            $code
        end
    end

    for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
        code = quote
            for $(esc(idxvar)) in $idxset
                $code
            end
        end
    end
    if hascond || hasdependentsets(idxvars,idxsets)
        # force a JuMPDict
        N = length(idxsets)
        clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? nothing : idxsets[i])
        mac = :($(esc(varname)) = JuMPDict{$(sym),$N}(Dict{NTuple{$N},$sym}(),
                                                        $(quot(varname)),
                                                        $(Expr(:tuple,map(clear_dependencies,1:N)...)),
                                                        $idxpairs,
                                                        $(quot(condition))))
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

getname(c::Symbol) = c
getname(c::Nothing) = ()
getname(c::Expr) = c.args[1]

function assert_validmodel(m, macrocode)
    # assumes m is already escaped
    quote
        isa($m, Model) || error("Expected ", $(quot(m.args[1])), " to be a JuMP model but it has type ", typeof($m))
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

# two-argument _construct_constraint! is used for one-sided constraints.
# Right-hand side is zero.
_construct_constraint!(v::Variable, sense::Symbol) = _construct_constraint(convert(AffExpr,v), sense)
function _construct_constraint!(aff::AffExpr, sense::Symbol)
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

function _construct_constraint!(aff::AffExpr, lb, ub)
    offset = aff.constant
    aff.constant = 0.0
    LinearConstraint(aff, lb-offset, ub-offset)
end

_construct_constraint!(quad::QuadExpr, sense::Symbol) = QuadConstraint(quad, sense)

_construct_constraint!(x::Array, sense::Symbol) = map(c->_construct_constraint!(c,sense), x)

_vectorize_like(x::Number, y::Array{AffExpr}) = fill(x, size(y))
function _vectorize_like{R<:Number}(x::Array{R}, y::Array{AffExpr})
    for i in 1:max(ndims(x),ndims(y))
        size(x,i) == size(y,i) || error("Unequal sizes for ranged constraint")
    end
    x
end

function _construct_constraint!(x::Array{AffExpr}, lb, ub)
    LB = _vectorize_like(lb,x)
    UB = _vectorize_like(ub,x)
    ret = similar(x, LinearConstraint)
    map!(ret, 1:length(ret)) do i
        _construct_constraint!(x[i], LB[i], UB[i])
    end
end

# three-argument _construct_constraint! is used for two-sided constraints.
function _construct_constraint!(aff::AffExpr, lb::Real, ub::Real)
    offset = aff.constant
    aff.constant = 0.0
    LinearConstraint(aff,lb-offset,ub-offset)
end

_construct_constraint!(aff::Variable, lb::Real, ub::Real) = _construct_constraint!(convert(AffExpr,v),lb,ub)

_construct_constraint!(q::QuadExpr, lb, ub) = error("Two-sided quadratic constraints not supported. (Try @addNLConstraint instead.)")

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
        constraintcall = :($addconstr($m, _construct_constraint!($newaff,$(quot(sense)))))
        for kw in kwargs.args
            @assert isexpr(kw, :kw)
            push!(constraintcall.args, esc(kw))
        end
        code = quote
            q = AffExpr()
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
        if VERSION < v"0.4-"
            newlb = esc(x.args[1])
            parselb = nothing
            newub = esc(x.args[5])
            parseub = nothing
        else
            newlb, parselb = parseExprToplevel(x.args[1],:lb)
            newub, parseub = parseExprToplevel(x.args[5],:ub)
        end
        constraintcall = :($addconstr($m, _construct_constraint!($newaff,$newlb,$newub)))
        for kw in kwargs.args
            @assert isexpr(kw, :kw)
            push!(constraintcall.args, esc(kw))
        end
        code = quote
            aff = AffExpr()
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
            c = _construct_constraint!(newaff,$(quot(sense)))
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
        error("Code block passed as constraint. Perhaps you meant to use @LinearConstraints instead?")
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
            q = _construct_constraint!(newaff,$(quot(sense)))
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

for (mac,sym) in [(:LinearConstraints, symbol("@LinearConstraint")),
                  (:QuadConstraints,   symbol("@QuadConstraint"))]
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
        q = AffExpr()
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
    if VERSION <= v"0.4-"
        code = quote
            q = AffExpr()
            $parsecode
        end
    else
        code = quote
            q = 0.0
            $parsecode
        end
    end
    if isa(c,Expr)
        code = quote
            $code
            isa($newaff,AffExpr) || error("Three argument form of @defExpr only supports linear expressions")
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

macro defVar(args...)
    length(args) <= 1 &&
        error("in @defVar ($var): expected model as first argument, then variable information.")
    m = args[1]
    x = args[2]
    extra = vcat(args[3:end]...)
    m = esc(m)

    t = :Cont
    gottype = 0
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
            gottype = 1
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

    if (obj != nothing || inconstraints != nothing || coefficients != nothing) &&
        (obj == nothing || inconstraints == nothing || coefficients == nothing)
        error("in @defVar ($var): Must provide 'objective', 'inconstraints', and 'coefficients' arguments all together for column-wise modeling")
    end


    # Determine variable type (if present).
    # Types: default is continuous (reals)
    if length(extra) > 0
        if t == :Fixed
            error("in @defVar ($var): unexpected extra arguments when declaring a fixed variable")
        end
        if extra[1] in [:Bin, :Int, :SemiCont, :SemiInt]
            gottype = 1
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

        # Handle the old column generation functionality
        if length(extra) - gottype == 3
            warn("in @defVar ($var): the syntax for column-wise modeling has changed. Use @defVar($(m.args[1]), $var, objective=$(extra[1+gottype]), inconstraints=$(extra[2+gottype]), coefficients=$(extra[3+gottype]))")
            objcoef = esc(extra[1+gottype])
            cols    = esc(extra[2+gottype])
            coeffs  = esc(extra[3+gottype])
            return assert_validmodel(m, quote
                $(esc(var)) = Variable($m,$lb,$ub,$(quot(t)),$objcoef,$cols,$coeffs,$(string(var)),$value)
                nothing
            end)
        end

        gottype == 0 &&
            error("in @defVar ($var): syntax error")
    end

    # Handle the column generation functionality
    if coefficients != nothing
        !isa(var,Symbol) &&
        error("in @defVar ($var): can only create one variable at a time when adding to existing constraints.")

        return assert_validmodel(m, quote
            $(esc(var)) = Variable($m,$lb,$ub,$(quot(t)),$obj,$inconstraints,$coefficients,$(string(var)),$value)
            nothing
        end)
    end


    if isa(var,Symbol)
        # Easy case - a single variable
        return assert_validmodel(m, quote
            $(esc(var)) = Variable($m,$lb,$ub,$(quot(t)),$(string(var)),$value)
            registervar($m, $(quot(var)), $(esc(var)))
        end)
    end
    isa(var,Expr) || error("in @defVar: expected $var to be a variable name")

    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(var)
    code = :( $(refcall) = Variable($m, $lb, $ub, $(quot(t)), "", $value) )
    looped = getloopedcode(var, code, condition, idxvars, idxsets, idxpairs, :Variable)
    varname = esc(getname(var))
    return assert_validmodel(m, quote
        $looped
        push!($(m).dictList, $varname)
        registervar($m, $(quot(getname(var))), $varname)
    end)
end

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
        $m.obj = QuadExpr()
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
