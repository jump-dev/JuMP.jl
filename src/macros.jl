#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta

issum(s::Symbol) = (s == :sum) || (s == :∑) || (s == :Σ)
isprod(s::Symbol) = (s == :prod) || (s == :∏)

function error_curly(x)
    Base.error("The curly syntax (sum{},prod{},norm2{}) is no longer supported. Expression: $x.")
end

include("parseexpr.jl")

function buildrefsets(expr::Expr, cname)
    c = copy(expr)
    idxvars = Any[]
    idxsets = Any[]
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
            end
        end
        if !parse_done # No index variable specified
            idxvar = gensym()
            idxset = esc(s)
        end
        push!(idxvars, idxvar)
        push!(idxsets, idxset)
        push!(refcall.args, esc(idxvar))
    end
    return refcall, idxvars, idxsets, condition
end

buildrefsets(c, cname)  = (cname, Any[], Any[], :())

"""
    JuMP.buildrefsets(expr::Expr)

Helper function for macros to construct container objects. Takes an `Expr` that specifies the container, e.g. `:(x[i=1:3,[:red,:blue]],k=S; i+k <= 6)`, and returns:

    1. `refcall`: Expr to reference a particular element in the container, e.g. `:(x[i,red,s])`
    2. `idxvars`: Names for the index variables, e.g. `[:i, gensym(), :k]`
    3. `idxsets`: Sets used for indexing, e.g. `[1:3, [:red,:blue], S]`
    4. `condition`: Expr containing any conditional imposed on indexing, or `:()` if none is present
"""
buildrefsets(c) = buildrefsets(c, getname(c))

"""
    JuMP.getloopedcode(varname, code, condition, idxvars, idxsets, sym, requestedcontainer::Symbol; lowertri=false)

Helper function for macros to transform expression objects containing kernel code, index sets, conditionals, etc. to an expression that performs the desired loops that iterate over the kernel code. Arguments to the function are:

    1. `varname`: name and appropriate indexing sets (if any) for container that is assigned to in the kernel code, e.g. `:myvar` or `:(x[i=1:3,[:red,:blue]])`
    2. `code`: `Expr` containing kernel code
    3. `condition`: `Expr` that is evaluated immediately before kernel code in each iteration. If none, pass `:()`.
    4. `idxvars`: Names for the index variables for each loop, e.g. `[:i, gensym(), :k]`
    5. `idxsets`: Sets used to define iteration for each loop, e.g. `[1:3, [:red,:blue], S]`
    6. `sym`: A `Symbol`/`Expr` containing the element type of the container that is being iterated over, e.g. `:AffExpr` or `:Variable`
    7. `requestedcontainer`: Argument that is passed through to `generatedcontainer`. Either `:Auto`, `:Array`, `:JuMPArray`, or `:Dict`.
    8. `lowertri`: `Bool` keyword argument that is `true` if the iteration is over a cartesian array and should only iterate over the lower triangular entries, filling upper triangular entries with copies, e.g. `x[1,3] === x[3,1]`, and `false` otherwise.
"""
function getloopedcode(varname, code, condition, idxvars, idxsets, sym, requestedcontainer::Symbol; lowertri=false)

    # if we don't have indexing, just return to avoid allocating stuff
    if isempty(idxsets)
        return code
    end

    hascond = (condition != :())

    requestedcontainer in [:Auto, :Array, :JuMPArray, :Dict] || return :(error("Invalid container type $container. Must be Auto, Array, JuMPArray, or Dict."))

    if hascond
        if requestedcontainer == :Auto
            requestedcontainer = :Dict
        elseif requestedcontainer == :Array || requestedcontainer == :JuMPArray
            return :(error("Requested container type is incompatible with conditional indexing. Use :Dict or :Auto instead."))
        end
    end
    containercode, autoduplicatecheck = generatecontainer(sym, idxvars, idxsets, requestedcontainer)

    if lowertri
        @assert !hascond
        @assert length(idxvars)  == 2
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
        if !autoduplicatecheck # we need to check for duplicate keys in the index set
            if length(idxvars) > 1
                keytuple = Expr(:tuple, esc.(idxvars)...)
            else
                keytuple = esc(idxvars[1])
            end
            code = quote
                if haskey($varname, $keytuple)
                    error(string("Repeated index ", $keytuple,". Index sets must have unique elements."))
                end
                $code
            end
        end
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


    return quote
        $varname = $containercode
        $code
        nothing
    end
end

# TODO: Remove all localvar calls for Julia 0.7. The scope of loop variables
# has changed to match the behavior we enforce here.
 localvar(x::Symbol) = _localvar(x)
localvar(x::Expr) = Expr(:block, _localvar(x)...)
_localvar(x::Symbol) = :(local $(esc(x)))
function _localvar(x::Expr)
    @assert x.head in (:escape, :tuple)
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

"""
    extract_kwargs(args)

Process the arguments to a macro, separating out the keyword arguments.
Return a tuple of (flat_arguments, keyword arguments, and requestedcontainer),
where `requestedcontainer` is a symbol to be passed to `getloopedcode`.
"""
function extract_kwargs(args)
    kwargs = filter(x -> isexpr(x, :(=)) && x.args[1] != :container , collect(args))
    flat_args = filter(x->!isexpr(x, :(=)), collect(args))
    requestedcontainer = :Auto
    for kw in args
        if isexpr(kw, :(=)) && kw.args[1] == :container
            requestedcontainer = kw.args[2]
        end
    end
    return flat_args, kwargs, requestedcontainer
end

function addkwargs!(call, kwargs)
    for kw in kwargs
        @assert isexpr(kw, :(=))
        push!(call.args, esc(Expr(:kw, kw.args...)))
    end
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

# two-argument constructconstraint! is used for one-sided constraints.
# Right-hand side is zero.
function sense_to_set(sense::Symbol)
    if sense == :(<=)
        return MOI.LessThan(0.0)
    elseif sense == :(>=)
        return MOI.GreaterThan(0.0)
    else
        @assert sense == :(==)
        return MOI.EqualTo(0.0)
    end
end
const ScalarPolyhedralSets = Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo,MOI.Interval}

constructconstraint!(v::Variable, set::MOI.AbstractScalarSet) = SingleVariableConstraint(v, set)
constructconstraint!(v::Vector{Variable}, set::MOI.AbstractVectorSet) = VectorOfVariablesConstraint(v, set)

function constructconstraint!(aff::AffExpr, set::S) where S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo}
    offset = aff.constant
    aff.constant = 0.0
    return AffExprConstraint(aff, S(MOIU.getconstant(set)-offset))
end

# function constructconstraint!(aff::AffExpr, lb, ub)
#     offset = aff.constant
#     aff.constant = 0.0
#     AffExprConstraint(aff, lb-offset, ub-offset)
# end

constructconstraint!(x::AbstractArray, set::MOI.AbstractScalarSet) = error("Unexpected vector in scalar constraint. Did you mean to use the dot comparison operators like .==, .<=, and .>= instead?")
constructconstraint!(x::Vector{AffExpr}, set::MOI.AbstractVectorSet) = VectorAffExprConstraint(x, set)

function constructconstraint!(quad::QuadExpr, set::S) where S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo}
    offset = quad.aff.constant
    quad.aff.constant = 0.0
    return QuadExprConstraint(quad, S(MOIU.getconstant(set)-offset))
end
#constructconstraint!(x::Vector{QuadExpr}, set::MOI.AbstractVectorSet) = VectorQuadExprConstraint(x, set)


# _vectorize_like(x::Number, y::AbstractArray{AffExpr}) = (ret = similar(y, typeof(x)); fill!(ret, x))
# function _vectorize_like{R<:Number}(x::AbstractArray{R}, y::AbstractArray{AffExpr})
#     for i in 1:max(ndims(x),ndims(y))
#         _size(x,i) == _size(y,i) || error("Unequal sizes for ranged constraint")
#     end
#     x
# end
#
# function constructconstraint!(x::AbstractArray{AffExpr}, lb, ub)
#     LB = _vectorize_like(lb,x)
#     UB = _vectorize_like(ub,x)
#     ret = similar(x, AffExprConstraint)
#     map!(ret, eachindex(ret)) do i
#         constructconstraint!(x[i], LB[i], UB[i])
#     end
# end

# three-argument constructconstraint! is used for two-sided constraints.
function constructconstraint!(aff::AffExpr, lb::Real, ub::Real)
    offset = aff.constant
    aff.constant = 0.0
    AffExprConstraint(aff,MOI.Interval(lb-offset,ub-offset))
end

# constructconstraint!(aff::Variable, lb::Real, ub::Real) = constructconstraint!(convert(AffExpr,v),lb,ub)

constructconstraint!(q::QuadExpr, lb, ub) = error("Two-sided quadratic constraints not supported. (Try @NLconstraint instead.)")

constraint_error(args, str) = error("In @constraint($(join(args,","))): ", str)

# TODO: update 3-argument @constraint macro to pass through names like @variable

"""
    @constraint(m::Model, con)

add linear or quadratic constraints.

    @constraint(m::Model, ref, con)

add groups of linear or quadratic constraints.

"""
macro constraint(args...)
    _error(str) = constraint_error(args, str)

    args, kwargs, requestedcontainer = extract_kwargs(args)

    if length(args) < 2
        if length(kwargs) > 0
            _error("Not enough positional arguments")
        else
            _error("Not enough arguments")
        end
    end
    m = args[1]
    x = args[2]
    extra = args[3:end]

    m = esc(m)
    # Two formats:
    # - @constraint(m, a*x <= 5)
    # - @constraint(m, myref[a=1:5], a*x <= 5)
    length(extra) > 1 && _error("Too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : gensym()
    x = length(extra) == 1 ? extra[1] : x

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat) || length(extra) != 1
    variable = gensym()
    quotvarname = quot(getname(c))
    escvarname  = anonvar ? variable : esc(getname(c))
    basename = anonvar ? "" : string(getname(c))
    # TODO: support the basename keyword argument

    if isa(x, Symbol)
        _error("Incomplete constraint specification $x. Are you missing a comparison (<=, >=, or ==)?")
    end

    (x.head == :block) &&
        _error("Code block passed as constraint. Perhaps you meant to use @constraints instead?")

    # Strategy: build up the code for non-macro addconstraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs
    refcall, idxvars, idxsets, condition = buildrefsets(c, variable)
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
            constraintcall = :(addconstraint($m, constructconstraint!($newaff,$(esc(x.args[3]))), $(namecall(basename, idxvars))))
        else
            # Simple comparison - move everything to the LHS
            @assert length(x.args) == 3
            (sense,vectorized) = _canonicalize_sense(x.args[1])
            set = sense_to_set(sense)
            lhs = :($(x.args[2]) - $(x.args[3]))
            newaff, parsecode = parseExprToplevel(lhs, :q)
            # `set` is an MOI.AbstractScalarSet, if `newaff` is not scalar, vectorized should be true.
            # Otherwise, `constructconstraint!(::AbstractArray, ::MOI.AbstractScalarSet)` throws an helpful error
            if vectorized
                # TODO: Pass through names here.
                constraintcall = :(addconstraint.($m, constructconstraint!.($newaff,$set)))
            else
                constraintcall = :(addconstraint($m, constructconstraint!($newaff,$set), $(namecall(basename, idxvars))))
            end
        end
        addkwargs!(constraintcall, kwargs)
        code = quote
            q = Val{false}()
            $parsecode
            $(refcall) = $constraintcall
        end
    elseif isexpr(x, :comparison)
        # Ranged row
        (lsign,lvectorized) = _canonicalize_sense(x.args[2])
        (rsign,rvectorized) = _canonicalize_sense(x.args[4])
        if (lsign != :(<=)) || (rsign != :(<=))
            _error("Only two-sided rows of the form lb <= expr <= ub are supported.")
        end
        ((vectorized = lvectorized) == rvectorized) || _error("Signs are inconsistently vectorized")
        x_str = string(x)
        lb_str = string(x.args[1])
        ub_str = string(x.args[5])
        newaff, parsecode = parseExprToplevel(x.args[3],:aff)

        newlb, parselb = parseExprToplevel(x.args[1],:lb)
        newub, parseub = parseExprToplevel(x.args[5],:ub)

        if lvectorized
            # TODO: Pass through names here.
            constraintcall = :(addconstraint.($m, constructconstraint!.($newaff,$newlb,$newub)))
        else
            constraintcall = :(addconstraint($m, constructconstraint!($newaff,$newlb,$newub), $(namecall(basename, idxvars))))
        end
        addkwargs!(constraintcall, kwargs)
        code = quote
            aff = Val{false}()
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
                    _error(string("Expected ",$lb_str," to be a ", CoefType, "."))
                end
                try
                    ubval = convert(CoefType, $newub)
                catch
                    _error(string("Expected ",$ub_str," to be a ", CoefType, "."))
                end
            end
        end
        code = quote
            $code
            $(refcall) = $constraintcall
        end
    else
        # Unknown
        _error(string("Constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2\n" * "       lb <= expr <= ub"))
    end
    return assert_validmodel(m, quote
        $(getloopedcode(variable, code, condition, idxvars, idxsets, :ConstraintRef, requestedcontainer))
        $(if anonvar
            variable
        else
            quote
                registercon($m, $quotvarname, $variable)
                $escvarname = $variable
            end
        end)
    end)
end


"""
    @SDconstraint(m, x)

Adds a semidefinite constraint to the `Model m`. The expression `x` must be a square, two-dimensional array.
"""
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
    assert_validmodel(m, quote
        q = Val{false}()
        $parsecode
        addconstraint($m, constructconstraint!($newaff, PSDCone()))
    end)
end


# """
#     @LinearConstraint(x)
#
# Constructs a `LinearConstraint` instance efficiently by parsing the `x`. The same as `@constraint`, except it does not attach the constraint to any model.
# """
# macro LinearConstraint(x)
#     (x.head == :block) &&
#         error("Code block passed as constraint. Perhaps you meant to use @LinearConstraints instead?")
#
#     if isexpr(x, :call) && length(x.args) == 3
#         (sense,vectorized) = _canonicalize_sense(x.args[1])
#         # Simple comparison - move everything to the LHS
#         vectorized &&
#             error("in @LinearConstraint ($(string(x))): Cannot add vectorized constraints")
#         lhs = :($(x.args[2]) - $(x.args[3]))
#         return quote
#             newaff = @Expression($(esc(lhs)))
#             c = constructconstraint!(newaff,$(quot(sense)))
#             isa(c, LinearConstraint) ||
#                 error("Constraint in @LinearConstraint is really a $(typeof(c))")
#             c
#         end
#     elseif isexpr(x, :comparison)
#         # Ranged row
#         (lsense,lvectorized) = _canonicalize_sense(x.args[2])
#         (rsense,rvectorized) = _canonicalize_sense(x.args[4])
#         if (lsense != :<=) || (rsense != :<=)
#             error("in @constraint ($(string(x))): only ranged rows of the form lb <= expr <= ub are supported.")
#         end
#         (lvectorized || rvectorized) &&
#             error("in @LinearConstraint ($(string(x))): Cannot add vectorized constraints")
#         lb = x.args[1]
#         ub = x.args[5]
#         return quote
#             if !isa($(esc(lb)),Number)
#                 error(string("in @LinearConstraint (",$(string(x)),"): expected ",$(string(lb))," to be a number."))
#             elseif !isa($(esc(ub)),Number)
#                 error(string("in @LinearConstraint (",$(string(x)),"): expected ",$(string(ub))," to be a number."))
#             end
#             newaff = @Expression($(esc(x.args[3])))
#             offset = newaff.constant
#             newaff.constant = 0.0
#             isa(newaff,AffExpr) || error("Ranged quadratic constraints are not allowed")
#             LinearConstraint(newaff,$(esc(lb))-offset,$(esc(ub))-offset)
#         end
#     else
#         # Unknown
#         error("in @LinearConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
#               "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
#               "       expr1 == expr2\n" * "       lb <= expr <= ub")
#     end
# end

# """
#     @QuadConstraint(x)
#
# Constructs a `QuadConstraint` instance efficiently by parsing the `x`. The same as `@constraint`, except it does not attach the constraint to any model.
# """
# macro QuadConstraint(x)
#     (x.head == :block) &&
#         error("Code block passed as constraint. Perhaps you meant to use @QuadConstraints instead?")
#
#     if isexpr(x, :call) && length(x.args) == 3
#         (sense,vectorized) = _canonicalize_sense(x.args[1])
#         # Simple comparison - move everything to the LHS
#         vectorized &&
#             error("in @QuadConstraint ($(string(x))): Cannot add vectorized constraints")
#         lhs = :($(x.args[2]) - $(x.args[3]))
#         return quote
#             newaff = @Expression($(esc(lhs)))
#             q = constructconstraint!(newaff,$(quot(sense)))
#             isa(q, QuadConstraint) || error("Constraint in @QuadConstraint is really a $(typeof(q))")
#             q
#         end
#     elseif isexpr(x, :comparison)
#         error("Ranged quadratic constraints are not allowed")
#     else
#         # Unknown
#         error("in @QuadConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
#               "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
#               "       expr1 == expr2")
#     end
# end

# macro SOCConstraint(x)
#     (x.head == :block) &&
#         error("Code block passed as constraint. Perhaps you meant to use @SOCConstraints instead?")
#
#     if isexpr(x, :call) && length(x.args) == 3
#         (sense,vectorized) = _canonicalize_sense(x.args[1])
#         # Simple comparison - move everything to the LHS
#         vectorized &&
#             error("in @SOCConstraint ($(string(x))): Cannot add vectorized constraints")
#         lhs = :($(x.args[2]) - $(x.args[3]))
#         return quote
#             newaff = @Expression($(esc(lhs)))
#             q = constructconstraint!(newaff,$(quot(sense)))
#             isa(q, SOCConstraint) || error("Constraint in @SOCConstraint is really a $(typeof(q))")
#             q
#         end
#     elseif isexpr(x, :comparison)
#         error("Ranged second-order cone constraints are not allowed")
#     else
#         # Unknown
#         error("in @SOCConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
#               "       expr1 <= expr2\n" * "       expr1 >= expr2")
#     end
# end

# for (mac,sym) in [(:LinearConstraints, Symbol("@LinearConstraint")),
#                   (:QuadConstraints,   Symbol("@QuadConstraint")),
#                   (:SOCConstraints,    Symbol("@SOCConstraint"))]
#     @eval begin
#         macro $mac(x)
#             x.head == :block || error(string("Invalid syntax for @", $(string(mac))))
#             @assert x.args[1].head == :line
#             code = Expr(:vect)
#             for it in x.args
#                 if it.head == :line
#                     # do nothing
#                 elseif it.head == :comparison || (it.head == :call && it.args[1] in (:<=,:≤,:>=,:≥,:(==))) # regular constraint
#                     push!(code.args, Expr(:macrocall, $sym, esc(it)))
#                 elseif it.head == :tuple # constraint ref
#                     if all([isexpr(arg,:comparison) for arg in it.args]...)
#                         # the user probably had trailing commas at end of lines, e.g.
#                         # @LinearConstraints(m, begin
#                         #     x <= 1,
#                         #     x >= 1
#                         # end)
#                         error(string("Invalid syntax in @", $(string(mac)), ". Do you have commas at the end of a line specifying a constraint?"))
#                     end
#                     error("@", string($(string(mac)), " does not currently support the two argument syntax for specifying groups of constraints in one line."))
#                 else
#                     error("Unexpected constraint expression $it")
#                 end
#             end
#             return code
#         end
#     end
# end

for (mac,sym) in [(:constraints,  Symbol("@constraint")),
                  (:NLconstraints,Symbol("@NLconstraint")),
                  (:SDconstraints,Symbol("@SDconstraint")),
                  (:variables,Symbol("@variable")),
                  (:expressions, Symbol("@expression")),
                  (:NLexpressions, Symbol("@NLexpression"))]
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
                        if isexpr(ex, :(=)) && VERSION < v"0.6.0-dev.1934"
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


# Doc strings for the auto-generated macro pluralizations
@doc """
    @constraints(m, args...)

adds groups of constraints at once, in the same fashion as @constraint. The model must be the first argument, and multiple constraints can be added on multiple lines wrapped in a `begin ... end` block. For example:

    @constraints(m, begin
      x >= 1
      y - w <= 2
      sum_to_one[i=1:3], z[i] + y == 1
    end)
""" :(@constraints)

@doc """
    @LinearConstraints(m, args...)

Constructs a vector of `LinearConstraint` objects. Similar to `@LinearConstraint`, except it accepts multiple constraints as input as long as they are separated by newlines.
""" :(@LinearConstraints)

@doc """
    @QuadConstraints(m, args...)

Constructs a vector of `QuadConstraint` objects. Similar to `@QuadConstraint`, except it accepts multiple constraints as input as long as they are separated by newlines.
""" :(@QuadConstraints)






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
        q = Val{false}()
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
        q = Val{false}()
        $parsecode
        $newaff
    end
end


"""
    @expression(args...)

efficiently builds a linear, quadratic, or second-order cone expression but does not add to model immediately. Instead, returns the expression which can then be inserted in other constraints. For example:

```julia
@expression(m, shared, sum(i*x[i] for i=1:5))
@constraint(m, shared + y >= 5)
@constraint(m, shared + z <= 10)
```

The `ref` accepts index sets in the same way as `@variable`, and those indices can be used in the construction of the expressions:

```julia
@expression(m, expr[i=1:3], i*sum(x[j] for j=1:3))
```

Anonymous syntax is also supported:

```julia
expr = @expression(m, [i=1:3], i*sum(x[j] for j=1:3))
```
"""
macro expression(args...)

    args, kwargs, requestedcontainer = extract_kwargs(args)
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
    length(kwargs) == 0 || error("@expression: unrecognized keyword argument")

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    variable = gensym()
    escvarname  = anonvar ? variable : esc(getname(c))

    refcall, idxvars, idxsets, condition = buildrefsets(c, variable)
    newaff, parsecode = parseExprToplevel(x, :q)
    code = quote
        q = Val{false}()
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
    code = getloopedcode(variable, code, condition, idxvars, idxsets, :AffExpr, requestedcontainer)
    # don't do anything with the model, but check that it's valid anyway
    return assert_validmodel(m, quote
        $code
        $(anonvar ? variable : :($escvarname = $variable))
    end)
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
function constructvariable!(m::Model, _error::Function, haslb::Bool, lowerbound::Number, hasub::Bool, upperbound::Number,
                            hasfix::Bool, fixedvalue::Number, binary::Bool, integer::Bool, name::String,
                            hasstart::Bool, start::Number; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    v = Variable(m)
    if haslb
        setlowerbound(v, lowerbound)
    end
    if hasub
        setupperbound(v, upperbound)
    end
    if hasfix
        fix(v, fixedvalue)
    end
    if binary
        setbinary(v)
    end
    if integer
        setinteger(v)
    end
    if hasstart
        setstartvalue(v, start)
    end
    if name != EMPTYSTRING
        setname(v, name)
    end
    return v
end

const EMPTYSTRING = ""

variable_error(args, str) = error("In @variable($(join(args,","))): ", str)

# Given a basename and idxvars, returns an expression that constructs the name
# of the object. For use within macros only.
function namecall(basename, idxvars)
    if length(idxvars) == 0 || basename == ""
        return basename
    end
    ex = Expr(:call,:string,basename,"[")
    for i in 1:length(idxvars)
        push!(ex.args, esc(idxvars[i]))
        i < length(idxvars) && push!(ex.args,",")
    end
    push!(ex.args,"]")
    return ex
end

# @variable(m, expr, extra...; kwargs...)
# where `extra` is a list of extra positional arguments and `kwargs` is a list of keyword arguments.
#
# It creates a new variable (resp. a container of new variables) belonging to the model `m` using `constructvariable!` to create the variable (resp. each variable of the container).
# The following modifications will be made to the arguments before they are passed to `constructvariable!`:
# * The `expr` argument will not be passed but the expression will be parsed to determine the kind of container needed (if one is needed) and
#   additional information that will alter what is passed with the keywords `lowerbound`, `upperbound`, `basename`, `start`, `binary`, and `integer`.
# * The `PSD` and `Symmetric` positional arguments in `extra` will not be passed to `constructvariable!`. Instead,
#    * the `Symmetric` argument will check that the container is symmetric and only allocate one variable for each pair of non-diagonal entries.
#    * the `PSD` argument will do the same as `Symmetric` but in addition it will specify that the variables created belongs to the PSD cone in the `varCones` field of the model.
#   Moreover, Int and Bin are special keywords that are equivalent to `integer=true` and `binary=true`.
# * The keyword arguments start, basename, lowerbound, upperbound, binary, and integer category may not be passed as is to
#   `constructvariable!` since they may be altered by the parsing of `expr` and we may need to pass it pointwise if it is a container since
#   `constructvariable!` is called separately for each variable of the container. Moreover it will be passed as positional argument to `constructvariable!`.
# * A custom error function is passed as positional argument to print the full @variable call before the error message.
#
# Examples (... is the custom error function):
# * `@variable(m, x >= 0)` is equivalent to `x = constructvariable!(m, msg -> error("In @variable(m, x >= 0): ", msg), true, 0, false, NaN, false, NaN, false, false, "x", false, NaN)
# * `@variable(m, x[1:N,1:N], Symmetric, Poly(X))` is equivalent to
#   ```
#   x = Matrix{...}(N, N)
#   for i in 1:N
#       for j in 1:N
#           x[i,j] = x[j,i] = constructvariable!(m, Poly(X), msg -> error("In @variable(m, x[1:N,1:N], Symmetric, Poly(X)): ", msg), false, NaN, false, NaN, false, NaN, false, false, "", false, NaN)
#       end
#   end
#   ```
macro variable(args...)
    _error(str) = variable_error(args, str)

    m = esc(args[1])

    extra, kwargs, requestedcontainer = extract_kwargs(args[2:end])

    # if there is only a single non-keyword argument, this is an anonymous
    # variable spec and the one non-kwarg is the model
    if length(extra) == 0
        x = gensym()
        anon_singleton = true
    else
        x = shift!(extra)
        if x in [:Int,:Bin,:PSD]
            _error("Ambiguous variable name $x detected. Use the \"category\" keyword argument to specify a category for an anonymous variable.")
        end
        anon_singleton = false
    end

    haslb = false
    hasub = false
    hasfix = false
    hasstart = false
    fixedvalue = NaN
    var = x
    lb = NaN
    ub = NaN
    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    explicit_comparison = false
    if isexpr(x,:comparison) # two-sided
        explicit_comparison = true
        hasub = true
        haslb = true
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
            fixedvalue = esc(x.args[3])
            hasfix = true
        else
            # Its a comparsion, but not using <= ... <=
            _error("Unexpected syntax $(string(x)).")
        end
    end

    anonvar = isexpr(var, :vect) || isexpr(var, :vcat) || anon_singleton
    anonvar && explicit_comparison && error("Cannot use explicit bounds via >=, <= with an anonymous variable")
    variable = gensym()
    quotvarname = anonvar ? :(:__anon__) : quot(getname(var))
    escvarname  = anonvar ? variable     : esc(getname(var))
    # TODO: Should we generate non-empty default names for variables?
    basename = anonvar ? "" : string(getname(var))

    if !isa(getname(var),Symbol) && !anonvar
        Base.error("Expression $(getname(var)) should not be used as a variable name. Use the \"anonymous\" syntax $(getname(var)) = @variable(m, ...) instead.")
    end

    # process keyword arguments
    value = NaN
    obj = nothing
    binary = false
    integer = false
    extra_kwargs = []
    for ex in kwargs
        kwarg = ex.args[1]
        if kwarg == :start
            hasstart = true
            value = esc(ex.args[2])
        elseif kwarg == :basename
            basename = esc(ex.args[2])
        elseif kwarg == :lowerbound
            haslb && _error("Cannot specify variable lowerbound twice")
            lb = esc_nonconstant(ex.args[2])
            haslb = true
        elseif kwarg == :upperbound
            hasub && _error("Cannot specify variable upperbound twice")
            ub = esc_nonconstant(ex.args[2])
            hasub = true
        elseif kwarg == :integer
            integer = esc_nonconstant(ex.args[2])
        elseif kwarg == :binary
            binary = esc_nonconstant(ex.args[2])
        else
            push!(extra_kwargs, ex)
        end
    end

    sdp = any(t -> (t == :PSD), extra)
    symmetric = (sdp || any(t -> (t == :Symmetric), extra))
    extra = filter(x -> (x != :PSD && x != :Symmetric), extra) # filter out PSD and sym tag
    for ex in extra
        if ex == :Int
            if integer != false
                _error("'Int' and 'integer' keyword argument cannot both be specified.")
            end
            integer = true
        elseif ex == :Bin
            if binary != false
                _error("'Bin' and 'binary' keyword argument cannot both be specified.")
            end
            binary = true
        end
    end
    extra = esc.(filter(ex -> !(ex in [:Int,:Bin]), extra))

    if isa(var,Symbol)
        # Easy case - a single variable
        sdp && _error("Cannot add a semidefinite scalar variable")
        variablecall = :( constructvariable!($m, $(extra...), $_error, $haslb, $lb, $hasub, $ub, $hasfix, $fixedvalue, $binary, $integer, $basename, $hasstart, $value) )
        addkwargs!(variablecall, extra_kwargs)
        code = :($variable = $variablecall)
        if !anonvar
            code = quote
                $code
                registervar($m, $quotvarname, $variable)
                $escvarname = $variable
            end
        end
        return assert_validmodel(m, code)
    end
    isa(var,Expr) || _error("Expected $var to be a variable name")

    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall, idxvars, idxsets, condition = buildrefsets(var, variable)
    clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? () : idxsets[i])

    # Code to be used to create each variable of the container.
    variablecall = :( constructvariable!($m, $(extra...), $_error, $haslb, $lb, $hasub, $ub, $hasfix, $fixedvalue, $binary, $integer, $(namecall(basename, idxvars)), $hasstart, $value) )
    addkwargs!(variablecall, extra_kwargs)
    code = :( $(refcall) = $variablecall )
    # Determine the return type of constructvariable!. This is needed to create the container holding them.
    vartype = :( variabletype($m, $(extra...)) )

    if symmetric
        # Sanity checks on PSD input stuff
        condition == :() ||
            _error("Cannot have conditional indexing for PSD variables")
        length(idxvars) == length(idxsets) == 2 ||
            _error("PSD variables must be 2-dimensional")
        !symmetric || (length(idxvars) == length(idxsets) == 2) ||
            _error("Symmetric variables must be 2-dimensional")
        hasdependentsets(idxvars, idxsets) &&
            _error("Cannot have index dependencies in symmetric variables")
        for _rng in idxsets
            isexpr(_rng, :escape) ||
                _error("Internal error 1")
            rng = _rng.args[1] # undo escaping
            (isexpr(rng,:(:)) && rng.args[1] == 1 && length(rng.args) == 2) ||
                _error("Index sets for PSD variables must be ranges of the form 1:N")
        end

        if haslb || hasub
            _error("Semidefinite or symmetric variables cannot be provided bounds")
        end
        return assert_validmodel(m, quote
            $(esc(idxsets[1].args[1].args[2])) == $(esc(idxsets[2].args[1].args[2])) || error("Cannot construct symmetric variables with nonsquare dimensions")
            $(getloopedcode(variable, code, condition, idxvars, idxsets, vartype, requestedcontainer; lowertri=symmetric))
            $(if sdp
                quote
                    JuMP.addconstraint($m, JuMP.constructconstraint!(Symmetric($variable), JuMP.PSDCone()))
                end
            end)
            !$anonvar && registervar($m, $quotvarname, $variable)
            $(anonvar ? variable : :($escvarname = Symmetric($variable)))
        end)
    else
        return assert_validmodel(m, quote
            $(getloopedcode(variable, code, condition, idxvars, idxsets, vartype, requestedcontainer))
            !$anonvar && registervar($m, $quotvarname, $variable)
            $(anonvar ? variable : :($escvarname = $variable))
        end)
    end
end

# TODO: replace with a general macro that can construct any container type
# macro constraintref(var)
#     if isa(var,Symbol)
#         # easy case
#         return esc(:(local $var))
#     else
#         if !isexpr(var,:ref)
#             error("Syntax error: Expected $var to be of form var[...]")
#         end
#
#         varname = var.args[1]
#         idxsets = var.args[2:end]
#
#         code = quote
#             $(esc(gendict(varname, :ConstraintRef, idxsets...)))
#             nothing
#         end
#         return code
#     end
# end

macro NLobjective(m, sense, x)
    m = esc(m)
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    return assert_validmodel(m, quote
        ex = @processNLExpr($m, $(esc(x)))
        setobjective($m, $(esc(sense)), ex)
    end)
end

macro NLconstraint(m, x, extra...)
    m = esc(m)
    # Two formats:
    # - @NLconstraint(m, a*x <= 5)
    # - @NLconstraint(m, myref[a=1:5], sin(x^a) <= 5)
    extra, kwargs, requestedcontainer = extract_kwargs(extra)
    (length(extra) > 1 || length(kwargs) > 0) && error("in @NLconstraint: too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : gensym()
    x = length(extra) == 1 ? extra[1] : x

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat) || length(extra) != 1
    variable = gensym()
    quotvarname = anonvar ? :(:__anon__) : quot(getname(c))
    escvarname  = anonvar ? variable : esc(getname(c))

    # Strategy: build up the code for non-macro addconstraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs
    refcall, idxvars, idxsets, condition = buildrefsets(c, variable)
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
            $(refcall) = ConstraintRef($m, NonlinearConstraintIndex(length($m.nlpdata.nlconstr)))
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
            $(refcall) = ConstraintRef($m, NonlinearConstraintIndex(length($m.nlpdata.nlconstr)))
        end
    else
        # Unknown
        error("in @NLconstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2")
    end
    looped = getloopedcode(variable, code, condition, idxvars, idxsets, :(ConstraintRef{Model,NonlinearConstraintIndex}), requestedcontainer)
    return assert_validmodel(m, quote
        initNLP($m)
        $looped
        $(if anonvar
            variable
        else
            quote
                registercon($m, $quotvarname, $variable)
                $escvarname = $variable
            end
        end)
    end)
end

macro NLexpression(args...)
    args, kwargs, requestedcontainer = extract_kwargs(args)
    if length(args) <= 1
        error("in @NLexpression: To few arguments ($(length(args))); must pass the model and nonlinear expression as arguments.")
    elseif length(args) == 2
        m, x = args
        m = esc(m)
        c = gensym()
    elseif length(args) == 3
        m, c, x = args
        m = esc(m)
    end
    if length(args) > 3 || length(kwargs) > 0
        error("in @NLexpression: To many arguments ($(length(args))).")
    end

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    variable = gensym()
    escvarname  = anonvar ? variable : esc(getname(c))

    refcall, idxvars, idxsets, condition = buildrefsets(c, variable)
    code = quote
        $(refcall) = NonlinearExpression($m, @processNLExpr($m, $(esc(x))))
    end
    return assert_validmodel(m, quote
        $(getloopedcode(variable, code, condition, idxvars, idxsets, :NonlinearExpression, requestedcontainer))
        $(anonvar ? variable : :($escvarname = $variable))
    end)
end

# syntax is @NLparameter(m, p[i=1] == 2i)
macro NLparameter(m, ex, extra...)

    extra, kwargs, requestedcontainer = extract_kwargs(extra)
    (length(extra) == 0 && length(kwargs) == 0) || error("in @NLperameter: too many arguments.")
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
    escvarname  = anonvar ? variable : esc(getname(c))

    refcall, idxvars, idxsets, condition = buildrefsets(c, variable)
    code = quote
        $(refcall) = newparameter($m, $(esc(x)))
    end
    return assert_validmodel(m, quote
        $(getloopedcode(variable, code, condition, idxvars, idxsets, :NonlinearParameter, :Auto))
        $(anonvar ? variable : :($escvarname = $variable))
    end)
end
