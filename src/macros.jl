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

include("parse_expr.jl")

function buildrefsets(expr::Expr, cname)
    c = copy(expr)
    idxvars = Any[]
    idxsets = Any[]
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
            parse_done, idxvar, _idxset = try_parse_idx_set(s::Expr)
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
    6. `sym`: A `Symbol`/`Expr` containing the element type of the container that is being iterated over, e.g. `:AffExpr` or `:VariableRef`
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
validmodel(m, name) = error("Expected $name to be a JuMP model, but it has type ", typeof(m))

function assert_validmodel(m, macrocode)
    # assumes m is already escaped
    quote
        validmodel($m, $(quot(m.args[1])))
        $macrocode
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
                            final_variable=variable,
                            model_for_registering=nothing)

Return runs `code` in a local scope which returns the value of `variable`
and then assign `final_variable` to `name`.
If `model_for_registering` is given, `register_object(model, name, variable)` is
called in the generated code.
"""
function macro_assign_and_return(code, variable, name;
                                 final_variable=variable,
                                 model_for_registering=nothing)
    macro_code = macro_return(code, variable)
    return quote
        $variable = $macro_code
        $(if model_for_registering !== nothing
              :(register_object($model_for_registering, $(quot(name)),
                $variable))
          end)
        # This assignment should be in the scope calling the macro
        $(esc(name)) = $final_variable
    end
end

function _check_vectorized(sense::Symbol)
    sense_str = string(sense)
    if sense_str[1] == '.'
        Symbol(sense_str[2:end]), true
    else
        sense, false
    end
end

# two-argument build_constraint is used for one-sided constraints.
# Right-hand side is zero.
sense_to_set(_error::Function, ::Union{Val{:(<=)}, Val{:(≤)}}) = MOI.LessThan(0.0)
sense_to_set(_error::Function, ::Union{Val{:(>=)}, Val{:(≥)}}) = MOI.GreaterThan(0.0)
sense_to_set(_error::Function, ::Val{:(==)}) = MOI.EqualTo(0.0)
sense_to_set(_error::Function, ::Val{S}) where S = _error("Unrecognized sense $S")

function parse_one_operator_constraint(_error::Function, vectorized::Bool, ::Val{:in}, aff, set)
    newaff, parseaff = parseExprToplevel(aff, :q)
    parsecode = :(q = Val{false}(); $parseaff)
    if vectorized
        buildcall = :(build_constraint.($_error, $newaff, Ref($(esc(set)))))
    else
        buildcall = :(build_constraint($_error, $newaff, $(esc(set))))
    end
    parsecode, buildcall
end

function parse_one_operator_constraint(_error::Function, vectorized::Bool, sense::Val, lhs, rhs)
    # Simple comparison - move everything to the LHS
    aff = :($lhs - $rhs)
    set = sense_to_set(_error, sense)
    parse_one_operator_constraint(_error, vectorized, Val(:in), aff, set)
end

function parse_constraint(_error::Function, sense::Symbol, lhs, rhs)
    (sense, vectorized) = _check_vectorized(sense)
    vectorized, parse_one_operator_constraint(_error, vectorized, Val(sense), lhs, rhs)...
end

function parse_ternary_constraint(_error::Function, vectorized::Bool, lb, ::Union{Val{:(<=)}, Val{:(≤)}}, aff, rsign::Union{Val{:(<=)}, Val{:(≤)}}, ub)
    newaff, parseaff = parseExprToplevel(aff, :aff)
    newlb, parselb = parseExprToplevel(lb, :lb)
    newub, parseub = parseExprToplevel(ub, :ub)
    if vectorized
        buildcall = :(build_constraint.($_error, $newaff, $newlb, $newub))
    else
        buildcall = :(build_constraint($_error, $newaff, $newlb, $newub))
    end
    parseaff, parselb, parseub, buildcall
end

function parse_ternary_constraint(_error::Function, vectorized::Bool, ub, ::Union{Val{:(>=)}, Val{:(≥)}}, aff, rsign::Union{Val{:(>=)}, Val{:(≥)}}, lb)
    parse_ternary_constraint(_error, vectorized, lb, Val(:(<=)), aff, Val(:(<=)), ub)
end

function parse_ternary_constraint(_error::Function, args...)
    _error("Only two-sided rows of the form lb <= expr <= ub or ub >= expr >= lb are supported.")
end

function parse_constraint(_error::Function, lb, lsign::Symbol, aff, rsign::Symbol, ub)
    (lsign, lvectorized) = _check_vectorized(lsign)
    (rsign, rvectorized) = _check_vectorized(rsign)
    ((vectorized = lvectorized) == rvectorized) || _error("Signs are inconsistently vectorized")
    parseaff, parselb, parseub, buildcall = parse_ternary_constraint(_error, vectorized, lb, Val(lsign), aff, Val(rsign), ub)
    parsecode = quote
        aff = Val{false}()
        $parseaff
        lb = 0.0
        $parselb
        ub = 0.0
        $parseub
    end
    vectorized, parsecode, buildcall
end

function parse_constraint(_error::Function, args...)
    # Unknown
    _error("Constraints must be in one of the following forms:\n" *
          "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
          "       expr1 == expr2\n" * "       lb <= expr <= ub")
end

function build_constraint(_error::Function, v::AbstractJuMPScalar,
                          set::MOI.AbstractScalarSet)
    return ScalarConstraint(v, set)
end
function build_constraint(_error::Function,
                          expr::Union{GenericAffExpr, GenericQuadExpr},
                          set::MOI.AbstractScalarSet)
    offset = constant(expr)
    add_to_expression!(expr, -offset)
    return ScalarConstraint(expr, MOIU.shift_constant(set, -offset))
end
function build_constraint(_error::Function, α::Number,
                          set::MOI.AbstractScalarSet)
    return build_constraint(_error, convert(AffExpr, α), set)
end

function build_constraint(_error::Function, x::Vector{<:AbstractJuMPScalar},
                          set::MOI.AbstractVectorSet)
    return VectorConstraint(x, set)
end
function build_constraint(_error::Function, x::AbstractArray,
                          set::MOI.AbstractScalarSet)
    return _error("Unexpected vector in scalar constraint. Did you mean to use",
                  " the dot comparison operators like .==, .<=, and .>=",
                  " instead?")
end

# _vectorize_like(x::Number, y::AbstractArray{AffExpr}) = (ret = similar(y, typeof(x)); fill!(ret, x))
# function _vectorize_like{R<:Number}(x::AbstractArray{R}, y::AbstractArray{AffExpr})
#     for i in 1:max(ndims(x),ndims(y))
#         _size(x,i) == _size(y,i) || error("Unequal sizes for ranged constraint")
#     end
#     x
# end
#
# function build_constraint(x::AbstractArray{AffExpr}, lb, ub)
#     LB = _vectorize_like(lb,x)
#     UB = _vectorize_like(ub,x)
#     ret = similar(x, AffExprConstraint)
#     map!(ret, eachindex(ret)) do i
#         build_constraint(x[i], LB[i], UB[i])
#     end
# end

# three-argument build_constraint is used for two-sided constraints.
build_constraint(_error::Function, func::AbstractJuMPScalar, lb::Real, ub::Real) = build_constraint(_error, func, MOI.Interval(lb, ub))

function build_constraint(_error::Function, expr, lb, ub)
    lb isa Number || _error(string("Expected $lb to be a number."))
    ub isa Number || _error(string("Expected $ub to be a number."))
    if lb isa Number && ub isa Number
        _error("Range constraint is not supported for $expr.")
    end
end

# TODO: update 3-argument @constraint macro to pass through names like @variable

"""
    constraint_macro(args, macro_name::Symbol, parsefun::Function)

Returns the code for the macro `@constraint_like args...` of syntax
```julia
@constraint_like con     # Single constraint
@constraint_like ref con # group of constraints
```
where `@constraint_like` is either `@constraint` or `@SDconstraint`.
The expression `con` is parsed by `parsefun` which returns a code that, when
executed, returns an `AbstractConstraint`. This `AbstractConstraint` is passed
to `add_constraint` with the macro keyword arguments (except the `container`
keyword argument which is used to determine the container type).
"""
function constraint_macro(args, macro_name::Symbol, parsefun::Function)
    _error(str...) = macro_error(macro_name, args, str...)

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
    # - @constraint_like(m, a*x <= 5)
    # - @constraint_like(m, myref[a=1:5], a*x <= 5)
    length(extra) > 1 && _error("Too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : gensym()
    x = length(extra) == 1 ? extra[1] : x

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat) || length(extra) != 1
    variable = gensym()
    name = getname(c)
    base_name = anonvar ? "" : string(name)
    # TODO: support the base_name keyword argument

    if isa(x, Symbol)
        _error("Incomplete constraint specification $x. Are you missing a comparison (<=, >=, or ==)?")
    end

    (x.head == :block) &&
        _error("Code block passed as constraint. Perhaps you meant to use @constraints instead?")

    # Strategy: build up the code for add_constraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs
    refcall, idxvars, idxsets, condition = buildrefsets(c, variable)

    vectorized, parsecode, buildcall = parsefun(_error, x.args...)
    if vectorized
        # TODO: Pass through names here.
        constraintcall = :(add_constraint.($m, $buildcall))
    else
        constraintcall = :(add_constraint($m, $buildcall, $(namecall(base_name, idxvars))))
    end
    addkwargs!(constraintcall, kwargs)
    code = quote
        $parsecode
        $(refcall) = $constraintcall
    end

    # Determine the return type of add_constraint. This is needed for JuMP extensions for which this is different than ConstraintRef
    if vectorized
        contype = :( AbstractArray{constraint_type($m)} ) # TODO use a concrete type instead of AbstractArray, see #525, #1310
    else
        contype = :( constraint_type($m) )
    end
    creationcode = getloopedcode(variable, code, condition, idxvars, idxsets, contype, requestedcontainer)

    if anonvar
        # Anonymous constraint, no need to register it in the model-level
        # dictionary nor to assign it to a variable in the user scope.
        # We simply return the constraint reference
        macro_code = macro_return(creationcode, variable)
    else
        # We register the constraint reference to its name and
        # we assign it to a variable in the local scope of this name
        macro_code = macro_assign_and_return(creationcode, variable, name,
                                             model_for_registering = m)
    end
    return assert_validmodel(m, macro_code)
end

# This function needs to be implemented by all `AbstractModel`s
constraint_type(m::Model) = ConstraintRef{typeof(m)}

"""
    @constraint(m::Model, expr)

Add a constraint described by the expression `expr`.

    @constraint(m::Model, ref[i=..., j=..., ...], expr)

Add a group of constraints described by the expression `expr` parametrized by
`i`, `j`, ...

The expression `expr` can either be

* of the form `func in set` constraining the function `func` to belong to the
  set `set`, e.g. `@constraint(m, [1, x-1, y-2] in MOI.SecondOrderCone(3))`
  constrains the norm of `[x-1, y-2]` be less than 1;
* of the form `a sign b`, where `sign` is one of `==`, `≥`, `>=`, `≤` and
  `<=` building the single constraint enforcing the comparison to hold for the
  expression `a` and `b`, e.g. `@constraint(m, x^2 + y^2 == 1)` constrains `x`
  and `y` to lie on the unit circle;
* of the form `a ≤ b ≤ c` or `a ≥ b ≥ c` (where `≤` and `<=` (resp. `≥` and
  `>=`) can be used interchangeably) constraining the paired the expression
  `b` to lie between `a` and `c`;
* of the forms `@constraint(m, a .sign b)` or
  `@constraint(m, a .sign b .sign c)` which broadcast the constraint creation to
  each element of the vectors.

## Note for extending the constraint macro

Each constraint will be created using
`add_constraint(m, build_constraint(_error, func, set))` where
* `_error` is an error function showing the constraint call in addition to the
  error message given as argument,
* `func` is the expression that is constrained
* and `set` is the set in which it is constrained to belong.

For `expr` of the first type (i.e. `@constraint(m, func in set)`), `func` and
`set` are passed unchanged to `build_constraint` but for the other types, they
are determined from the expressions and signs. For instance,
`@constraint(m, x^2 + y^2 == 1)` is transformed into
`add_constraint(m, build_constraint(_error, x^2 + y^2, MOI.EqualTo(1.0)))`.

To extend JuMP to accept new constraints of this form, it is necessary to add
the corresponding methods to `build_constraint`. Note that this will likely mean
that either `func` or `set` will be some custom type, rather than e.g. a
`Symbol`, since we will likely want to dispatch on the type of the function or
set appearing in the constraint.
"""
macro constraint(args...)
    constraint_macro(args, :constraint, parse_constraint)
end

function parseSDconstraint(_error::Function, sense::Symbol, lhs, rhs)
    # Simple comparison - move everything to the LHS
    aff = :()
    if sense == :⪰ || sense == :(≥) || sense == :(>=)
        aff = :($lhs - $rhs)
    elseif sense == :⪯ || sense == :(≤) || sense == :(<=)
        aff = :($rhs - $lhs)
    else
        _error("Invalid sense $sense in SDP constraint")
    end
    vectorized = false
    parsecode, buildcall = parse_one_operator_constraint(_error, false, Val(:in), aff, :(PSDCone()))
    vectorized, parsecode, buildcall
end

function parseSDconstraint(_error::Function, args...)
    _error("Constraints must be in one of the following forms:\n" *
           "       expr1 <= expr2\n" *
           "       expr1 >= expr2")
end

"""
    @SDconstraint(model::Model, expr)

Add a semidefinite constraint described by the expression `expr`.

    @SDconstraint(model::Model, ref[i=..., j=..., ...], expr)

Add a group of semidefinite constraints described by the expression `expr`
parametrized by `i`, `j`, ...

The expression `expr` needs to be of the form `a sign b` where `sign` is `⪰`,
`≥`, `>=`, `⪯`, `≤` or `<=` and `a` and `b` are `square` matrices. It
constrains the matrix `x = a - b` (or `x = b - a` if the sign is `⪯`, `≤` or
`<=`) to be symmetric and positive semidefinite.

If `x` is already symmetric, use `@constraint(model, Symmetric(x) in PSDCone())`
instead to remove the need to compare the corresponding off-diagonal entries
and add equality constraints if they contain different expression.

## Examples

The following constrains the matrix `[x-1 2x-2; -3 x-4]` to be symmetric and
positive semidefinite, that is, it constrains `2x-2` to be equal to `-3` and
constrains all eigenvalues of the matrix to be nonnegative.
```jldoctest SDconstraint
julia> using JuMP

julia> model = Model();

julia> @variable(model, x)
x

julia> a = [x 2x
            0  x];

julia> b = [1 2
            3 4];

julia> @SDconstraint(model, a ⪰ b)
[x - 1, -3, 2 x - 2, x - 4] ∈ MathOptInterface.PositiveSemidefiniteConeSquare(2)
```
In the set `PositiveSemidefiniteConeSquare(2)` in the last output, `Square`
means that the matrix is passed as a square matrix as the corresponding
off-diagonal entries need to be constrained to be equal. A similar set
`PositiveSemidefiniteConeTriangle` exists which only uses the upper triangular
part of the matrix assuming that it is symmetric. As mentioned above, use
`@constraint(model, Symmetric(x) in PSDCone())` to construct such constraint
where `x` is known to be symmetric as in the following example:
```jldoctest SDconstraint
julia> a = [ x 2x
            2x  x];

julia> b = [1 2
            2 4];

julia> @SDconstraint(model, a ⪰ b)
[x - 1, 2 x - 2, 2 x - 2, x - 4] ∈ MathOptInterface.PositiveSemidefiniteConeSquare(2)

julia> using LinearAlgebra # For Symmetric

julia> @constraint(model, Symmetric(a - b) in PSDCone())
[x - 1, 2 x - 2, x - 4] ∈ MathOptInterface.PositiveSemidefiniteConeTriangle(2)
```
As we see in the output of the `@constraint`, only the upper triangular part
of the matrix is passed.
"""
macro SDconstraint(args...)
    constraint_macro(args, :SDconstraint, parseSDconstraint)
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
#             c = build_constraint(newaff,$(quot(sense)))
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
#             q = build_constraint(newaff,$(quot(sense)))
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
#             q = build_constraint(newaff,$(quot(sense)))
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

add_JuMP_prefix(s::Symbol) = Expr(:., JuMP, :($(QuoteNode(s))))

for (mac,sym) in [(:constraints,  Symbol("@constraint")),
                  (:NLconstraints,Symbol("@NLconstraint")),
                  (:SDconstraints,Symbol("@SDconstraint")),
                  (:variables,Symbol("@variable")),
                  (:expressions, Symbol("@expression")),
                  (:NLexpressions, Symbol("@NLexpression"))]
    if VERSION >= v"0.7-"
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
                        # x, (start = 10, lower_bound = 5)
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

"""
    moi_sense(_error::Function, sense)

Return an expression whose value is an `MOI.OptimizationSense` corresponding
to `sense`. Sense is either the symbol `:Min` or `:Max`, corresponding
respectively to `MOI.MinSense` and `MOI.MaxSense` or it is another symbol,
which should be the name of a variable or expression whose value is an
`MOI.OptimizationSense`.
In the last case, the expression throws an error using the `_error`
function in case the value is not an `MOI.OptimizationSense`.
"""
function moi_sense(_error::Function, sense)
    if sense == :Min
        expr = MOI.MinSense
    elseif sense == :Max
        expr = MOI.MaxSense
    else
        # Refers to a variable that holds the sense.
        # TODO: Better document this behavior
        expr = esc(sense)
    end
    return :(throw_error_for_invalid_sense($_error, $expr))
end

function throw_error_for_invalid_sense(_error::Function,
                                       sense)
    _error("Unexpected sense `$value`. The sense must be an",
           " `MOI.OptimizatonSense`, `Min` or `Max`.")
end
function throw_error_for_invalid_sense(_error::Function,
                                       sense::MOI.OptimizationSense)
    return sense
end

# TODO: Add a docstring.
macro objective(model, args...)
    _error(str...) = macro_error(:objective, (model, args...), str...)

    # We don't overwrite `model` as it is used in `_error`
    esc_model = esc(model)
    if length(args) != 2
        # Either just an objective sense, or just an expression.
        _error("needs three arguments: model, objective sense (Max or Min) and expression.")
    end
    sense, x = args
    sense_expr = moi_sense(_error, sense)
    newaff, parsecode = parseExprToplevel(x, :q)
    code = quote
        q = Val{false}()
        $parsecode
        set_objective($esc_model, $sense_expr, $newaff)
    end
    return assert_validmodel(esc_model, macro_return(code, newaff))
end

# Return a standalone, unnamed expression
# ex = @Expression(2x + 3y)
# Currently for internal use only.
macro Expression(x)
    newaff, parsecode = parseExprToplevel(x, :q)
    code = quote
        q = Val{false}()
        $parsecode
    end
    return macro_return(code, newaff)
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

    refcall, idxvars, idxsets, condition = buildrefsets(c, variable)
    newaff, parsecode = parseExprToplevel(x, :q)
    code = quote
        q = Val{false}()
        $parsecode
    end
    if isa(c,Expr)
        code = quote
            $code
            (isa($newaff,AffExpr) || isa($newaff,Number) || isa($newaff,VariableRef)) || error("Collection of expressions with @expression must be linear. For quadratic expressions, use your own array.")
        end
    end
    code = quote
        $code
        $(refcall) = $newaff
    end
    code = getloopedcode(variable, code, condition, idxvars, idxsets, :AffExpr, requestedcontainer)
    # don't do anything with the model, but check that it's valid anyway
    if anonvar
        macro_code = macro_return(code, variable)
    else
        macro_code = macro_assign_and_return(code, variable, getname(c),
                                             model_for_registering = m)
    end
    return assert_validmodel(m, macro_code)
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

# Returns the type of what `add_variable(::Model, build_variable(...))` would return where `...` represents the positional arguments.
# Example: `@variable m [1:3] foo` will allocate an vector of element type `variable_type(m, foo)`
# Note: it needs to be implemented by all `AbstractModel`s
variable_type(m::Model) = VariableRef
# Returns a new variable. Additional positional arguments can be used to dispatch the call to a different method.
# The return type should only depends on the positional arguments for `variable_type` to make sense. See the @variable macro doc for more details.
# Example: `@variable m x` foo will call `build_variable(_error, info, foo)`
function build_variable(_error::Function, info::VariableInfo; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    return ScalarVariable(info)
end

const EMPTYSTRING = ""

function macro_error(macroname, args, str...)
    error("In @$macroname($(join(args,","))): ", str...)
end

# Given a base_name and idxvars, returns an expression that constructs the name
# of the object. For use within macros only.
function namecall(base_name, idxvars)
    if isempty(idxvars) || base_name == ""
        return base_name
    end
    ex = Expr(:call, :string, base_name, "[")
    for i in 1:length(idxvars)
        # Converting the arguments to strings before concatenating is faster:
        # https://github.com/JuliaLang/julia/issues/29550.
        esc_idxvar = esc(idxvars[i])
        push!(ex.args, :(string($esc_idxvar)))
        i < length(idxvars) && push!(ex.args, ",")
    end
    push!(ex.args, "]")
    return ex
end

reverse_sense(::Val{:<=})   = Val(:>=)
reverse_sense(::Val{:≤})    = Val(:≥)
reverse_sense(::Val{:>=})   = Val(:<=)
reverse_sense(::Val{:≥})    = Val(:≤)
reverse_sense(::Val{:(==)}) = Val(:(==))

"""
    parse_one_operator_variable(_error::Function, infoexpr::VariableInfoExpr, sense::Val{S}, value) where S

Update `infoexr` for a variable expression in the `@variable` macro of the form `variable name S value`.
"""
parse_one_operator_variable(_error::Function, infoexpr::VariableInfoExpr, ::Union{Val{:<=}, Val{:≤}}, upper) = set_upper_bound_or_error(_error, infoexpr, upper)
parse_one_operator_variable(_error::Function, infoexpr::VariableInfoExpr, ::Union{Val{:>=}, Val{:≥}}, lower) = set_lower_bound_or_error(_error, infoexpr, lower)
parse_one_operator_variable(_error::Function, infoexpr::VariableInfoExpr, ::Val{:(==)}, value) = fix_or_error(_error, infoexpr, value)
parse_one_operator_variable(_error::Function, infoexpr::VariableInfoExpr, ::Val{S}, value) where S = _error("Unknown sense $S.")

# There is not way to determine at parsing time which of lhs or rhs is the
# variable name and which is the value if both are symbols. For instance,
# lhs could be the Symbol `:x` and rhs could be the Symbol `:a` where a
# variable `a` is assigned to 1 in the local scope. Knowing this, we know
# that `x` is the variable name but at parse time there is now way to know
# that `a` has a value.
# In that case we assume the variable is the lhs.
function parse_variable(_error::Function, infoexpr::VariableInfoExpr, sense::Symbol, var, value)
    parse_one_operator_variable(_error, infoexpr, Val(sense), esc_nonconstant(value))
    var
end

# If the lhs is a number and not the rhs, we can deduce that the rhs is
# the variable.
function parse_variable(_error::Function, infoexpr::VariableInfoExpr, sense::Symbol, value::Number, var)
    parse_one_operator_variable(_error, infoexpr, reverse_sense(Val(sense)), esc_nonconstant(value))
    var
end

function parse_ternary_variable(_error::Function, infoexpr::VariableInfoExpr,
                              ::Union{Val{:<=}, Val{:≤}}, lower,
                              ::Union{Val{:<=}, Val{:≤}}, upper)
    set_lower_bound_or_error(_error, infoexpr, lower)
    set_upper_bound_or_error(_error, infoexpr, upper)
end
function parse_ternary_variable(_error::Function, infoexpr::VariableInfoExpr,
                              ::Union{Val{:>=}, Val{:≥}}, upper,
                              ::Union{Val{:>=}, Val{:≥}}, lower)
    parse_ternary_variable(_error, infoexpr, Val(:≤), lower, Val(:≤), upper)
end
function parse_ternary_variable(_error::Function, infoexpr::VariableInfoExpr,
                              ::Val, lvalue,
                              ::Val, rvalue)
    _error("Use the form lb <= ... <= ub.")
end
function parse_variable(_error::Function, infoexpr::VariableInfoExpr, lvalue, lsign::Symbol, var, rsign::Symbol, rvalue)
    # lvalue lsign var rsign rvalue
    parse_ternary_variable(_error, infoexpr, Val(lsign), esc_nonconstant(lvalue), Val(rsign), esc_nonconstant(rvalue))
    var
end

"""
    @variable(model, kwargs...)

Add an *anonymous* (see [Names](@ref)) variable to the model `model` described
by the keyword arguments `kwargs` and returns the variable.

    @variable(model, expr, args..., kwargs...)

Add a variable to the model `model` described by the expression `expr`, the
positional arguments `args` and the keyword arguments `kwargs`. The expression
`expr` can either be (note that in the following the symbol `<=` can be used
instead of `≤` and the symbol `>=`can be used instead of `≥`)

* of the form `varexpr` creating variables described by `varexpr`;
* of the form `varexpr ≤ ub` (resp. `varexpr ≥ lb`) creating variables described by
  `varexpr` with upper bounds given by `ub` (resp. lower bounds given by `lb`);
* of the form `varexpr == value` creating variables described by `varexpr` with
  fixed values given by `value`; or
* of the form `lb ≤ varexpr ≤ ub` or `ub ≥ varexpr ≥ lb` creating variables
  described by `varexpr` with lower bounds given by `lb` and upper bounds given
  by `ub`.

The expression `varexpr` can either be

* of the form `varname` creating a scalar real variable of name `varname`;
* of the form `varname[...]` or `[...]` creating a container of variables (see
  [Containers in macro](@ref).

The recognized positional arguments in `args` are the following:

* `Bin`: Sets the variable to be binary, i.e. either 0 or 1.
* `Int`: Sets the variable to be integer, i.e. one of ..., -2, -1, 0, 1, 2, ...
* `Symmetric`: Only available when creating a square matrix of variables, i.e.
  when `varexpr` is of the form `varname[1:n,1:n]` or `varname[i=1:n,j=1:n]`.
  It creates a symmetric matrix of variable, that is, it only creates a
  new variable for `varname[i,j]` with `i ≤ j` and sets `varname[j,i]` to the
  same variable as `varname[i,j]`.
* `PSD`: The square matrix of variable is both `Symmetric` and constrained to be
  positive semidefinite.

The recognized keyword arguments in `kwargs` are the following:

* `base_name`: Sets the name prefix used to generate variable names. It
  corresponds to the variable name for scalar variable, otherwise, the
  variable names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
* `lower_bound`: Sets the value of the variable lower bound.
* `upper_bound`: Sets the value of the variable upper bound.
* `start`: Sets the variable starting value used as initial guess in optimization.
* `binary`: Sets whether the variable is binary or not.
* `integer`: Sets whether the variable is integer or not.
* `variable_type`: See the "Note for extending the variable macro" section below.
* `container`: Specify the container type, see [Containers in macro](@ref).

## Examples

The following are equivalent ways of creating a variable `x` of name `x` with
lower bound 0:
```julia
# Specify everything in `expr`
@variable(model, x >= 0)
# Specify the lower bound using a keyword argument
@variable(model, x, lower_bound=0)
# Specify everything in `kwargs`
x = @variable(model, base_name="x", lower_bound=0)
```

The following are equivalent ways of creating a `JuMPArray` of index set
`[:a, :b]` and with respective upper bounds 2 and 3 and names `x[a]` and `x[b].
```julia
ub = Dict(:a => 2, :b => 3)
# Specify everything in `expr`
@variable(model, x[i=keys(ub)] <= ub[i])
# Specify the upper bound using a keyword argument
@variable(model, x[i=keys(ub)], upper_bound=ub[i])
```

## Note for extending the variable macro

The single scalar variable or each scalar variable of the container are created
using `add_variable(model, build_variable(_error, info, extra_args...;
extra_kwargs...))` where

* `model` is the model passed to the `@variable` macro;
* `_error` is an error function with a single `String` argument showing the
  `@variable` call in addition to the error message given as argument;
* `info` is the `VariableInfo` struct containing the information gathered in
  `expr`, the recognized keyword arguments (except `base_name` and
  `variable_type`) and the recognized positional arguments (except `Symmetric`
  and `PSD`);
* `extra_args` are the unrecognized positional arguments of `args` plus the
  value of the `variable_type` keyword argument if present. The `variable_type`
  keyword argument allows the user to pass a position argument to
  `build_variable` without the need to give a positional argument to
  `@variable`. In particular, this allows the user to give a positional
  argument to the `build_variable` call when using the anonymous single variable
  syntax `@variable(model, kwargs...)`; and
* `extra_kwargs` are the unrecognized keyword argument of `kwargs`.

## Examples

The following creates a variable `x` of name `x` with lower_bound 0 as with the first
example above but does it without using the `@variable` macro
```julia
info = VariableInfo(true, 0, false, NaN, false, NaN, false, NaN, false, false)
JuMP.add_variable(model, JuMP.build_variable(error, info), "x")
```

The following creates a `JuMPArray` of index set `[:a, :b]` and with respective
upper bounds 2 and 3 and names `x[a]` and `x[b]` as with the second example
above but does it without using the `@variable` macro
```julia
# Without the `@variable` macro
data = Vector{JuMP.variable_type(model)}(undef, length(keys(ub)))
x = JuMPArray(data, keys(ub))
for i in keys(ub)
    info = VariableInfo(false, NaN, true, ub[i], false, NaN, false, NaN, false, false)
    x[i] = JuMP.add_variable(model, JuMP.build_variable(error, info), "x[\$i]")
end
```

The following are equivalent ways of creating a `Matrix` of size
`N x N` with variables custom variables created with a JuMP extension using
the `Poly(X)` positional argument to specify its variables:
```julia
# Using the `@variable` macro
@variable(model, x[1:N,1:N], Symmetric, Poly(X))
# Without the `@variable` macro
x = Matrix{JuMP.variable_type(model, Poly(X))}(N, N)
info = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)
for i in 1:N, j in i:N
    x[i,j] = x[j,i] = JuMP.add_variable(model, build_variable(error, info, Poly(X)), "x[\$i,\$j]")
end
```
"""
macro variable(args...)
    _error(str...) = macro_error(:variable, args, str...)

    model = esc(args[1])

    extra, kwargs, requestedcontainer = extract_kwargs(args[2:end])

    # if there is only a single non-keyword argument, this is an anonymous
    # variable spec and the one non-kwarg is the model
    if length(extra) == 0
        x = gensym()
        anon_singleton = true
    else
        x = popfirst!(extra)
        if x in [:Int,:Bin,:PSD]
            _error("Ambiguous variable name $x detected. Use the \"category\" keyword argument to specify a category for an anonymous variable.")
        end
        anon_singleton = false
    end

    info_kwargs = filter(is_info_keyword, kwargs)
    extra_kwargs = filter(kw -> kw.args[1] != :base_name && kw.args[1] != :variable_type && !is_info_keyword(kw), kwargs)
    base_name_kwargs = filter(kw -> kw.args[1] == :base_name, kwargs)
    variable_type_kwargs = filter(kw -> kw.args[1] == :variable_type, kwargs)
    infoexpr = VariableInfoExpr(; keywordify.(info_kwargs)...)

    # There are four cases to consider:
    # x                                       | type of x | x.head
    # ----------------------------------------+-----------+------------
    # var                                     | Symbol    | NA
    # var[1:2]                                | Expr      | :ref
    # var <= ub or var[1:2] <= ub             | Expr      | :call
    # lb <= var <= ub or lb <= var[1:2] <= ub | Expr      | :comparison
    # In the two last cases, we call parse_variable
    explicit_comparison = isexpr(x, :comparison) || isexpr(x, :call)
    if explicit_comparison
        var = parse_variable(_error, infoexpr, x.args...)
    else
        var = x
    end

    anonvar = isexpr(var, :vect) || isexpr(var, :vcat) || anon_singleton
    anonvar && explicit_comparison && error("Cannot use explicit bounds via >=, <= with an anonymous variable")
    variable = gensym()
    # TODO: Should we generate non-empty default names for variables?
    name = getname(var)
    if isempty(base_name_kwargs)
        base_name = anonvar ? "" : string(name)
    else
        base_name = esc(base_name_kwargs[1].args[2])
    end

    if !isa(name, Symbol) && !anonvar
        Base.error("Expression $name should not be used as a variable name. Use the \"anonymous\" syntax $name = @variable(model, ...) instead.")
    end

    # process keyword arguments
    obj = nothing

    sdp = any(t -> (t == :PSD), extra)
    symmetric = (sdp || any(t -> (t == :Symmetric), extra))
    extra = filter(x -> (x != :PSD && x != :Symmetric), extra) # filter out PSD and sym tag
    for ex in extra
        if ex == :Int
            set_integer_or_error(_error, infoexpr)
        elseif ex == :Bin
            set_binary_or_error(_error, infoexpr)
        end
    end
    extra = esc.(filter(ex -> !(ex in [:Int,:Bin]), extra))
    if !isempty(variable_type_kwargs)
        push!(extra, esc(variable_type_kwargs[1].args[2]))
    end

    info = constructor_expr(infoexpr)
    if isa(var,Symbol)
        # Easy case - a single variable
        sdp && _error("Cannot add a semidefinite scalar variable")
        buildcall = :( build_variable($_error, $info, $(extra...)) )
        addkwargs!(buildcall, extra_kwargs)
        variablecall = :( add_variable($model, $buildcall, $base_name) )
        # The looped code is trivial here since there is a single variable
        creationcode = :($variable = $variablecall)
        final_variable = variable
    else
        isa(var,Expr) || _error("Expected $var to be a variable name")

        # We now build the code to generate the variables (and possibly the JuMPDict
        # to contain them)
        refcall, idxvars, idxsets, condition = buildrefsets(var, variable)
        clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? () : idxsets[i])

        # Code to be used to create each variable of the container.
        buildcall = :( build_variable($_error, $info, $(extra...)) )
        addkwargs!(buildcall, extra_kwargs)
        variablecall = :( add_variable($model, $buildcall, $(namecall(base_name, idxvars))) )
        code = :( $(refcall) = $variablecall )
        # Determine the return type of add_variable. This is needed to create the container holding them.
        vartype = :( variable_type($model, $(extra...)) )

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
                if VERSION >= v"0.7-"
                    (isexpr(rng,:call) && length(rng.args) == 3 && rng.args[1] == :(:) && rng.args[2] == 1) ||
                        _error("Index sets for SDP variables must be ranges of the form 1:N")
                else
                    (isexpr(rng,:(:)) && rng.args[1] == 1 && length(rng.args) == 2) ||
                        _error("Index sets for SDP variables must be ranges of the form 1:N")
                end
            end

            if infoexpr.has_lb || infoexpr.has_ub
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
            creationcode = quote
                $dimension_check || error("Cannot construct symmetric variables with nonsquare dimensions.")
                $(getloopedcode(variable, code, condition, idxvars, idxsets, vartype, requestedcontainer; lowertri=symmetric))
                $(if sdp
                    quote
                        JuMP.add_constraint($model, JuMP.build_constraint($_error, Symmetric($variable), JuMP.PSDCone()))
                    end
                end)
            end
            final_variable = :(Symmetric($variable))
        else
            creationcode = getloopedcode(variable, code, condition, idxvars, idxsets, vartype, requestedcontainer)
            final_variable = variable
        end
    end
    if anonvar
        # Anonymous variable, no need to register it in the model-level
        # dictionary nor to assign it to a variable in the user scope.
        # We simply return the variable
        macro_code = macro_return(creationcode, final_variable)
    else
        # We register the variable reference to its name and
        # we assign it to a variable in the local scope of this name
        macro_code = macro_assign_and_return(creationcode, variable, name,
                                             final_variable=final_variable,
                                             model_for_registering = model)
    end
    return assert_validmodel(model, macro_code)
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

# TODO: Add a docstring.
macro NLobjective(model, sense, x)
    _error(str...) = macro_error(:NLobjective, (model, sense, x), str...)
    sense_expr = moi_sense(_error, sense)
    ex = gensym()
    code = quote
        $ex = $(processNLExpr(model, x))
        set_objective($(esc(model)), $sense_expr, $ex)
    end
    return assert_validmodel(esc(model), macro_return(code, ex))
end

# TODO: Add a docstring.
macro NLconstraint(m, x, extra...)
    esc_m = esc(m)
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

    # Strategy: build up the code for non-macro add_constraint, and if needed
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
            c = NonlinearConstraint($(processNLExpr(m, lhs)), $lb, $ub)
            push!($esc_m.nlp_data.nlconstr, c)
            $(refcall) = ConstraintRef($esc_m, NonlinearConstraintIndex(length($esc_m.nlp_data.nlconstr)), ScalarShape())
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
            push!($esc_m.nlp_data.nlconstr, c)
            $(refcall) = ConstraintRef($esc_m, NonlinearConstraintIndex(length($esc_m.nlp_data.nlconstr)), ScalarShape())
        end
    else
        # Unknown
        error("in @NLconstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2")
    end
    looped = getloopedcode(variable, code, condition, idxvars, idxsets, :(ConstraintRef{Model,NonlinearConstraintIndex}), requestedcontainer)
    creation_code = quote
        initNLP($esc_m)
        $looped
    end
    if anonvar
        macro_code = macro_return(creation_code, variable)
    else
        macro_code = macro_assign_and_return(creation_code, variable,
                                             getname(c),
                                             model_for_registering = esc_m)
    end
    return assert_validmodel(esc_m, macro_code)
end

# TODO: Add a docstring.
macro NLexpression(args...)
    args, kwargs, requestedcontainer = extract_kwargs(args)
    if length(args) <= 1
        error("in @NLexpression: To few arguments ($(length(args))); must pass the model and nonlinear expression as arguments.")
    elseif length(args) == 2
        m, x = args
        c = gensym()
    elseif length(args) == 3
        m, c, x = args
    end
    if length(args) > 3 || length(kwargs) > 0
        error("in @NLexpression: To many arguments ($(length(args))).")
    end

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    variable = gensym()

    refcall, idxvars, idxsets, condition = buildrefsets(c, variable)
    code = quote
        $(refcall) = NonlinearExpression($(esc(m)), $(processNLExpr(m, x)))
    end
    creation_code = getloopedcode(variable, code, condition, idxvars, idxsets, :NonlinearExpression, requestedcontainer)
    if anonvar
        macro_code = macro_return(creation_code, variable)
    else
        macro_code = macro_assign_and_return(creation_code, variable,
                                             getname(c),
                                             model_for_registering = esc(m))
    end
    return assert_validmodel(esc(m), macro_code)
end

"""
    @NLparameter(model, param == value)

Create and return a nonlinear parameter `param` attached to the model `model`
with initial value set to `value`. Nonlinear parameters may be used only in
nonlinear expressions.

# Example
```jldoctest
model = Model()
@NLparameter(model, x == 10)
JuMP.value(x)

# output
10.0
```

    @NLparameter(model, param_collection[...] == value_expr)

Create and return a collection of nonlinear parameters `param_collection`
attached to the model `model` with initial value set to `value_expr` (may
depend on index sets).
Uses the same syntax for specifying index sets as [`@variable`](@ref).

# Example
```jldoctest
model = Model()
@NLparameter(model, y[i = 1:10] == 2 * i)
JuMP.value(y[9])

# output
18.0
```
"""
macro NLparameter(m, ex, extra...)

    extra, kwargs, requestedcontainer = extract_kwargs(extra)
    (length(extra) == 0 && length(kwargs) == 0) || error("in @NLperameter: too many arguments.")
    if !isexpr(ex, :call) || length(ex.args) != 3 || ex.args[1] != :(==)
        error("In @NLparameter($m, $ex): syntax error.")
    end
    c = ex.args[2]
    x = ex.args[3]

    anonvar = isexpr(c, :vect) || isexpr(c, :vcat)
    if anonvar
        error("In @NLparameter($m, $ex): Anonymous nonlinear parameter syntax is not currently supported")
    end
    m = esc(m)
    variable = gensym()

    refcall, idxvars, idxsets, condition = buildrefsets(c, variable)
    code = quote
        if !isa($(esc(x)), Number)
            error(string("in @NLparameter (", $(string(ex)), "): expected ",
                         $(string(x))," to be a number."))
        end
        $(refcall) = newparameter($m, $(esc(x)))
    end
    creation_code = getloopedcode(variable, code, condition, idxvars, idxsets, :NonlinearParameter, :Auto)

    # TODO: NLparameters are not registered in the model because we don't yet
    # have an anonymous version.
    macro_code = macro_assign_and_return(creation_code, variable,
                                         getname(c))
    return assert_validmodel(m, macro_code)
end
