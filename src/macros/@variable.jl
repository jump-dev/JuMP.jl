#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    @variable(model, expr, args..., kw_args...)

Add a variable to the model `model` described by the expression `expr`, the
positional arguments `args` and the keyword arguments `kw_args`.

## Anonymous and named variables

`expr` must be one of the forms:

 * Omitted, like `@variable(model)`, which creates an anonymous variable
 * A single symbol like `@variable(model, x)`
 * A container expression like `@variable(model, x[i=1:3])`
 * An anoymous container expression like `@variable(model, [i=1:3])`

## Bounds

In addition, the expression can have bounds, such as:

 * `@variable(model, x >= 0)`
 * `@variable(model, x <= 0)`
 * `@variable(model, x == 0)`
 * `@variable(model, 0 <= x <= 1)`

and bounds can depend on the indices of the container expressions:

 * `@variable(model, -i <= x[i=1:3] <= i)`

## Sets

You can explicitly specify the set to which the variable belongs:

 * `@variable(model, x in MOI.Interval(0.0, 1.0))`

 For more information on this syntax, read
[Variables constrained on creation](@ref).

## Positional arguments

The recognized positional arguments in `args` are the following:

 * `Bin`: restricts the variable to the [`MOI.ZeroOne`](@ref) set, that is,
   `{0, 1}`. For example, `@variable(model, x, Bin)`. Note: you cannot use
   `@variable(model, Bin)`, use the `binary` keyword instead.
 * `Int`: restricts the variable to the set of integers, that is, ..., -2, -1,
    0, 1, 2, ... For example, `@variable(model, x, Int)`. Note: you cannot use
    `@variable(model, Int)`, use the `integer` keyword instead.
 * `Symmetric`: Only available when creating a square matrix of variables, i.e.,
   when `expr` is of the form `varname[1:n,1:n]` or `varname[i=1:n,j=1:n]`,
   it creates a symmetric matrix of variables.
 * `PSD`: A restrictive extension to `Symmetric` which constraints a square
   matrix of variables to `Symmetric` and constrains to be positive
   semidefinite.

## Keyword arguments

Four keyword arguments are useful in all cases:

 * `base_name`: Sets the name prefix used to generate variable names. It
   corresponds to the variable name for scalar variable, otherwise, the
   variable names are set to `base_name[...]` for each index `...` of the axes
   `axes`.
 * `start::Float64`: specify the value passed to `set_start_value` for each
   variable
 * `container`: specify the container type. See
   [Forcing the container type](@ref variable_forcing) for more information.
 * `set_string_name::Bool = true`: control whether to set the
   [`MOI.VariableName`](@ref) attribute. Passing `set_string_name = false` can
   improve performance.

Other keyword arguments are needed to disambiguate sitations with anonymous
variables:

 * `lower_bound::Float64`: an alternative to `x >= lb`, sets the value of the
   variable lower bound.
 * `upper_bound::Float64`: an alternative to `x <= ub`, sets the value of the
   variable upper bound.
 * `binary::Bool`: an alternative to passing `Bin`, sets whether the variable
   is binary or not.
 * `integer::Bool`: an alternative to passing `Int`, sets whether the variable
   is integer or not.
 * `set::MOI.AbstractSet`: an alternative to using `x in set`
 * `variable_type`: used by JuMP extensions. See
   [Extend `@variable`](@ref extend_variable_macro) for more information.

## Example

The following are equivalent ways of creating a variable `x` of name `x` with
lower bound 0:
```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0)
x
```

```jldoctest
julia> model = Model();

julia> @variable(model, x, lower_bound = 0)
x
```

```jldoctest
julia> model = Model();

julia> x = @variable(model, base_name = "x", lower_bound = 0)
x
```

Other examples:

```jldoctest
julia> model = Model();

julia> @variable(model, x[i=1:3] <= i, Int, start = sqrt(i), lower_bound = -i)
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> @variable(model, y[i=1:3], container = DenseAxisArray, set = MOI.ZeroOne())
1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
    Dimension 1, Base.OneTo(3)
And data, a 3-element Vector{VariableRef}:
 y[1]
 y[2]
 y[3]

julia> @variable(model, z[i=1:3], set_string_name = false)
3-element Vector{VariableRef}:
 _[7]
 _[8]
 _[9]
```
"""
macro variable(args...)
    error_fn(str...) = _macro_error(:variable, args, __source__, str...)
    # We need to re-order the parameters here to account for cases like
    # `@variable(model; integer = true)`, since Julia handles kwargs by placing
    # them first(!) in the list of arguments.
    args = _reorder_parameters(args)
    model = esc(args[1])
    if length(args) >= 2 && Meta.isexpr(args[2], :block)
        error_fn("Invalid syntax. Did you mean to use `@variables`?")
    end
    extra, kw_args, requested_container =
        Containers._extract_kw_args(args[2:end])

    # if there is only a single non-keyword argument, this is an anonymous
    # variable spec and the one non-kwarg is the model
    if length(extra) == 0
        x = nothing
        anon_singleton = true
    else
        x = popfirst!(extra)
        if x == :Int
            error_fn(
                "Ambiguous variable name $x detected. To specify an anonymous integer " *
                "variable, use `@variable(model, integer = true)` instead.",
            )
        elseif x == :Bin
            error_fn(
                "Ambiguous variable name $x detected. To specify an anonymous binary " *
                "variable, use `@variable(model, binary = true)` instead.",
            )
        elseif x == :PSD
            error_fn(
                "Size of anonymous square matrix of positive semidefinite anonymous variables is not specified. To specify size of square matrix " *
                "use `@variable(model, [1:n, 1:n], PSD)` instead.",
            )
        end
        anon_singleton = false
    end

    info_kw_args = filter(_is_info_keyword, kw_args)
    extra_kw_args = filter(
        kw -> begin
            kw.args[1] != :base_name &&
                kw.args[1] != :variable_type &&
                kw.args[1] != :set &&
                kw.args[1] != :set_string_name &&
                !_is_info_keyword(kw)
        end,
        kw_args,
    )
    base_name_kw_args = filter(kw -> kw.args[1] == :base_name, kw_args)
    variable_type_kw_args = filter(kw -> kw.args[1] == :variable_type, kw_args)
    set_string_name_kw_args =
        filter(kw -> kw.args[1] == :set_string_name, kw_args)
    info_expr = _VariableInfoExpr(; _keywordify.(info_kw_args)...)

    # There are four cases to consider:
    # x                                       | type of x | x.head
    # ----------------------------------------+-----------+------------
    # var                                     | Symbol    | NA
    # var[1:2]                                | Expr      | :ref
    # var <= ub or var[1:2] <= ub             | Expr      | :call
    # var in set or var[1:2] in set           | Expr      | :call
    # lb <= var <= ub or lb <= var[1:2] <= ub | Expr      | :comparison
    # In the three last cases, we call parse_variable
    explicit_comparison = Meta.isexpr(x, :comparison) || Meta.isexpr(x, :call)
    if explicit_comparison
        var, set = parse_variable(error_fn, info_expr, x.args...)
    else
        var = x
        set = nothing
    end
    anonvar =
        Meta.isexpr(var, :vect) || Meta.isexpr(var, :vcat) || anon_singleton
    if anonvar && explicit_comparison && set === nothing
        error_fn(
            "Cannot use explicit bounds via >=, <= with an anonymous variable",
        )
    end
    variable = gensym()
    # TODO: Should we generate non-empty default names for variables?
    name = something(Containers._get_name(var), "")
    if isempty(base_name_kw_args)
        base_name = anonvar ? "" : string(name)
    else
        base_name = esc(base_name_kw_args[1].args[2])
    end
    set_kw_args = filter(kw -> kw.args[1] == :set, kw_args)
    if length(set_kw_args) == 1
        if set !== nothing
            error_fn(
                "Cannot use set keyword because the variable is already " *
                "constrained to `$set`.",
            )
        end
        set = esc(set_kw_args[1].args[2])
    elseif length(set_kw_args) > 1
        error_fn(
            "`set` keyword argument was given $(length(set_kw_args)) times.",
        )
    end
    for (sym, cone) in (
        :PSD => PSDCone(),
        :Symmetric => SymmetricMatrixSpace(),
        :Hermitian => HermitianMatrixSpace(),
    )
        if any(isequal(sym), extra)
            if set !== nothing
                error_fn(
                    "Cannot pass `$sym` as a positional argument because the " *
                    "variable is already constrained to `$set`.",
                )
            end
            set = cone
            filter!(!isequal(sym), extra)
        end
    end
    for ex in extra
        if ex == :Int
            _set_integer_or_error(error_fn, info_expr)
        elseif ex == :Bin
            _set_binary_or_error(error_fn, info_expr)
        end
    end
    extra = esc.(filter(ex -> !(ex in [:Int, :Bin]), extra))
    if !isempty(variable_type_kw_args)
        push!(extra, esc(variable_type_kw_args[1].args[2]))
    end

    info = _constructor_expr(info_expr)
    if isa(var, Symbol) || var === nothing
        # Easy case - a single variable
        name_code = base_name
    else
        isa(var, Expr) || error_fn("Expected $var to be a variable name")
        # We now build the code to generate the variables (and possibly the
        # SparseAxisArray to contain them)
        idxvars, indices = Containers.build_ref_sets(error_fn, var)
        if args[1] in idxvars
            error_fn(
                "Index $(args[1]) is the same symbol as the model. Use a " *
                "different name for the index.",
            )
        end
        name_code = _name_call(base_name, idxvars)
        if set !== nothing
            name_code = Containers.container_code(
                idxvars,
                indices,
                name_code,
                requested_container,
            )
        end
    end

    # Code to be used to create each variable of the container.
    buildcall = :(build_variable($error_fn, $info, $(extra...)))
    _add_kw_args(buildcall, extra_kw_args)
    if set !== nothing
        if isa(var, Symbol) || var === nothing
            scalar_variables = buildcall
        else
            scalar_variables = Containers.container_code(
                idxvars,
                indices,
                buildcall,
                requested_container,
            )
            if any(Base.Fix2(Containers.depends_on, set), idxvars)
                set = Containers.container_code(
                    idxvars,
                    indices,
                    set,
                    requested_container,
                )
            end
        end
        buildcall = :(build_variable($error_fn, $scalar_variables, $set))
    end
    buildcall = :(model_convert($model, $buildcall))
    new_name_code = if isempty(set_string_name_kw_args)
        Expr(:if, :(set_string_names_on_creation($model)), name_code, "")
    else
        Expr(:if, esc(set_string_name_kw_args[1].args[2]), name_code, "")
    end
    variablecall = :(add_variable($model, $buildcall, $new_name_code))
    if isa(var, Symbol) || var === nothing || set !== nothing
        # The looped code is trivial here since there is a single variable
        creation_code = variablecall
    else
        creation_code = Containers.container_code(
            idxvars,
            indices,
            variablecall,
            requested_container,
        )
    end
    return _finalize_macro(
        model,
        creation_code,
        __source__;
        register_name = Containers._get_name(var),
        wrap_let = true,
    )
end

"""
    @variables(model, args...)

Adds multiple variables to model at once, in the same fashion as the
[`@variable`](@ref) macro.

The model must be the first argument, and multiple variables can be added on
multiple lines wrapped in a `begin ... end` block.

The macro returns a tuple containing the variables that were defined.

## Example

```jldoctest
julia> model = Model();

julia> @variables(model, begin
           x
           y[i = 1:2] >= 0, (start = i)
           z, Bin, (start = 0, base_name = "Z")
       end)
(x, VariableRef[y[1], y[2]], Z)
```

!!! note
    Keyword arguments must be contained within parentheses (refer to the example
    above).
"""
macro variables(model, block)
    return _plural_macro_code(model, block, Symbol("@variable"))
end

"""
    parse_variable(error_fn::Function, ::_VariableInfoExpr, args...)

A hook for extensions to intercept the parsing of inequality constraints in the
[`@variable`](@ref) macro.
"""
function parse_variable(error_fn::Function, ::_VariableInfoExpr, args...)
    return error_fn(
        "Invalid syntax: your syntax is wrong, but we don't know why. " *
        "Consult the documentation for various ways to create variables in " *
        "JuMP.",
    )
end

# There is not way to determine at parsing time which of lhs or rhs is the
# variable name and which is the value if both are symbols. For instance,
# lhs could be the Symbol `:x` and rhs could be the Symbol `:a` where a
# variable `a` is assigned to 1 in the local scope. Knowing this, we know
# that `x` is the variable name but at parse time there is now way to know
# that `a` has a value.
# In that case we assume the variable is the lhs.
function parse_variable(
    error_fn::Function,
    info_expr::_VariableInfoExpr,
    sense::Symbol,
    var,
    value,
)
    set = parse_one_operator_variable(
        error_fn,
        info_expr,
        Val(sense),
        _esc_non_constant(value),
    )
    return var, set
end

# If the lhs is a number and not the rhs, we can deduce that the rhs is
# the variable.
function parse_variable(
    error_fn::Function,
    info_expr::_VariableInfoExpr,
    sense::Symbol,
    value::Number,
    var,
)
    set = parse_one_operator_variable(
        error_fn,
        info_expr,
        reverse_sense(Val(sense)),
        _esc_non_constant(value),
    )
    return var, set
end

function parse_variable(
    error_fn::Function,
    info_expr::_VariableInfoExpr,
    lvalue,
    lsign::Symbol,
    var,
    rsign::Symbol,
    rvalue,
)
    set = parse_ternary_variable(
        error_fn,
        info_expr,
        Val(lsign),
        _esc_non_constant(lvalue),
        Val(rsign),
        _esc_non_constant(rvalue),
    )
    return var, set
end

"""
    parse_one_operator_variable(
        error_fn::Function,
        info_expr::_VariableInfoExpr,
        sense::Val{S},
        value,
    ) where {S}

Update `infoexr` for a variable expression in the `@variable` macro of the form
`variable name S value`.
"""
function parse_one_operator_variable(
    error_fn::Function,
    ::_VariableInfoExpr,
    ::Val{S},
    ::Any,
) where {S}
    return error_fn("unsupported operator $S")
end

function parse_one_operator_variable(
    ::Function,
    ::_VariableInfoExpr,
    ::Union{Val{:in},Val{:∈}},
    set,
)
    return set
end

function parse_one_operator_variable(
    error_fn::Function,
    info_expr::_VariableInfoExpr,
    ::Union{Val{:<=},Val{:≤},Val{:.<=},Val{:.≤}},
    upper,
)
    _set_upper_bound_or_error(error_fn, info_expr, upper)
    return
end

function parse_one_operator_variable(
    error_fn::Function,
    info_expr::_VariableInfoExpr,
    ::Union{Val{:>=},Val{:≥},Val{:.>=},Val{:.≥}},
    lower,
)
    _set_lower_bound_or_error(error_fn, info_expr, lower)
    return
end

function parse_one_operator_variable(
    error_fn::Function,
    info_expr::_VariableInfoExpr,
    ::Union{Val{:(==)},Val{:.==}},
    value,
)
    _fix_or_error(error_fn, info_expr, value)
    return
end

function parse_one_operator_variable(
    error_fn::Function,
    ::_VariableInfoExpr,
    ::Val{:>},
    ::Any,
)
    return error_fn(
        "unsupported operator `>`.\n\n" *
        "JuMP does not support strict inequalities, use `>=` instead.\n\n" *
        "If you require a strict inequality, you will need to use a " *
        "tolerance. For example, instead of `x > 1`, do `x >= 1 + 1e-4`. " *
        "If the variable must take integer values, use a tolerance of " *
        "`1.0`. If the variable may take continuous values, note that this " *
        "work-around can cause numerical issues, and your bound may not " *
        "hold exactly.",
    )
end

function parse_one_operator_variable(
    error_fn::Function,
    ::_VariableInfoExpr,
    ::Val{:<},
    ::Any,
)
    return error_fn(
        "unsupported operator `<`.\n\n" *
        "JuMP does not support strict inequalities, use `<=` instead.\n\n" *
        "If you require a strict inequality, you will need to use a " *
        "tolerance. For example, instead of `x < 1`, do `x <= 1 - 1e-4`. " *
        "If the variable must take integer values, use a tolerance of " *
        "`1.0`. If the variable may take continuous values, note that this " *
        "work-around can cause numerical issues, and your bound may not " *
        "hold exactly.",
    )
end

"""
    parse_ternary_variable(error_fn, info_expr, lhs_sense, lhs, rhs_sense, rhs)

A hook for JuMP extensiosn to intercept the parsing of a `:comparison`
expression, which has the form `lhs lhs_sense variable rhs_sense rhs`.
"""
function parse_ternary_variable(
    error_fn::Function,
    info_expr::_VariableInfoExpr,
    ::Val{A},
    lb,
    ::Val{B},
    ub,
) where {A,B}
    return error_fn(
        "unsupported mix of comparison operators `$lb $A ... $B $ub`.\n\n" *
        "Two-sided variable bounds must of the form `$lb <= ... <= $ub` or " *
        "`$ub >= ... >= $lb`.",
    )
end

function parse_ternary_variable(
    error_fn::Function,
    info_expr::_VariableInfoExpr,
    ::Union{Val{:<=},Val{:≤},Val{:.<=},Val{:.≤}},
    lower,
    ::Union{Val{:<=},Val{:≤},Val{:.<=},Val{:.≤}},
    upper,
)
    _set_lower_bound_or_error(error_fn, info_expr, lower)
    _set_upper_bound_or_error(error_fn, info_expr, upper)
    return
end

function parse_ternary_variable(
    error_fn::Function,
    info_expr::_VariableInfoExpr,
    ::Union{Val{:>=},Val{:≥},Val{:.>=},Val{:.≥}},
    upper,
    ::Union{Val{:>=},Val{:≥},Val{:.>=},Val{:.≥}},
    lower,
)
    return parse_ternary_variable(
        error_fn,
        info_expr,
        Val(:≤),
        lower,
        Val(:≤),
        upper,
    )
end

"""
    build_variable(
        error_fn::Function,
        info::VariableInfo,
        args...;
        kwargs...,
    )

Return a new [`AbstractVariable`](@ref) object.

This method should only be implemented by developers creating JuMP extensions.
It should never be called by users of JuMP.

## Arguments

 * `error_fn`: a function to call instead of `error`. `error_fn` annotates the
   error message with additional information for the user.
 * `info`: an instance of [`VariableInfo`](@ref). This has a variety of fields
   relating to the variable such as `info.lower_bound` and `info.binary`.
 * `args`: optional additional positional arguments for extending the
   [`@variable`](@ref) macro.
 * `kwargs`: optional keyword arguments for extending the [`@variable`](@ref)
   macro.

See also: [`@variable`](@ref)

!!! warning
    Extensions should define a method with ONE positional argument to dispatch
    the call to a different method. Creating an extension that relies on
    multiple positional arguments leads to `MethodError`s if the user passes the
    arguments in the wrong order.

## Example

```julia
@variable(model, x, Foo)
```
will call
```julia
build_variable(error_fn::Function, info::VariableInfo, ::Type{Foo})
```

Passing special-case positional arguments such as `Bin`, `Int`, and `PSD` is
okay, along with keyword arguments:
```julia
@variable(model, x, Int, Foo(), mykwarg = true)
# or
@variable(model, x, Foo(), Int, mykwarg = true)
```
will call
```julia
build_variable(error_fn::Function, info::VariableInfo, ::Foo; mykwarg)
```
and `info.integer` will be true.

Note that the order of the positional arguments does not matter.
"""
function build_variable(
    error_fn::Function,
    info::VariableInfo,
    args...;
    kwargs...,
)
    if length(args) > 0
        error_fn(
            "Unrecognized positional arguments: $(args). (You may have " *
            "passed it as a positional argument, or as a keyword value to " *
            "`variable_type`.)\n\nIf you're trying to create a JuMP " *
            "extension, you need to implement `build_variable`. Read the " *
            "docstring for more details.",
        )
    end
    for (key, _) in kwargs
        if key == :Bool
            error_fn(
                "Unsupported keyword argument: $key.\n\nIf you intended to " *
                "create a `{0, 1}` decision variable, use the `binary` keyword " *
                "argument instead: `@variable(model, x, binary = true)`.",
            )
        end
        error_fn(
            "Unrecognized keyword argument: $key.\n\nIf you're trying " *
            "to create a JuMP extension, you need to implement " *
            "`build_variable`. Read the docstring for more details.",
        )
    end
    if info.lower_bound isa AbstractArray
        error_fn(
            """
            Passing arrays as variable bounds without indexing them is not supported.

            Instead of:
            ```julia
            @variable(model, x[1:2] >= lb)
            ```
            use
            ```julia
            @variable(model, x[i=1:2] >= lb[i])
            ```
            or
            ```julia
            @variable(model, x[1:2])
            set_lower_bound.(x, lb)
            ```
            """,
        )
    elseif info.upper_bound isa AbstractArray
        error_fn(
            """
            Passing arrays as variable bounds without indexing them is not supported.

            Instead of:
            ```julia
            @variable(model, x[1:2] <= ub)
            ```
            use
            ```julia
            @variable(model, x[i=1:2] <= ub[i])
            ```
            or
            ```julia
            @variable(model, x[1:2])
            set_upper_bound.(x, ub)
            ```
            """,
        )
    elseif info.fixed_value isa AbstractArray
        error_fn(
            """
            Passing arrays as variable bounds without indexing them is not supported.

            Instead of:
            ```julia
            @variable(model, x[1:2] == fx)
            ```
            use
            ```julia
            @variable(model, x[i=1:2] == fx[i])
            ```
            or
            ```julia
            @variable(model, x[1:2])
            fix.(x, fx)
            ```
            """,
        )
    elseif info.start isa AbstractArray
        error_fn(
            """
            Passing arrays as variable starts without indexing them is not supported.

            Instead of:
            ```julia
            @variable(model, x[1:2], start = x0)
            ```
            use
            ```julia
            @variable(model, x[i=1:2], start = x0[i])
            ```
            or
            ```julia
            @variable(model, x[1:2])
            set_start_value.(x, x0)
            ```
            """,
        )
    end
    return ScalarVariable(info)
end

function build_variable(
    error_fn::Function,
    info::VariableInfo,
    ::Type{Bool};
    kwargs...,
)
    return error_fn(
        "Unsupported positional argument `Bool`. If you intended to create a " *
        "`{0, 1}` decision variable, use `Bin` instead. For example, " *
        "`@variable(model, x, Bin)` or `@variable(model, x, binary = true)`.",
    )
end

function build_variable(
    ::Function,
    variable::AbstractVariable,
    set::MOI.AbstractScalarSet,
)
    return VariableConstrainedOnCreation(variable, set)
end

function build_variable(
    error_fn::Function,
    variables::AbstractArray{<:ScalarVariable},
    sets::AbstractArray{<:MOI.AbstractScalarSet},
)
    if length(variables) != length(sets)
        return error_fn(
            "Dimensions must match. Got a vector of scalar variables with" *
            "$(length(variables)) elements and a vector of " *
            "scalar sets with $(length(sets)).",
        )
    end
    return VariableConstrainedOnCreation.(variables, sets)
end

function build_variable(
    ::Function,
    variables::AbstractArray{<:ScalarVariable},
    set::MOI.AbstractScalarSet,
)
    return VariableConstrainedOnCreation.(variables, Ref(set))
end

function build_variable(
    error_fn::Function,
    ::ScalarVariable,
    sets::AbstractArray{<:MOI.AbstractScalarSet},
)
    return error_fn(
        "It is not possible to add a scalar variable in an Array of " *
        "sets. Either add an Array of scalar variables in a scalar set or " *
        "add an Array of scalar variables in an Array of scalar sets of " *
        "the same dimension.",
    )
end

function build_variable(
    ::Function,
    variables::Vector{<:AbstractVariable},
    set::MOI.AbstractVectorSet,
)
    return VariablesConstrainedOnCreation(variables, set)
end
