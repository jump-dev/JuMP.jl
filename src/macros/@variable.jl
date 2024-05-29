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
 * An anonymous container expression like `@variable(model, [i=1:3])`

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
 * `Symmetric`: Only available when creating a square matrix of variables, that
   is when `expr` is of the form `varname[1:n,1:n]` or `varname[i=1:n,j=1:n]`,
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
macro variable(input_args...)
    error_fn = Containers.build_error_fn(:variable, input_args, __source__)
    args, kwargs = Containers.parse_macro_arguments(error_fn, input_args)
    if length(args) >= 2 && Meta.isexpr(args[2], :block)
        error_fn("Invalid syntax. Did you mean to use `@variables`?")
    end
    model_sym = popfirst!(args)
    model = esc(model_sym)
    x = isempty(args) ? nothing : popfirst!(args)
    if x == :Int
        error_fn(
            "Ambiguous variable name $x detected. To specify an anonymous " *
            "integer variable, use `@variable(model, integer = true)` instead.",
        )
    elseif x == :Bin
        error_fn(
            "Ambiguous variable name $x detected. To specify an anonymous " *
            "binary variable, use `@variable(model, binary = true)` instead.",
        )
    elseif x == :PSD
        error_fn(
            "Size of anonymous square matrix of positive semidefinite " *
            "anonymous variables is not specified. To specify size of square " *
            "matrix use `@variable(model, [1:n, 1:n], PSD)` instead.",
        )
    end
    info_kwargs =
        [(k, _esc_non_constant(v)) for (k, v) in kwargs if k in _INFO_KWARGS]
    info_expr = _VariableInfoExpr(; info_kwargs...)
    # There are four cases to consider:
    # x                                       | type of x | x.head
    # ----------------------------------------+-----------+------------
    # var                                     | Symbol    | NA
    # var[1:2]                                | Expr      | :ref
    # var <= ub or var[1:2] <= ub             | Expr      | :call
    # var in set or var[1:2] in set           | Expr      | :call
    # lb <= var <= ub or lb <= var[1:2] <= ub | Expr      | :comparison
    # In the three last cases, we call parse_variable
    var, set = x, nothing
    if Meta.isexpr(x, (:comparison, :call))
        var, set = parse_variable(error_fn, info_expr, x.args...)
        is_anonymous = Meta.isexpr(var, (:vect, :vcat)) || x === nothing
        if is_anonymous && set === nothing
            error_fn(
                "Cannot use explicit bounds via >=, <= with an anonymous variable",
            )
        end
    end
    # if var === nothing, then the variable is anonymous
    if !(var isa Union{Nothing,Symbol,Expr})
        error_fn("Expected $var to be a variable name")
    end
    name, index_vars, indices = Containers.parse_ref_sets(
        error_fn,
        var;
        invalid_index_variables = [model_sym],
    )
    # Handle special keyword arguments
    # ; set
    set_kw = get(kwargs, :set, nothing)
    if set_kw !== nothing
        if set !== nothing
            error_fn(
                "Cannot use set keyword because the variable is already " *
                "constrained to `$set`.",
            )
        end
        set = esc(set_kw)
    end
    # ; set_string_name
    name_expr = Containers.build_name_expr(name, index_vars, kwargs)
    if name_expr != ""
        set_string_name = if haskey(kwargs, :set_string_name)
            esc(kwargs[:set_string_name])
        else
            :(set_string_names_on_creation($model))
        end
        name_expr = :($set_string_name ? $name_expr : "")
    end
    # ; variable_type
    variable_type_kw = get(kwargs, :variable_type, nothing)
    if variable_type_kw !== nothing
        push!(args, variable_type_kw)
    end
    # Handle positional arguments
    for ex in args
        if ex == :Int
            _set_integer_or_error(error_fn, info_expr)
        elseif ex == :Bin
            _set_binary_or_error(error_fn, info_expr)
        elseif ex == :PSD
            if set !== nothing
                error_fn(
                    "Cannot pass `$ex` as a positional argument because the " *
                    "variable is already constrained to `$set`.",
                )
            end
            set = PSDCone()
        elseif ex == :Symmetric
            if set !== nothing
                error_fn(
                    "Cannot pass `$ex` as a positional argument because the " *
                    "variable is already constrained to `$set`.",
                )
            end
            set = SymmetricMatrixSpace()
        elseif ex == :Hermitian
            if set !== nothing
                error_fn(
                    "Cannot pass `$ex` as a positional argument because the " *
                    "variable is already constrained to `$set`.",
                )
            end
            set = HermitianMatrixSpace()
        end
    end
    filter!(ex -> !(ex in (:Int, :Bin, :PSD, :Symmetric, :Hermitian)), args)
    build_code = :(build_variable($error_fn, $(_constructor_expr(info_expr))))
    Containers.add_additional_args(
        build_code,
        args,
        kwargs;
        kwarg_exclude = vcat(
            _INFO_KWARGS,
            [:base_name, :container, :variable_type, :set, :set_string_name],
        ),
    )
    code = if set === nothing
        # This is for calls like:
        #   @variable(model, x)
        #   @variable(model, x[i in 1:2] >= i)
        Containers.container_code(
            index_vars,
            indices,
            quote
                variable = model_convert($model, $build_code)
                add_variable($model, variable, $name_expr)
            end,
            kwargs,
        )
    elseif any(Base.Fix2(Containers.depends_on, set), index_vars)
        # This is for calls in which the set depends on the indices.
        #   @variable(model, x[i in 1:2] in MOI.GreaterThan(i))
        Containers.container_code(
            index_vars,
            indices,
            quote
                build = build_variable($error_fn, $build_code, $set)
                add_variable($model, model_convert($model, build), $name_expr)
            end,
            kwargs,
        )
    else
        # This is for calls in which the set does not depend on the indices.
        #   @variable(model, x in MOI.GreaterThan(1.0))
        #   @variable(model, x[1:2, 1:2] in PSDCone())
        # We can't containerize the full expression; we need to keep the set
        # apart.
        build_code = Containers.container_code(
            index_vars,
            indices,
            build_code,
            kwargs,
        )
        name_expr = Containers.container_code(
            index_vars,
            indices,
            name_expr,
            kwargs,
        )
        quote
            build = build_variable($error_fn, $build_code, $set)
            add_variable($model, model_convert($model, build), $name_expr)
        end
    end
    return _finalize_macro(
        model,
        code,
        __source__;
        register_name = name,
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

_esc_non_constant(x::Number) = x
_esc_non_constant(x::Expr) = Meta.isexpr(x, :quote) ? x : esc(x)
_esc_non_constant(x) = esc(x)

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

"""
    reverse_sense(::Val{T}) where {T}

Given an (in)equality symbol `T`, return a new `Val` object with the opposite
(in)equality symbol.

This function is intended for use in JuMP extensions.

## Example

```jldoctest
julia> reverse_sense(Val(:>=))
Val{:<=}()
```
"""
function reverse_sense end
reverse_sense(::Val{:<=}) = Val(:>=)
reverse_sense(::Val{:≤}) = Val(:≥)
reverse_sense(::Val{:>=}) = Val(:<=)
reverse_sense(::Val{:≥}) = Val(:≤)
reverse_sense(::Val{:(==)}) = Val(:(==))

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

A hook for JuMP extensions to intercept the parsing of a `:comparison`
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
            """
            Unrecognized keyword argument: $key.

            The supported keyword arguments are:

             * `base_name`
             * `binary`
             * `container`
             * `integer`
             * `lower_bound`
             * `set`
             * `set_string_name`
             * `start`
             * `upper_bound`
             * `variable_type`

            If you're trying to create a JuMP extension, you need to implement `JuMP.build_variable`.
            See the docstring for more details.
            """,
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
