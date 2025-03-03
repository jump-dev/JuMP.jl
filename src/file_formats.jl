#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function _throw_write_to_file_explanatory_message(
    ::MOI.UnsupportedConstraint{F,S},
) where {F,S}
    return error(
        "Unable to write problem to file because the chosen file format " *
        "doesn't support constraints of the type $F-in-$S.",
    )
end

_throw_write_to_file_explanatory_message(err) = rethrow(err)

function _copy_to_bridged_model(f::Function, model::GenericModel{T}) where {T}
    inner = MOI.instantiate(f; with_bridge_type = T)
    try
        MOI.copy_to(inner, model)
    catch err
        _throw_write_to_file_explanatory_message(err)
    end
    @assert inner isa MOI.Bridges.LazyBridgeOptimizer
    if inner.model isa MOI.Utilities.CachingOptimizer
        MOI.Utilities.attach_optimizer(inner.model)
    end
    return unsafe_backend(inner)
end

"""
    write_to_file(
        model::GenericModel,
        filename::String;
        format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC,
        kwargs...,
    )

Write the JuMP model `model` to `filename` in the format `format`.

See [`MOI.FileFormats.FileFormat`](@ref) for a list of supported formats.

## Compression

If the filename ends in `.gz`, the file will be compressed using GZip.

If the filename ends in `.bz2`, the file will be compressed using BZip2.

## Keyword arguments

Other `kwargs` are passed to the `Model` constructor of the chosen format.

For details, see the docstring each file format's `Model` constructor. For
example, [`MOI.FileFormats.MPS.Model`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0);

julia> @objective(model, Min, 2 * x + 1);

julia> filename = joinpath(mktempdir(), "model.mps");

julia> write_to_file(model, filename; generic_names = true)

julia> print(read(filename, String))
NAME
ROWS
 N  OBJ
COLUMNS
    C1        OBJ       2
RHS
    rhs       OBJ       -1
RANGES
BOUNDS
 LO bounds    C1        0
 PL bounds    C1
ENDATA
```

## Solver-specific formats

[`write_to_file`](@ref) calls a Julia function from the MathOptInterface package
that is independent of the choice of solver. That is, it does not call a
solver's internal API to write a model to disk.

For MPS files in particular, [`write_to_file`](@ref) may not support the full
range of features that a solver's internal API supports. This is because some
solvers have defined solver-specific extensions to the MPS format, whereas our
Julia implementation supports only features which are standardized across
multiple solvers.

To write a file to disk using the solver's internal API, use
[`direct_model`](@ref) and call the solver's C API. For example:
```jldoctest
julia> import HiGHS

julia> model = direct_model(HiGHS.Optimizer());

julia> set_silent(model)

julia> @variable(model, x >= 0);

julia> @objective(model, Min, 2 * x + 1);

julia> filename = joinpath(mktempdir(), "model.mps");

julia> HiGHS.Highs_writeModel(backend(model), filename);

julia> print(read(filename, String))
NAME
ROWS
 N  Obj
COLUMNS
    c0        Obj       2
RHS
    RHS_V     Obj       -1
ENDATA
```

"""
function write_to_file(
    model::GenericModel,
    filename::String;
    format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC,
    kwargs...,
)
    dest = _copy_to_bridged_model(model) do
        return MOI.FileFormats.Model(;
            format = format,
            filename = filename,
            kwargs...,
        )
    end
    MOI.write_to_file(dest, filename)
    return
end

"""
    Base.write(
        io::IO,
        model::GenericModel;
        format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_MOF,
        kwargs...,
    )

Write the JuMP model `model` to `io` in the format `format`.

See [`MOI.FileFormats.FileFormat`](@ref) for a list of supported formats.

Other `kwargs` are passed to the `Model` constructor of the chosen format.

## Keyword arguments

Other `kwargs` are passed to the `Model` constructor of the chosen format.

For details, see the docstring each file format's `Model` constructor. For
example, [`MOI.FileFormats.MPS.Model`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0);

julia> @objective(model, Min, 2 * x + 1);

julia> io = IOBuffer();

julia> write(io, model; format = MOI.FileFormats.FORMAT_MPS);

julia> seekstart(io);

julia> print(read(io, String))
NAME
ROWS
 N  OBJ
COLUMNS
    x         OBJ       2
RHS
    rhs       OBJ       -1
RANGES
BOUNDS
 LO bounds    x         0
 PL bounds    x
ENDATA
```
"""
function Base.write(
    io::IO,
    model::GenericModel;
    format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_MOF,
    kwargs...,
)
    if format == MOI.FileFormats.FORMAT_AUTOMATIC
        error("Unable to infer the file format from an IO stream.")
    end
    dest = _copy_to_bridged_model(model) do
        return MOI.FileFormats.Model(; format = format, kwargs...)
    end
    write(io, dest)
    return
end

_value_type(::MOI.Utilities.AbstractModelLike{T}) where {T} = T

# This fallback may not get the correct value type. However, since
# all models defined in `MOI.FileFormats` return a
# `MOI.Utilities.GenericModel` except `NL` and `MOF` which only supports
# `Float64`, this does the job for now.
_value_type(::MOI.ModelLike) = Float64

"""
    read_from_file(
        filename::String;
        format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC,
        kwargs...,
    )

Return a JuMP model read from `filename` in the format `format`.

See [`MOI.FileFormats.FileFormat`](@ref) for a list of supported formats.

## Compression

If the filename ends in `.gz`, the file will be uncompressed using GZip.

If the filename ends in `.bz2`, the file will be uncompressed using BZip2.

## Keyword arguments

Other `kwargs` are passed to the `Model` constructor of the chosen format.

For details, see the docstring each file format's `Model` constructor. For
example, [`MOI.FileFormats.MPS.Model`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0);

julia> @objective(model, Min, 2 * x + 1);

julia> filename = joinpath(mktempdir(), "model.mps");

julia> write_to_file(model, filename; generic_names = true)

julia> new_model = read_from_file(filename)
A JuMP Model
├ solver: none
├ objective_sense: MIN_SENSE
│ └ objective_function_type: AffExpr
├ num_variables: 1
├ num_constraints: 1
│ └ VariableRef in MOI.GreaterThan{Float64}: 1
└ Names registered in the model: none

julia> print(new_model)
Min 2 C1 + 1
Subject to
 C1 ≥ 0
```
"""
function read_from_file(
    filename::String;
    format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC,
    kwargs...,
)
    src = MOI.FileFormats.Model(; format, filename, kwargs...)
    MOI.read_from_file(src, filename)
    model = GenericModel{_value_type(src)}()
    MOI.copy_to(model, src)
    return model
end

"""
    Base.read(
        io::IO,
        ::Type{<:GenericModel};
        format::MOI.FileFormats.FileFormat,
        kwargs...,
    )

Return a JuMP model read from `io` in the format `format`.

See [`MOI.FileFormats.FileFormat`](@ref) for a list of supported formats.

## Keyword arguments

Other `kwargs` are passed to the `Model` constructor of the chosen format.

For details, see the docstring each file format's `Model` constructor. For
example, [`MOI.FileFormats.MPS.Model`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0);

julia> @objective(model, Min, 2 * x + 1);

julia> io = IOBuffer();

julia> write(io, model; format = MOI.FileFormats.FORMAT_MPS);

julia> seekstart(io);

julia> new_model = read(io, Model; format = MOI.FileFormats.FORMAT_MPS)
A JuMP Model
├ solver: none
├ objective_sense: MIN_SENSE
│ └ objective_function_type: AffExpr
├ num_variables: 1
├ num_constraints: 1
│ └ VariableRef in MOI.GreaterThan{Float64}: 1
└ Names registered in the model: none

julia> print(new_model)
Min 2 x + 1
Subject to
 x ≥ 0
```
"""
function Base.read(
    io::IO,
    ::Type{GenericModel{T}};
    format::MOI.FileFormats.FileFormat,
    kwargs...,
) where {T}
    if format == MOI.FileFormats.FORMAT_AUTOMATIC
        error("Unable to infer the file format from an IO stream.")
    end
    src = MOI.FileFormats.Model(; format = format, kwargs...)
    read!(io, src)
    model = GenericModel{T}()
    MOI.copy_to(model, src)
    return model
end
