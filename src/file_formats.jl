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

If the filename ends in `.gz`, it will be compressed using Gzip.
If the filename ends in `.bz2`, it will be compressed using BZip2.

Other `kwargs` are passed to the `Model` constructor of the chosen format.
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

Other `kwargs` are passed to the `Model` constructor of the chosen format.
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

"""
    read_from_file(
        filename::String;
        format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC,
        kwargs...,
    )

Return a JuMP model read from `filename` in the format `format`.

If the filename ends in `.gz`, it will be uncompressed using Gzip.
If the filename ends in `.bz2`, it will be uncompressed using BZip2.

Other `kwargs` are passed to the `Model` constructor of the chosen format.
"""
function read_from_file(
    filename::String;
    format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC,
    kwargs...,
)
    src =
        MOI.FileFormats.Model(; format = format, filename = filename, kwargs...)
    MOI.read_from_file(src, filename)
    # TODO(odow): what number type to choose? Are there any non-Float64 file
    # formats?
    model = GenericModel{Float64}()
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

Other `kwargs` are passed to the `Model` constructor of the chosen format.
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
