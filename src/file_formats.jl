#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function _copy_to_bridged_model(f::Function, model::Model)
    inner = MOI.instantiate(f; with_bridge_type = Float64)
    MOI.copy_to(inner, model)
    @assert inner isa MOI.Bridges.LazyBridgeOptimizer
    if inner.model isa MOI.Utilities.CachingOptimizer
        MOI.Utilities.attach_optimizer(inner.model)
    end
    return unsafe_backend(inner)
end

"""
    write_to_file(
        model::Model,
        filename::String;
        format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC
    )

Write the JuMP model `model` to `filename` in the format `format`.

If the filename ends in `.gz`, it will be compressed using Gzip.
If the filename ends in `.bz2`, it will be compressed using BZip2.
"""
function write_to_file(
    model::Model,
    filename::String;
    format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC,
)
    dest = _copy_to_bridged_model(model) do
        return MOI.FileFormats.Model(format = format, filename = filename)
    end
    MOI.write_to_file(dest, filename)
    return
end

"""
    Base.write(
        io::IO,
        model::Model;
        format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_MOF
    )

Write the JuMP model `model` to `io` in the format `format`.
"""
function Base.write(
    io::IO,
    model::Model;
    format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_MOF,
)
    if format == MOI.FileFormats.FORMAT_AUTOMATIC
        error("Unable to infer the file format from an IO stream.")
    end
    dest = _copy_to_bridged_model(model) do
        return MOI.FileFormats.Model(format = format)
    end
    write(io, dest)
    return
end

"""
    read_from_file(
        filename::String;
        format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC
    )

Return a JuMP model read from `filename` in the format `format`.

If the filename ends in `.gz`, it will be uncompressed using Gzip.
If the filename ends in `.bz2`, it will be uncompressed using BZip2.
"""
function read_from_file(
    filename::String;
    format::MOI.FileFormats.FileFormat = MOI.FileFormats.FORMAT_AUTOMATIC,
)
    src = MOI.FileFormats.Model(format = format, filename = filename)
    MOI.read_from_file(src, filename)
    model = Model()
    MOI.copy_to(model, src)
    return model
end

"""
    Base.read(io::IO, ::Type{Model}; format::MOI.FileFormats.FileFormat)

Return a JuMP model read from `io` in the format `format`.
"""
function Base.read(io::IO, ::Type{Model}; format::MOI.FileFormats.FileFormat)
    if format == MOI.FileFormats.FORMAT_AUTOMATIC
        error("Unable to infer the file format from an IO stream.")
    end
    src = MOI.FileFormats.Model(format = format)
    read!(io, src)
    model = Model()
    MOI.copy_to(model, src)
    return model
end
