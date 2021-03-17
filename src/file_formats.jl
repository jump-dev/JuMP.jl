#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
    dest = MOI.FileFormats.Model(format = format, filename = filename)
    # We add a `full_bridge_optimizer` here because MOI.FileFormats models may not
    # support all constraint types in a JuMP model.
    bridged_dest = MOI.Bridges.full_bridge_optimizer(dest, Float64)
    MOI.copy_to(bridged_dest, model)
    # `dest` will contain the underlying model, with constraints bridged if
    # necessary.
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
    dest = MOI.FileFormats.Model(format = format)
    # We add a `full_bridge_optimizer` here because MOI.FileFormats models may not
    # support all constraint types in a JuMP model.
    bridged_dest = MOI.Bridges.full_bridge_optimizer(dest, Float64)
    MOI.copy_to(bridged_dest, model)
    # `dest` will contain the underlying model, with constraints bridged if
    # necessary.
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
