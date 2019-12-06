#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    FileFormat

A list of supported file formats for reading and writing JuMP models.
"""
@enum(
    FileFormat,
    FILE_FORMAT_AUTOMATIC,
    FILE_FORMAT_CBF,
    FILE_FORMAT_MOF,
    FILE_FORMAT_MPS,
    FILE_FORMAT_LP,
)

const _FORMAT_MAPPING = Dict(
    FILE_FORMAT_AUTOMATIC => MathOptFormat.AUTOMATIC_FILE_FORMAT,
    FILE_FORMAT_CBF => MathOptFormat.FORMAT_CBF,
    FILE_FORMAT_MOF => MathOptFormat.FORMAT_MOF,
    FILE_FORMAT_MPS => MathOptFormat.FORMAT_MPS,
    FILE_FORMAT_LP => MathOptFormat.FORMAT_LP,
)

"""
    write_to_file(
        model::Model,
        filename::String;
        format::FileFormat = FILE_FORMAT_AUTOMATIC,
        compression::MathOptFormat.AbstractCompressionScheme =
            MathOptFormat.AutomaticCompression(),
    )

Write the JuMP model `model` to `filename` in the format `format` and compress
the file using the compression scheme `compression`.

Supported compression formats include:
 - `compression = MathOptFormat.NoCompression()`
 - `compression = MathOptFormat.Gzip()`. Inferred if filename ends in `.gz`.
 - `compression = MathOptFormat.Bzip2()`. Inferred if filename ends in `.bz2`.

Default arguments will infer the compression type and format from the filename.
"""
function write_to_file(
    model::Model,
    filename::String;
    format::FileFormat = FILE_FORMAT_AUTOMATIC,
    compression::MathOptFormat.AbstractCompressionScheme =
        MathOptFormat.AutomaticCompression()
)
    if format == FILE_FORMAT_AUTOMATIC
        format = MathOptFormat._detect_file_format(filename)
    else
        format = _FORMAT_MAPPING[format]
    end
    dest = MathOptFormat._FILE_FORMATS[format][2]()
    bridged_dest = MOI.Bridges.full_bridge_optimizer(dest, Float64)
    MOI.copy_to(bridged_dest, backend(model))
    MOI.write_to_file(dest, filename; compression = compression)
    return
end

"""
    Base.write(io::IO, model::Model; format::FileFormat)

Write the JuMP model `model` to `io` in the format `format`.
"""
function Base.write(io::IO, model::Model; format::FileFormat)
    if format == FILE_FORMAT_AUTOMATIC
        error("Unable to infer the file format from an IO stream.")
    end
    dest = MathOptFormat._FILE_FORMATS[_FORMAT_MAPPING[format]][2]()
    bridged_dest = MOI.Bridges.full_bridge_optimizer(dest, Float64)
    MOI.copy_to(bridged_dest, backend(model))
    MOI.write_to_file(dest, io)
    return
end

"""
    read_from_file(
        filename::String;
        format::FileFormat = FILE_FORMAT_AUTOMATIC,
        compression::MathOptFormat.AbstractCompressionScheme =
            MathOptFormat.AutomaticCompression(),
    )

Return a JuMP model read from `filename` in the format `format` and de-compress
the file using the compression scheme `compression`.

Supported compression formats include:
 - `compression = MathOptFormat.NoCompression()`
 - `compression = MathOptFormat.Gzip()`. Inferred if filename ends in `.gz`.
 - `compression = MathOptFormat.Bzip2()`. Inferred if filename ends in `.bz2`.

Default arguments will infer the compression type and format from the filename.
"""
function read_from_file(
    filename::String;
    format::FileFormat = FILE_FORMAT_AUTOMATIC,
    compression::MathOptFormat.AbstractCompressionScheme =
        MathOptFormat.AutomaticCompression(),
)
    src = MathOptFormat.read_from_file(
        filename;
        file_format = _FORMAT_MAPPING[format],
        compression = compression
    )
    model = Model()
    MOI.copy_to(backend(model), src)
    return model
end

"""
    Base.read(io::IO, ::Type{Model}; format::FileFormat)

Return a JuMP model read from `io` in the format `format`.
"""
function Base.read(io::IO, ::Type{Model}; format::FileFormat)
    if format == FILE_FORMAT_AUTOMATIC
        error("Unable to infer the file format from an IO stream.")
    end
    src = MathOptFormat._FILE_FORMATS[_FORMAT_MAPPING[format]][2]()
    MOI.read_from_file(src, io)
    model = Model()
    MOI.copy_to(backend(model), src)
    return model
end
