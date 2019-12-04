#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

const MOF = MathOptFormat

"""
    write_to_file(
        model::Model,
        filename::String;
        compression::MOF.AbstractCompressionScheme = MOF.AutomaticCompression(),
        format::MOF.FileFormat == MOF.AUTOMATIC_FILE_FORMAT
    )

Write the JuMP model `model` to `filename` in the format `format` and compress
the file using the compression scheme `compression`.

Supported formats include:
 - `format = MOF.FORMAT_CBF`. Inferred if filename ends in `.cbf`.
 - `format = MOF.FORMAT_MOF`. Inferred if filename ends in `.mof.json`.
 - `format = MOF.FORMAT_MPS`. Inferred if filename ends in `.mps`.
 - `format = MOF.FORMAT_LP`. Inferred if filename ends in `.lp`.

 See `MathOptFormat.FileFormat` for all supported formats.

Supported compression formats include:
 - `compression = MOF.NoCompression()`
 - `compression = MOF.Gzip()`. Inferred if filename ends in `.gz`.
 - `compression = MOF.Bzip2()`. Inferred if filename ends in `.bz2`.

Default arguments will infer the compression type and format from the filename.
"""
function write_to_file(
    model::Model,
    filename::String;
    compression::MOF.AbstractCompressionScheme = MOF.AutomaticCompression(),
    format::MOF.FileFormat = MOF.AUTOMATIC_FILE_FORMAT,
)
    if format == MOF.AUTOMATIC_FILE_FORMAT
        format = MOF._detect_file_format(filename)
    end
    dest = MOF._FILE_FORMATS[format][2]()
    bridged_dest = MOI.Bridges.full_bridge_optimizer(dest, Float64)
    MOI.copy_to(bridged_dest, backend(model))
    MOI.write_to_file(dest, filename; compression = compression)
    return
end

"""
    write_to_file(model::Model, io::IO; format::MOF.FileFormat)

Write the JuMP model `model` to `io` in the format `format`.

Supported formats include:
 - `format = MOF.FORMAT_CBF`
 - `format = MOF.FORMAT_MOF`
 - `format = MOF.FORMAT_MPS`
 - `format = MOF.FORMAT_LP`

See `MathOptFormat.FileFormat` for all supported formats.
"""
function write_to_file(model::Model, io::IO; format::MOF.FileFormat)
    if format == MOF.AUTOMATIC_FILE_FORMAT
        error("Unable to infer the file format from an IO stream.")
    end
    dest = MOF._FILE_FORMATS[format][2]()
    bridged_dest = MOI.Bridges.full_bridge_optimizer(dest, Float64)
    MOI.copy_to(bridged_dest, backend(model))
    MOI.write_to_file(dest, io)
    return
end

"""
    read_from_file(
        filename::String;
        compression::MOF.AbstractCompressionScheme = MOF.AutomaticCompression(),
        format::MOF.FileFormat == MOF.AUTOMATIC_FILE_FORMAT
    )

Return a JuMP model read from `filename` in the format `format` and de-compress
the file using the compression scheme `compression`.

Supported formats include:
 - `format = MOF.FORMAT_CBF`. Inferred if filename ends in `.cbf`.
 - `format = MOF.FORMAT_MOF`. Inferred if filename ends in `.mof.json`.
 - `format = MOF.FORMAT_MPS`. Inferred if filename ends in `.mps`.
 - `format = MOF.FORMAT_LP`. Inferred if filename ends in `.lp`.

See `MathOptFormat.FileFormat` for all supported formats.

Supported compression formats include:
 - `compression = MOF.NoCompression()`
 - `compression = MOF.Gzip()`. Inferred if filename ends in `.gz`.
 - `compression = MOF.Bzip2()`. Inferred if filename ends in `.bz2`.

Default arguments will infer the compression type and format from the filename.
"""
function read_from_file(
    filename::String;
    compression::MOF.AbstractCompressionScheme = MOF.AutomaticCompression(),
    format::MOF.FileFormat = MOF.AUTOMATIC_FILE_FORMAT,
)
    src = MathOptFormat.read_from_file(
        filename; compression = compression, file_format = format
    )
    model = Model()
    MOI.copy_to(backend(model), src)
    return model
end

"""
    read_from_file(io::IO; format::MOF.FileFormat)

Return a JuMP model read from `io` in the format `format`.

Supported formats include:
 - `format = MOF.FORMAT_CBF`
 - `format = MOF.FORMAT_MOF`
 - `format = MOF.FORMAT_MPS`
 - `format = MOF.FORMAT_LP`

 See `MathOptFormat.FileFormat` for all supported formats.
"""
function read_from_file(io::IO; format::MOF.FileFormat)
    if format == MOF.AUTOMATIC_FILE_FORMAT
        error("Unable to infer the file format from an IO stream.")
    end
    src = MOF._FILE_FORMATS[format][2]()
    MOI.read_from_file(src, io)
    model = Model()
    MOI.copy_to(backend(model), src)
    return model
end
