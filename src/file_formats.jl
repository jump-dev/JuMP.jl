function write_to_file(model::Model, io::IO, format::MathOptFormat.FileFormat; kwargs...)
    @assert format != AUTOMATIC_FILE_FORMAT
    dest = MathOptFormat._file_formats[format]()
    MOI.copy_to(dest, backend(model))
    MOI.write_to_file(dest, io)
    return
end

"""
    write_to_file(model::Model, filename::String; format::FileFormat=AUTOMATIC_FILE_FORMAT, compression::FileCompression=AUTOMATIC_FILE_COMPRESSION, kwargs...)

Write `model` to the file called `filename` using the format `format`.

Valid formats are given by the enum [`FileFormat`](@ref). Valid compression algorithms
are given by the enum [`FileCompression`](@ref).

For keyword options, see [MathOptFormat.jl](https://github.com/odow/MathOptFormat.jl).
"""
function write_to_file(model::Model, filename::String;
                       format::MathOptFormat.FileFormat=MathOptFormat.AUTOMATIC_FILE_FORMAT,
                       compression::MathOptFormat.FileCompression=AUTOMATIC_FILE_COMPRESSION,
                       kwargs...)
    if format == AUTOMATIC_FILE_FORMAT
        format = _filename_to_format(filename)
    end

    MathOptFormat.gzip_open(filename, "w", compression=compression) do io
        write_to_file(model, io, format; kwargs...)
    end
    return
end
