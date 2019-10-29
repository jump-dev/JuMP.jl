"List of accepted export formats. `Automatic` corresponds to a detection from the file name."
@enum(FileFormat, CBF, LP, MOF, MPS, Automatic)

function _filename_to_format(filename::String)
    return if endswith(filename, ".mof.json.gz") || endswith(filename, ".mof.json")
        MOF
    elseif endswith(filename, ".cbf.gz") || endswith(filename, ".cbf")
        CBF
    elseif endswith(filename, ".mps.gz") || endswith(filename, ".mps")
        MPS
    elseif endswith(filename, ".lp.gz") || endswith(filename, ".lp")
        LP
    else
        error("File type of $(filename) not recognized by JuMP.")
    end
end

"List of accepted export compression formats. `Automatic` corresponds to a detection from the file name."
@enum(FileCompression, None, GZip, Automatic)

function _filename_to_compression(filename::String)
    return if endswith(filename, ".gz")
        GZip
    else
        None
    end
end

function _open(f::Function, filename::String, mode::String; compression::FileCompression=Automatic)
    if compression == Automatic
        compression = _filename_to_compression(filename)
    end

    open_f = if compression == GZip
        GZip.open
    else
        open
    end
    return open_f(f, filename, mode)
end

function write_to_file(model::Model, io::IO, format::FileFormat=MPS; kwargs...)
    dest = if format == CBF
        MathOptFormat.CBF.Model(; kwargs...)
    elseif format == LP
        MathOptFormat.LP.Model(; kwargs...)
    elseif format == MOF
        MathOptFormat.MOF.Model(; kwargs...)
    elseif format == MPS
        MathOptFormat.MPS.Model(; kwargs...)
    else
        error("When passing an IO object to write_to_file, the file format cannot be guessed from the file name.")
    end
    MOI.copy_to(dest, backend(model))
    MOI.write_to_file(dest, io)
    return
end

"""
    write_to_file(model::Model, filename::String; format::FileFormat=Automatic, compression::FileCompression=Automatic, kwargs...)

Write `model` to the file called `filename` using the format `format`.

Valid formats are given by the enum [`FileFormat`](@ref). Valid compression algorithms
are given by the enum [`FileCompression`](@ref).

For keyword options, see [MathOptFormat.jl](https://github.com/odow/MathOptFormat.jl).
"""
function write_to_file(model::Model, filename::String; format::FileFormat=Automatic, compression::FileCompression=Automatic, kwargs...)
    if format == Automatic
        format = _filename_to_format(filename)
    end

    _open(filename, "w", compression) do io
        write_to_file(model, io, format; kwargs...)
    end
    return
end
