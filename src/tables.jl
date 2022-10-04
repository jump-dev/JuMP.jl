abstract type _SolutionTable end

Tables.istable(::Type{<:_SolutionTable}) = true
Tables.rowaccess(::Type{<:_SolutionTable}) = true

_column_names(t::_SolutionTable) = getfield(t, :column_names)
_lookup(t::_SolutionTable) = getfield(t, :lookup)

Base.eltype(::_SolutionTable) = _SolutionRow
Base.length(t::_SolutionTable) = length(t.var)

struct _SolutionRow <: Tables.AbstractRow
    index_vals::Any
    sol_val::Number
    source::_SolutionTable
end

function Tables.getcolumn(s::_SolutionRow, i::Int)
    if i > length(getfield(s, :index_vals))
        return getfield(s, :sol_val)
    end
    return getfield(s, :index_vals)[i]
end

function Tables.getcolumn(s::_SolutionRow, nm::Symbol)
    i = _lookup(getfield(s, :source))[nm]
    if i > length(getfield(s, :index_vals))
        return getfield(s, :sol_val)
    end
    return getfield(s, :index_vals)[i]
end

Tables.columnnames(s::_SolutionRow) = _column_names(getfield(s, :source))

struct _SolutionTableDense{C} <: _SolutionTable
    column_names::Vector{Symbol}
    lookup::Dict{Symbol,Int}
    index_lookup::Dict
    var::C
end

function Base.iterate(t::_SolutionTableDense, state = nothing)
    next =
        isnothing(state) ? iterate(CartesianIndices(t.var)) :
        iterate(CartesianIndices(t.var), state)
    next === nothing && return nothing
    index = next[1]
    index_vals = [t.index_lookup[i][index[i]] for i in 1:length(index)]
    return _SolutionRow(index_vals, JuMP.value(t.var[next[1]]), t), next[2]
end

function _SolutionTableDense(
    v::Containers.DenseAxisArray{T,N,Ax,L},
    name,
    colnames...,
) where {T<:AbstractVariableRef,N,Ax,L}
    if length(colnames) < N
        error("Not enough column names provided")
    end
    if length(v) > 0 && !has_values(owner_model(first(v)))
        error("No solution values available for variable")
    end
    all_names = vcat(colnames..., name)
    lookup = Dict(nm => i for (i, nm) in enumerate(all_names))
    index_lookup = Dict()
    for (i, ax) in enumerate(axes(v))
        index_lookup[i] = collect(ax)
    end
    return _SolutionTableDense(all_names, lookup, index_lookup, v)
end

"""
    solution_table(var::DenseAxisArray, name, colnames...)

Returns the solution values of the variable container `var` as a table
that implements the `Tables.jl` interface. 

The table will have one column for each index and a column with the
corresponding solution value. The name of the column with the solution 
value is provided by `name`, while `colnames` provides the name of the 
index columns.

## Example
```julia
model = Model()
@variable(model, x[1:10, 2000:2020] >= 0)
[...]
optimize!(model)
tbl = solution_table(x, :value, :car, :year)
```
"""
function solution_table(
    var::Containers.DenseAxisArray{T,N,Ax,L},
    name,
    colnames...,
) where {T<:AbstractVariableRef,N,Ax,L}
    return _SolutionTableDense(var, name, colnames...)
end

function _SolutionTableDense(
    v::Array{T},
    name,
    colnames...,
) where {T<:AbstractVariableRef}
    if length(colnames) < length(axes(v))
        error("Not enough column names provided")
    end
    if length(v) > 0 && !has_values(owner_model(first(v)))
        error("No solution values available for variable")
    end
    all_names = vcat(colnames..., name)
    lookup = Dict(nm => i for (i, nm) in enumerate(all_names))
    index_lookup = Dict()
    for (i, ax) in enumerate(axes(v))
        index_lookup[i] = collect(ax)
    end
    return _SolutionTableDense(all_names, lookup, index_lookup, v)
end

function solution_table(
    var::Array{T},
    name,
    colnames...,
) where {T<:AbstractVariableRef}
    return _SolutionTableDense(var, name, colnames...)
end

struct _SolutionTableSparse <: _SolutionTable
    column_names::Vector{Symbol}
    lookup::Dict{Symbol,Int}
    var::Containers.SparseAxisArray
end

function _SolutionTableSparse(
    v::Containers.SparseAxisArray{T,N,K},
    name,
    colnames...,
) where {T<:AbstractVariableRef,N,K}
    if length(colnames) < N
        error("Not enough column names provided")
    end
    if length(v) > 0 && !has_values(first(v).model)
        error("No solution values available for variable")
    end
    all_names = vcat(colnames..., name)
    lookup = Dict(nm => i for (i, nm) in enumerate(all_names))
    return _SolutionTableSparse(all_names, lookup, v)
end

function Base.iterate(t::_SolutionTableSparse, state = nothing)
    next =
        isnothing(state) ? iterate(eachindex(t.var)) :
        iterate(eachindex(t.var), state)
    next === nothing && return nothing
    return _SolutionRow(next[1], JuMP.value(t.var[next[1]]), t), next[2]
end

function solution_table(
    var::Containers.SparseAxisArray{T,N,K},
    name,
    colnames...,
) where {T<:AbstractVariableRef,N,K}
    return _SolutionTableSparse(var, name, colnames...)
end
