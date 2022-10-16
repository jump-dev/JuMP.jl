function _row_iterator(x::Union{Array,Containers.DenseAxisArray})
    return zip(eachindex(x), Iterators.product(axes(x)...))
end

function _row_iterator(x::Containers.SparseAxisArray)
    return zip(eachindex(x.data), keys(x.data))
end

_columns(x::Union{Array, Containers.DenseAxisArray}) = length(axes(x))
_columns(x::Containers.SparseAxisArray{T,N,K}) where {T,N,K} = N


function table(x, name::Symbol, col_names::Symbol...)
    return table(identity, x, name, col_names...)
end

"""
    table(f::Function, x, name, colnames...)

Applies the function `f` to all elements of the variable container `x`
and returns the result as a `Vector` of `NamedTuple`s using the provided
names.

A `Vector` of `NamedTuple`s implements the 'Tables.jl' interface
and can be used as input for anything that consumes a 'Tables.jl' 
compatible source. 

## Example
```julia
model = Model()
@variable(model, x[1:10, 2000:2020] >= 0)
[...]
optimize!(model)
tbl = table(value, x, :solution_value, :car, :year)
```
"""
function table(
    f::Function, 
    x::Union{Array,Containers.DenseAxisArray,Containers.SparseAxisArray}, 
    name::Symbol, 
    col_names::Symbol...,
)
    if length(col_names) < _columns(x)
        error("Not enough column names provided")
    end

    C = (col_names..., name)
    return vec([NamedTuple{C}((args..., f(x[i]))) for (i, args) in _row_iterator(x)])
end
