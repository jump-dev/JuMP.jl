using JuMP
# using Clp
# # s = ['b', 'c']
#
# m = Model(with_optimizer(Clp.Optimizer))

m = Model()
@variable(m, 0 <= x[ 2:3, [:a, :b] ] <= 1)
@show typeof(x.data)
@show typeof(x.axes)
@show typeof(x.lookup)
@show supertype(typeof(x))

@show keys(x)

for k in keys(x)
    @show x[k]
end


#     struct CartesianIndex{N} <: AbstractCartesianIndex{N}
#         I::NTuple{N,Int}
#         CartesianIndex{N}(index::NTuple{N,Integer}) where {N} = new(index)
#     end
#
# struct CartesianIndices{N,R<:NTuple{N,AbstractUnitRange{Int}}} <: AbstractArray{CartesianIndex{N},N}
#     indices::R
# end

# struct CartesianIndex{N}
#     I::NTuple{N,Int}
# end
#
# struct GCartesianIndex{N, T}
#     I::NTuple{N,T}
# end
#
# @show CartesianIndex((1,2))
# @show CartesianIndex((:a,:b,:c))
