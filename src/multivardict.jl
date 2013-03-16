# multivarate "dictionary" used for collections of variables/constraints
type MultivarDict{T,N}
  innerArray::Array{T,N}
  name::String
end

getindex(d::MultivarDict, vals...) = getindex(d.innerArray, vals...)

setindex!(d::MultivarDict, vals...) = setindex!(d.innerArray, vals...)
