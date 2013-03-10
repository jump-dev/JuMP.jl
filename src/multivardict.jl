# multivarate "dictionary" used for collections of variables/constraints
type MultivarDict{T,N}
  innerArray::Array{T,N}
  name::String
end

ref(d::MultivarDict, vals...) = ref(d.innerArray, vals...)

assign(d::MultivarDict, vals...) = assign(d.innerArray, vals...)
