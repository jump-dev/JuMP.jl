# multivarate "dictionary" used for collections of variables/constraints
type MultivarDict{T,N}
  innerArray::Array{T,N}
  name::String
end

function ref(d::MultivarDict, vals...)
    d.innerArray[vals...]
end

function assign(d::MultivarDict, vals...)
    assign(d.innerArray,vals...)
end
