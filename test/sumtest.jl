type Variable
  i::Integer
end

type AffExpr
  data::Array{(Variable,Float64),1}
end

vars = [Variable(i) for i = 1:10]

macro SumExpr(expr)
  #:(AffExpr( [ ($(expr.args[1].args[3]), $(expr.args[1].args[2])) for $(expr.args[2]) ] ))
  Expr(:comprehension,convert(Vector{Any},[:($(expr.args[1].args[3]), $(expr.args[1].args[2]) ),expr.args[2] ]),Any)
end


#macro rewrite(expr)
#  :[$(expr.args[1]) for $(expr.args[2])]
#end

macro rewrite2(expr)
  #Expr(:comprehension,expr.args,Any)
  Expr(:comprehension,convert(Vector{Any},[expr.args[1],expr.args[2]]),Any)
end 

@rewrite2 [i*i for i in 1:10]

#x = AffExpr( [ (vars[i],2.) for i in 1:10 ] ) 
y = @SumExpr([2.*vars[i] for i in 1:10])
show(y)
