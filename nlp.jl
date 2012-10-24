# NLP?


require("julp.jl")

m = Julp.Model("min")
x = Julp.Variable(m,"x",-Inf,Inf)
y = Julp.Variable(m,"y",-Inf,Inf)

con = :($x^4 + $x^2 <= 4)

println(con)
println(con.args[1])
println("Less than or equal to")
println(con.args[3])

lhs = con.args[1]
println(lhs.args[2].args)

function ChainRule(expr,var)
  println("In ChainRule")
  println(expr)
  if typeof(expr) == Julp.Variable
    println("Variable")
    if expr == var
      println("Symbol==Var, so is just 1")
      return 1
    else
      println("Symbol!=Var, so zero")
      return 0
    end
  end
  if typeof(expr) <: Number
    println("Constant, so zero")
    return 0
  end
  if expr.args[1] == :^
    println("^ operator, diff it")
    return :( $(expr.args[3]) * $(ChainRule(expr.args[2],var)) * ($(expr.args[2]) ^ ($(expr.args[3] - 1))) )
  end
  if expr.args[1] == :+
    println("+ operator, add diff of terms")
    return :( $(ChainRule(expr.args[2],var)) + $(ChainRule(expr.args[3],var)) )
  end
  if expr.args[1] == :-
    println("- operator, subtract diff of terms")
    return :( $(ChainRule(expr.args[2],var)) - $(ChainRule(expr.args[3],var)) )
  end
  if expr.args[1] == :*
    println("* operation")
    @assert length(expr.args) == 3
    return :( $(ChainRule(expr.args[2],var))*$(expr.args[3]) + $(expr.args[2])*$(ChainRule(expr.args[3],var)))
  end

  return 1
end

function substituteVariables(expr,values)
  if typeof(expr) == Expr
    for i in 1:length(expr.args)
      if typeof(expr.args[i]) == Julp.Variable
        expr.args[i] = values[expr.args[i].col]
      elseif typeof(expr.args[i]) == Expr
        substituteVariables(expr.args[i],values)
      end
    end
  end
  return expr
end

function evalAt(expr,values)
  sub = substituteVariables(copy(expr),values)
  eval(sub)
end

out = ChainRule(lhs,x)
println(out)
println("Rosenbrockin it")
rosenbrock = :( (1-$x)^2 + 100($y-$x^2)^2 )
dx = ChainRule(rosenbrock, x)
dy = ChainRule(rosenbrock, y)
println(rosenbrock)
println(dx)
println(dy)

# (1,1) is global minimizer
at = [1,1]
println("Rosenbrock function at $at:")
println("Obj: ",evalAt(rosenbrock,at))
println("Grad: ",[evalAt(dx,at),evalAt(dy,at)])

