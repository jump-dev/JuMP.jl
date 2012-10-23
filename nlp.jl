# NLP?

type Variable
  name::String
end

x = Variable("x")

con = :(x^4 + x^2 <= 4)

println(con)
println(con.args[1])
println("Less than or equal to")
println(con.args[3])

lhs = con.args[1]
println(lhs.args[2].args)

function ChainRule(expr,var)
  println("In ChainRule")
  println(expr)
  if typeof(expr) == Symbol
    println("Symbol")
    if eval(expr) == var
      println("Symbol==Var, so is just 1")
      return 1
    else
      println("Symbol!=Var, so return it")
      return expr
    end
  end
  if expr.args[1] == :^
    println("^ operator, diff it")
    return :( $(expr.args[3]) * $(ChainRule(expr.args[2],var)) * ($(expr.args[2]) ^ ($(expr.args[3]) - 1)) )
  end
  if expr.args[1] == :+
    println("+ operator, add diff of terms")
    return :( $(ChainRule(expr.args[2],var)) + $(ChainRule(expr.args[3],var)) )
  end
  return 1
end
out = ChainRule(lhs,x)
println(out)

