using JuMP
# rosenbrock

let

    m = Model()

    @defVar(m, x)
    @defVar(m, y)

    @setNLObjective(m, Min, (1-x)^2 + 100(y-x^2)^2)

    solve(m)

    println("x = ", getValue(x), " y = ", getValue(y))

end
