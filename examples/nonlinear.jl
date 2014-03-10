using JuMP
using Base.Test
# rosenbrock

let

    m = Model()

    @defVar(m, x)
    @defVar(m, y)

    @setNLObjective(m, Min, (1-x)^2 + 100(y-x^2)^2)

    solve(m)

    println("x = ", getValue(x), " y = ", getValue(y))

end

# clnlbeam
let
    N     = 1000
    ni    = N
    h     = 1/ni
    alpha = 350

    m = Model()

    @defVar(m, -1 <= t[1:(ni+1)] <= 1)
    @defVar(m, -0.05 <= x[1:(ni+1)] <= 0.05)
    @defVar(m, u[1:(ni+1)])

    @setNLObjective(m, Min, sum{ 0.5*h*(u[i+1]^2 + u[i]^2) + 0.5*alpha*h*(cos(t[i+1]) + cos(t[i])), i = 1:ni})
    
    # cons1
    for i in 1:ni
        @addNLConstraint(m, x[i+1] - x[i] - (0.5h)*(sin(t[i+1])+sin(t[i])) == 0)
    end
    # cons2
    for i in 1:ni
        @addConstraint(m, t[i+1] - t[i] - (0.5h)*u[i+1] - (0.5h)*u[i] == 0)
    end

    solve(m)

end
