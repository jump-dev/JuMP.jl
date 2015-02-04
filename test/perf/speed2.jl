#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# speed2.jl
#
# Runs some JuMP benchmarks to test for speed-related regressions.
# Test time to build model in memory.
#############################################################################

using JuMP

function pMedian(numFacility::Int,numCustomer::Int,numLocation::Int)
    srand(10)
    customerLocations = [rand(1:numLocation) for a = 1:numCustomer ]

    tic()
    m = Model()

    # Facility locations
    @defVar(m, 0 <= s[1:numLocation] <= 1)

    # Aux. variable: x_a,i = 1 iff the closest facility to a is at i
    @defVar(m, 0 <= x[1:numLocation,1:numCustomer] <= 1)

    # Objective: min distance
    @setObjective(m, Max, sum{abs(customerLocations[a]-i)*x[i,a], a = 1:numCustomer, i = 1:numLocation} )

    # Constraints
    for a in 1:numCustomer
        # Subject to linking x with s
        for i in 1:numLocation
            @addConstraint(m, x[i,a] - s[i] <= 0)
        end
        # Subject to one of x must be 1
        @addConstraint(m, sum{x[i,a],i=1:numLocation} == 1 )
    end

    # Subject to must allocate all facilities
    @addConstraint(m, sum{s[i],i=1:numLocation} == numFacility )
    buildTime = toq()

    tic()
    buildInternalModel(m)
    writeTime = toq()

    return buildTime, writeTime
end

function cont5(n)
    m = n
    n1 = n-1
    m1 = m-1
    dx = 1/n
    T = 1.58
    dt = T/m
    h2 = dx^2
    a = 0.001
    yt = [0.5*(1 - (j*dx)^2) for j=0:n]

    tic()
    mod = Model()
    @defVar(mod,  0 <= y[0:m,0:n] <= 1)
    @defVar(mod, -1 <= u[1:m] <= 1)
    @setObjective(mod, Min, 0.25*dx*( (y[m,0] - yt[1])^2 +
       2*sum{ (y[m,j]-yt[j+1])^2, j=1:n1} + (y[m,n]-yt[n+1])^2) +
       0.25*a*dt*(2*sum{u[i]^2,i=1:m1} + u[m]^2))

    # PDE
    for i = 0:m1
        for j = 1:n1
            @addConstraint(mod, h2*(y[i+1,j] - y[i,j]) == 0.5*dt*(y[i,j-1] - 2*y[i,j] + y[i,j+1] + y[i+1,j-1] - 2*y[i+1,j] + y[i+1,j+1]) )
        end
    end

    # IC
    for j = 0:n
        @addConstraint(mod, y[0,j] == 0)
    end

    # BC
    for i = 1:m
        @addConstraint(mod, y[i,2]   - 4*y[i,1]  + 3*y[i,0] == 0)
        @addConstraint(mod, y[i,n-2] - 4*y[i,n1] + 3*y[i,n] == (2*dx)*(u[i] - y[i,n]))
    end
    buildTime = toq()

    tic()
    buildInternalModel(mod)
    writeTime = toq()

    return buildTime, writeTime
end


function RunTests()
    # Pmedian
    pmedian_build = Float64[]
    pmedian_write = Float64[]
    for runs = 1:9
        bt, wt = pMedian(100,100,5000)
        push!(pmedian_build, bt)
        push!(pmedian_write, wt)
    end
    sort!(pmedian_build)
    sort!(pmedian_write)
    print("PMEDIAN BUILD MIN=",minimum(pmedian_build),"  MED=",pmedian_build[5],"\n")
    print("PMEDIAN INTRN MIN=",minimum(pmedian_write),"  MED=",pmedian_write[5],"\n")

    # Cont5
    cont5_build = Float64[]
    cont5_write = Float64[]
    for runs = 1:9
        bt, wt = cont5(500)
        push!(cont5_build, bt)
        push!(cont5_write, wt)
    end
    sort!(cont5_build)
    sort!(cont5_write)
    print("CONT5 BUILD   MIN=",minimum(cont5_build),"  MED=",cont5_build[5],"\n")
    print("CONT5 INTRN   MIN=",minimum(cont5_write),"  MED=",cont5_write[5],"\n")

end

RunTests()

