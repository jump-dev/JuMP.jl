ispsd(x::Matrix) = minimum(eigvals(x)) ≥ -1e-3
ispsd(x::JuMP.JuMPArray) = ispsd(x.innerArray)

facts("[sdp] Test simple SDP") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, X[1:3,1:3], SDP)
    @addSDPConstraint(m, X <= 1/2*eye(3,3))
    @defVar(m, Y[1:5,1:5], Symmetric)
    @addSDPConstraint(m, -ones(5,5) <= Y)
    @addSDPConstraint(m, Y <= 2*ones(5,5))
    @defVar(m, Z[1:4,1:4], Symmetric)
    @addSDPConstraint(m, ones(4,4) >= Z)

    @addConstraint(m, trace(X) == 1)
    @addConstraint(m, trace(Y) == 3)
    @addConstraint(m, trace(Z) == -1)
    @setObjective(m, Max, X[1,2] + Y[1,2] + Z[1,2])
    solve(m)

    XX, YY, ZZ, = getValue(X), getValue(Y), getValue(Z)
    Xtrue = [0.25 0.25 0
             0.25 0.25 0
             0    0    0.5]
    Ytrue = fill(0.6, 5, 5)
    Ztrue = [-1.5  3.5 1 1
              3.5 -1.5 1 1
              1    1   1 1
              1    1   1 1]

    @fact norm(XX-Xtrue) --> roughly(0, 1e-2)
    @fact norm(YY-Ytrue) --> roughly(0, 1e-2)
    @fact norm(ZZ-Ztrue) --> roughly(0, 1e-2)
    @fact ispsd(XX) --> true
    @fact ispsd(1/2*eye(3,3)-XX) --> true
    @fact ispsd(YY+ones(5,5)) --> true
    @fact ispsd(2*ones(5,5)-YY) --> true
    @fact ispsd(ones(4,4)-ZZ) --> true
    @fact trace(XX) --> roughly( 1, 1e-4)
    @fact trace(YY) --> roughly( 3, 1e-4)
    @fact trace(ZZ) --> roughly(-1, 1e-4)
    @fact getObjectiveValue(m) --> roughly(4.35, 1e-2)

    @setObjective(m, Min, X[1,2] + Y[1,2] + Z[1,2])
    solve(m)

    XX, YY, ZZ, = getValue(X), getValue(Y), getValue(Z)
    Xtrue = [ 0.25 -0.25 0
             -0.25  0.25 0
              0     0    0.5]
    Ytrue = fill(0.6, 5, 5)
    Ztrue = [-1.5 -1.5 1 1
             -1.5 -1.5 1 1
              1    1   1 1
              1    1   1 1]

    @fact norm(XX-Xtrue) --> roughly(0, 1e-2)
    @fact norm(YY-Ytrue) --> roughly(0, 1e-2)
    @fact norm(ZZ-Ztrue) --> roughly(0, 1e-2)
    @fact ispsd(XX) --> true
    @fact ispsd(1/2*eye(3,3)-XX) --> true
    @fact ispsd(YY+ones(5,5)) --> true
    @fact ispsd(2*ones(5,5)-YY) --> true
    @fact ispsd(ones(4,4)-ZZ) --> true
    @fact trace(XX) --> roughly( 1, 1e-4)
    @fact trace(YY) --> roughly( 3, 1e-4)
    @fact trace(ZZ) --> roughly(-1, 1e-4)

    # Test SDP constraints
    m = Model(solver=solver)
    @defVar(m, X[1:3,1:3], SDP)

    @addSDPConstraint(m, ones(3,3) <= X)
    @setObjective(m, Min, trace(ones(3,3)*X))
    stat = solve(m)
    @fact stat --> :Optimal
    XX = getValue(X)
    @fact ispsd(XX) --> true
    @fact ispsd(XX - ones(3,3)) --> true
    @fact getObjectiveValue(m) --> roughly(9,1e-4)

    # Another test SDP
    m = Model(solver=solver)
    @defVar(m, X[1:3,1:3], SDP)
    @defVar(m, Y[1:2,1:2], SDP)

    C = eye(3,3)
    A1 = zeros(3,3)
    A1[1,1] = 1.0
    A2 = zeros(3,3)
    A2[2,2] = 1.0
    A3 = zeros(3,3)
    A3[3,3] = 1.0
    D = eye(2,2)
    B1 = ones(2,2)
    B2 = zeros(2,2)
    B2[1,1] = 1
    B3 = zeros(2,2)
    B3[1,2] = B3[2,1] = 2

    @setObjective(m, Min, trace(C*X)+1+trace(D*Y))
    @addConstraint(m, trace(A1*X-eye(3,3)/3) == 0)
    @addConstraint(m, 2*trace(A2*X) == 1)
    @addConstraint(m, trace(A3*X) >= 2)
    @addConstraint(m, trace(B1*Y) == 1)
    @addConstraint(m, trace(B2*Y) == 0)
    @addConstraint(m, trace(B3*Y) <= 0)
    @addConstraint(m, trace(A1*X)+trace(B1*Y) >= 1)
    @addConstraint(m, Y[2,2] == 1)

    stat = solve(m)
    @fact stat --> :Optimal
    XX, YY = getValue(X), getValue(Y)
    @fact trace(A1*XX-eye(3,3)/3) --> roughly(0, 1e-5)
    @fact 2*trace(A2*XX) --> roughly(1, 1e-5)
    @fact (trace(A3*XX) >= 2 - 1e-5) --> true
    @fact trace(B1*YY) --> roughly(1, 1e-5)
    @fact trace(B2*YY) --> roughly(0, 1e-5)
    @fact (trace(B3*YY) <= 1e-3) --> true
    @fact (trace(A1*XX)+trace(B1*YY) >= 1) --> true
    @fact YY[2,2] --> roughly(1, 1e-5)
    @fact norm(XX - diagm([1,.5,2])) --> roughly(0, 1e-3)
    @fact norm(YY - [0 0;0 1]) --> roughly(0, 1e-3)
end; end; end

facts("[sdp] Nonsensical SDPs") do
    m = Model()
    @fact_throws @defVar(m, unequal[1:5,1:6], SDP)
    # Some of these errors happen at compile time, so we can't use @fact_throws
    @fact macroexpand(:(@defVar(m, notone[1:5,2:6], SDP))).head --> :error
    @fact macroexpand(:(@defVar(m, oneD[1:5], SDP))).head --> :error
    @fact macroexpand(:(@defVar(m, threeD[1:5,1:5,1:5], SDP))).head --> :error
    @fact macroexpand(:(@defVar(m, psd[2] <= rand(2,2), SDP))).head --> :error
    @fact macroexpand(:(@defVar(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), SDP))).head --> :error
    @fact_throws @defVar(m, -ones(3,4) <= foo[1:4,1:4] <= ones(4,4), Symmetric)
    @fact_throws @defVar(m, -ones(4,4) <= foo[1:4,1:4] <= ones(4,5), Symmetric)
    @fact_throws @defVar(m, -rand(5,5) <= nonsymmetric[1:5,1:5] <= rand(5,5), Symmetric)
end

facts("[sdp] SDP with quadratics") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, X[1:2,1:2], SDP)
    @defVar(m, y[0:2])
    @addConstraint(m, y[0] >= 0)
    @addConstraint(m, y[1]^2 + y[2]^2 <= y[0]^2)
    @addSDPConstraint(m, X <= eye(2))
    @addConstraint(m, X[1,1] + X[1,2] == y[1] + y[2])
    @setObjective(m, Max, trace(X) - y[0])
    stat = solve(m)

    @fact stat --> :Optimal
    XX, yy = getValue(X), getValue(y)
    @fact ispsd(XX) --> true
    @fact (yy[0] >= 0) --> true
    @fact (yy[1]^2 + yy[2]^2 <= yy[0]^2 + 1e-4) --> true
    @fact ispsd(eye(2)-XX) --> true
    @fact ((XX[1,1] + XX[1,2]) - (yy[1] + yy[2])) --> roughly(0,1e-4)
    @fact norm(XX - eye(2)) --> roughly(0, 1e-4)
    @fact norm(yy[:] - [1/sqrt(2), 0.5, 0.5]) --> roughly(0, 1e-4)
    @fact getObjectiveValue(m) --> roughly(1.293, 1e-2)
end; end; end

# Adapt SDP atom tests from Convex.jl:
#   https://github.com/JuliaOpt/Convex.jl/blob/master/test/test_sdp.jl
# facts("[sdp] Test problem #1") do
# for solver in sdp_solvers
# context("With solver $(typeof(solver))") do
#     m = Model(solver=solver)
#     @defVar(m, X[1:2,1:2], SDP)
#     @setObjective(m, Max, X[1,1])
#     stat = solve(m)
#     @fact stat --> :Unbounded

#     setObjectiveSense(m, :Min)
#     stat = solve(m)
#     @fact stat --> :Optimal
#     @fact getObjectiveValue(m) --> roughly(0, 1e-4)

#     @addConstraint(m, X[1,1] == 1)
#     stat = solve(m)
#     @fact stat --> :Optimal
#     @fact getObjectiveValue(m) --> roughly(1, 1e-5)
#     @fact getValue(X[1,1]) --> roughly(1, 1e-5)

#     @setObjective(m, Min, sum(diag(X)))
#     stat = solve(m)
#     @fact stat --> :Optimal
#     @fact getObjectiveValue(m) --> roughly(1, 1e-5)
# end; end; end

facts("[sdp] Test problem #2") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, Y[1:3,1:3], SDP)
    @addConstraint(m, Y[2,1] <= 4)
    @addConstraint(m, Y[2,2] >= 3)
    @setObjective(m, Min, trace(Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getObjectiveValue(m) --> roughly(3, 1e-5)
end; end; end

facts("[sdp] Test problem #3") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, x >= 0)
    @defVar(m, Y[1:3,1:3], SDP)
    @addConstraint(m, Y[2,1] == 1)
    @setObjective(m, Min, Y[1,2])
    stat = solve(m)
    # @fact stat --> :Optimal
    @fact getObjectiveValue(m) --> roughly(1, 1e-4)
end; end; end

facts("[sdp] Test problem #4") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, x >= 0)
    @defVar(m, Y[1:3,1:3], SDP)
    @addConstraint(m, x >= 1)
    @addConstraint(m, Y[2,1] == 1)
    @setObjective(m, Min, x + Y[1,1])
    stat = solve(m)
    # @fact stat --> :Optimal # TODO: remove this once SCS starts behaving
    @fact getObjectiveValue(m) --> roughly(1, 1e-3)
end; end; end

function nuclear_norm(model, A)
    m, n = size(A,1), size(A,2)
    @defVar(model, U[1:m,1:m])
    @defVar(model, V[1:n,1:n])
    @addSDPConstraint(model, 0 ⪯ [U A; A' V])
    return 0.5(trace(U) + trace(V'))
end

facts("[sdp] Test problem #5") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, Y[1:3,1:3], SDP)
    @addConstraint(m, Y[2,1] <= 4)
    @addConstraint(m, Y[2,2] >= 3)
    @addConstraint(m, Y[3,3] <= 2)
    @setObjective(m, Min, nuclear_norm(m, Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getObjectiveValue(m) --> roughly(3, 1e-5)
end; end; end

function operator_norm(model, A)
    m, n = size(A,1), size(A,2)
    @defVar(model, t >= 0)
    @addSDPConstraint(model, [t*eye(n) A; A' eye(n)*t] >= 0)
    return t
end

facts("[sdp] Test problem #6") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, Y[1:3,1:3])
    @addConstraint(m, Y[2,1] <= 4)
    @addConstraint(m, Y[2,2] >= 3)
    @addConstraint(m, sum(Y) >= 12)
    @setObjective(m, Min, operator_norm(m, Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getObjectiveValue(m) --> roughly(4, 1e-5)
end; end; end

function lambda_max(model, A)
    m, n = size(A,1), size(A,2)
    @defVar(model, t)
    @addSDPConstraint(model, speye(n)*t - A ⪰ 0)
    @addSDPConstraint(model, A >= 0)
    return t
end

facts("[sdp] Test problem #7") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, Y[1:3,1:3])
    @addConstraint(m, Y[1,1] >= 4)
    @setObjective(m, Min, lambda_max(m, Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getObjectiveValue(m) --> roughly(4, 1e-5)
end; end; end

function lambda_min(model, A)
    m, n = size(A,1), size(A,2)
    @defVar(model, t)
    @addSDPConstraint(model, A - eye(n)*t >= 0)
    @addSDPConstraint(model, A >= 0)
    return t
end

facts("[sdp] Test problem #8") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, Y[1:3,1:3], SDP)
    @addConstraint(m, trace(Y) <= 6)
    @setObjective(m, Max, lambda_min(m, Y))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getObjectiveValue(m) --> roughly(2, 1e-5)
end; end; end

function matrix_frac(model, x, P)
    n = size(P,1)
    @defVar(model, t)
    @defVar(model, M[1:(n+1),1:(n+1)], SDP)
    @addConstraint(model, M[1:n,1:n] .== P)
    @addConstraint(model, M[1:n,n+1] .== x)
    @addConstraint(model, M[n+1,n+1] == t)
    return t
end

facts("[sdp] Test problem #9") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    n = 3
    x = [1,2,3]
    lb = 0.5eye(n)
    ub = 2eye(n)
    @defVar(m, lb[i,j] <= P[i=1:n,j=1:n] <= ub[i,j])
    @setObjective(m, Min, matrix_frac(m, x, P))
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getObjectiveValue(m) --> roughly(7, 1e-5)
end; end; end

facts("[sdp] Correlation example") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)

    @defVar(m, X[1:3,1:3], SDP)

    # Diagonal is 1s
    @addConstraint(m, X[1,1] == 1)
    @addConstraint(m, X[2,2] == 1)
    @addConstraint(m, X[3,3] == 1)

    # Bounds on the known correlations
    @addConstraint(m, X[1,2] >= -0.2)
    @addConstraint(m, X[1,2] <= -0.1)
    @addConstraint(m, X[2,3] >=  0.4)
    @addConstraint(m, X[2,3] <=  0.5)

    # Find upper bound
    @setObjective(m, Max, X[1,3])
    stat = solve(m)
    @fact stat --> :Optimal
    @fact (+0.8719 <= getValue(X)[1,3] <= +0.8720) --> true

    # Find lower bound
    @setObjective(m, Min, X[1,3])
    stat = solve(m)
    @fact stat --> :Optimal
    @fact (-0.9779 >= getValue(X)[1,3] >= -0.9799) --> true
end; end; end

facts("[sdp] Robust uncertainty example") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    include(joinpath("data","robust_uncertainty.jl"))
    R = 1
    d = 3
    𝛿 = 0.05
    ɛ = 0.05
    N = ceil((2+2log(2/𝛿))^2) + 1

    Γ1(𝛿,N) = (R/sqrt(N))*(2+sqrt(2*log(1/𝛿)))
    Γ2(𝛿,N) = (2R^2/sqrt(N))*(2+sqrt(2*log(2/𝛿)))

    for d in [3,5,8]; context("d = $d") do

        μhat = μhats[d]
        M = Ms[d]
        Σhat = 1/(d-1)*(M-ones(d)*μhat')'*(M-ones(d)*μhat')

        m = Model(solver=solver)

        @defVar(m, Σ[1:d,1:d], SDP)
        @defVar(m, u[1:d])
        @defVar(m, μ[1:d])

        @defVar(m, t1 >= 0)
        @defVar(m, L1[1:d])
        @addConstraint(m, L1 .== (μ-μhat))
        @addConstraint(m, sum{L1[i]^2, i=1:d} <= t1^2)
        @addConstraint(m, t1 <= Γ1(𝛿/2,N))

        @defVar(m, t2 >= 0)
        @defVar(m, L2[1:d,1:d])
        @addConstraint(m, L2 .== (Σ-Σhat))
        @addConstraint(m, sum{L2[i,j]^2, i=1:d, j=1:d} <= t2^2)
        @addConstraint(m, t2 <= Γ2(𝛿/2,N))

        A = [(1-ɛ)/ɛ (u-μ)';
             (u-μ)     Σ   ]
        @addSDPConstraint(m, A >= 0)

        c = cs[d]
        @setObjective(m, Max, dot(c,u))

        stat = solve(m)

        object = getObjectiveValue(m)
        exact = dot(μhat,c) + Γ1(𝛿/2,N)*norm(c) + sqrt((1-ɛ)/ɛ)*sqrt(dot(c,(Σhat+Γ2(𝛿/2,N)*eye(d,d))*c))
        @fact stat --> :Optimal
        @fact abs(object - exact) --> roughly(0, 1e-5)
    end; end
end; end; end

facts("[sdp] Robust uncertainty example (with norms)") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    include(joinpath("data","robust_uncertainty.jl"))
    R = 1
    d = 3
    𝛿 = 0.05
    ɛ = 0.05
    N = ceil((2+2log(2/𝛿))^2) + 1

    Γ1(𝛿,N) = (R/sqrt(N))*(2+sqrt(2*log(1/𝛿)))
    Γ2(𝛿,N) = (2R^2/sqrt(N))*(2+sqrt(2*log(2/𝛿)))

    for d in [3,5,8]; context("d = $d") do

        μhat = μhats[d]
        M = Ms[d]
        Σhat = 1/(d-1)*(M-ones(d)*μhat')'*(M-ones(d)*μhat')

        m = Model(solver=solver)

        @defVar(m, Σ[1:d,1:d], SDP)
        @defVar(m, u[1:d])
        @defVar(m, μ[1:d])

        @addConstraint(m, norm(μ-μhat) <= Γ1(𝛿/2,N))
        @addConstraint(m, vecnorm(Σ-Σhat) <= Γ2(𝛿/2,N))

        A = [(1-ɛ)/ɛ (u-μ)';
             (u-μ)     Σ   ]
        @addSDPConstraint(m, A >= 0)

        c = cs[d]
        @setObjective(m, Max, dot(c,u))

        stat = solve(m)

        object = getObjectiveValue(m)
        exact = dot(μhat,c) + Γ1(𝛿/2,N)*norm(c) + sqrt((1-ɛ)/ɛ)*sqrt(dot(c,(Σhat+Γ2(𝛿/2,N)*eye(d,d))*c))
        @fact stat --> :Optimal
        @fact abs(object - exact) --> roughly(0, 1e-5)
    end; end
end; end; end

facts("[sdp] Can't mix SDP and QP (#665)") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    model = Model(solver=solver)
    @defVar(model, Q[1:2, 1:2], SDP)
    @addConstraint(model, Q[1,1] - 1 == Q[2,2])
    @setObjective(model, Min, Q[1,1] * Q[1,1])
    @fact_throws solve(model)
end; end; end

facts("[sdp] Just another SDP") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
    model = Model(solver=solver)
    @defVar(model, Q[1:2, 1:2], SDP)
    @addConstraint(model, Q[1,1] - 1 == Q[2,2])
    @defVar(model, objective)
    T = [1 Q[1,1]; Q[1,1] objective]
    @addSDPConstraint(model, T ⪰ 0)
    @setObjective(model, Min, objective)
    @fact solve(model) --> :Optimal
    @fact getValue(Q) --> roughly([1 0;0 0], 1e-3)
    @fact getObjectiveValue(model) --> roughly(1, TOL)
end; end; end
