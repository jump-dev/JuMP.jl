#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# maxcut_sdp.jl
#
# Solves the SDP relaxation of the classic MAXCUT problem:
# max   L•X
# s.t.  diag(X) == e
#       X ≽ 0
# where
#  L = 1/4(Diag(W e) - W)
#  W = edge-weight matrix
#  e = vector of 1s
#
# then applies the Goemans-Williamson algorithm
#############################################################################
using JuMP
using SCS

solver = SCSSolver(eps=1e-6)

function solve_maxcut_sdp(n, W)
    L = 0.25 * (diagm(W*ones(n)) - W)

    # Solve the SDP relaxation
    m = Model(solver=solver)
    @variable(m, X[1:n,1:n], SDP)
    @objective(m, Max, vecdot(L,X))
    @constraint(m, diag(X) .== 1)
    solve(m)

    # Cholesky the result
    F = cholfact(getvalue(X)[:,:], :U, Val{true})
    V = (F[:P]*F[:L])'

    # Normalize columns
    for i = 1:n
        V[:,i] ./= norm(V[:,i])
    end

    # Generate "random" vector
    # - seeded on problem size for repeatability
    # - for all the problems in this file, the
    #   solutions are "integral" anyway so there
    #   isn't really a need for this
    r = rand(n)
    cut = ones(n)
    for i = 1:n
        if sum(r' * V[:,i]) <= 0
            cut[i] = -1
        end
    end

    return cut, sum(L.*(cut*cut'))
end

function test0()
    #   [1] --- 5 --- [2]
    #
    # Solution:
    #  S  = {1}
    #  S' = {2}
    n = 2
    W = [0.0 5.0;
         5.0 0.0]
    cut, cutval = solve_maxcut_sdp(n, W)

    @assert cut[1] != cut[2]

    println("Solution for Graph 0 = $cutval")
    println(cut)
end

function test1()
    #   [1] --- 5 --- [2]
    #    |  \          |
    #    |    \        |
    #    7      6      1
    #    |        \    |
    #    |          \  |
    #   [3] --- 1 --- [4]
    #
    # Solution:
    #  S  = {1}
    #  S' = {2,3,4}
    n = 4
    W = [0.0 5.0 7.0 6.0;
         5.0 0.0 0.0 1.0;
         7.0 0.0 0.0 1.0;
         6.0 1.0 1.0 0.0]
    cut, cutval = solve_maxcut_sdp(n, W)
    @assert (v = cut[2]) == cut[3] == cut[4]
    @assert cut[1] != v

    println("Solution for Graph 1 = $cutval")
    println(cut)
end

function test2()
    #   [1] --- 1 --- [2]
    #    |             |
    #    |             |
    #    9             9
    #    |             |
    #    |             |
    #   [3] --- 1 --- [4]
    #
    # Solution:
    #  S  = {2,3}
    #  S' = {1,4}
    n = 4
    W = [0.0 1.0 9.0 0.0;
         1.0 0.0 0.0 9.0;
         9.0 0.0 0.0 1.0;
         0.0 9.0 1.0 0.0]
    cut, cutval = solve_maxcut_sdp(n, W)
    @assert (v = cut[1]) == cut[4]
    @assert (w = cut[2]) == cut[3]
    @assert v != w

    println("Solution for Graph 2 = $cutval")
    println(cut)
end

test0()
test1()
test2()
