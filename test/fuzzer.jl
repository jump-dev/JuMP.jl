#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function random_aff_expr(N, vars::Vector{Symbol})
    ex = Expr(:call, :+)
    for _ in 1:N
        vl, vr = Any[], Any[]
        for i in randperm(length(vars))
            v = vars[i]
            if rand(Bool)
                push!(vl, :($(rand()-0.5)))
            end
            if rand(Bool)
                push!(vr, :($(rand()-0.5)))
            end
            if rand(Bool)
                push!(vl, :($(rand()-0.5) * $v))
            end
            if rand(Bool)
                push!(vr, :($(rand()-0.5) * $v))
            end
            if rand(Bool)
                push!(vl, :($v * $(rand()-0.5)))
            end
            if rand(Bool)
                push!(vl, :($v * $(rand()-0.5)))
            end
        end
        if rand(Bool)
            push!(vl, :($(rand()-0.5)))
        end
        if rand(Bool)
            push!(vr, :($(rand()-0.5)))
        end
        if !isempty(vl) || !isempty(vr)
            if isempty(vl)
                tmp = Expr(:call, :+, vr...)
            elseif isempty(vr)
                tmp = Expr(:call, :+, vl...)
            else
                tmp = Expr(:call, :*,
                        Expr(:call, :+, vl...),
                        Expr(:call, :+, vr...))
            end
            push!(ex.args, tmp)
        end
    end
    return ex
end

m = Model()
@defVar(m, x)
@defVar(m, y)
@defVar(m, z)
@defVar(m, w)
@defVar(m, v)

N = 5
vars = [:x, :y, :z, :w, :v]


const ε = 5eps()
nvars = length(vars)

function test_approx_equal_exprs(ex1, ex2)
    # test constant term
    abs(ex1.aff.constant - ex2.aff.constant) < ε || return false

    # test aff terms
    vals = zeros(nvars)
    for i in 1:length(ex1.aff.vars)
        vals[ex1.aff.vars[i].col] += ex1.aff.coeffs[i]
    end
    for i in 1:length(ex2.aff.vars)
        vals[ex2.aff.vars[i].col] -= ex2.aff.coeffs[i]
    end
    for v in vals
        abs(v) < ε || return false
    end

    # test quad terms
    qvals = zeros(nvars,nvars)
    for i in 1:length(ex1.qcoeffs)
        j,k = sort([ex1.qvars1[i].col, ex1.qvars2[i].col])
        qvals[j,k] += ex1.qcoeffs[i]
    end
    for i in 1:length(ex2.qcoeffs)
        j,k = sort([ex2.qvars1[i].col, ex2.qvars2[i].col])
        qvals[j,k] -= ex2.qcoeffs[i]
    end
    for v in vals
        abs(v) < ε || return false
    end
    return true
end

println("[fuzzer] Check macros for expression construction")

for _ in 1:100
    raff = random_aff_expr(N, vars)
    ex = @eval @defExpr($raff)
    test_approx_equal_exprs(ex, eval(raff)) ||
        error("The following expression did not pass the fuzzer:\n    $raff")
end
