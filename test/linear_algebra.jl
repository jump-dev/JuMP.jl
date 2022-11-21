#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    TestLinearAlgebra

The purpose of this file is to test various fallbacks that MutableArithmetics
implements in order to fix compatibility between JuMP types and the various
wrappers in LinearAlgebra.

Historically, these have been a big source of issues for users, so let's track
what works, what doesnt, and make sure that we don't regress over time.
"""
module TestLinearAlgebra

using Test

using JuMP
import LinearAlgebra

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

for T in (
    :Adjoint,
    :Diagonal,
    :Hermitian,
    :Symmetric,
    :Transpose,
    :LowerTriangular,
    :UpperTriangular,
    # :UnitLowerTriangular,
    # :UnitUpperTriangular,
)
    f = getfield(LinearAlgebra, T)
    macro_matvec = Symbol("test_$(T)_macro_matvec")
    @eval begin
        function $macro_matvec()
            model = Model()
            @variable(model, x[1:2, 1:2])
            for X in (x, [1.0 2.0; 3.0 1.0])
                skips = (
                    LinearAlgebra.Hermitian,
                    LinearAlgebra.UnitUpperTriangular,
                    LinearAlgebra.UnitUpperTriangular,
                )
                if X isa Matrix{VariableRef} && $f in skips
                    # Cannot call f(X) for this wrapper type.
                    continue
                end
                A = $f(X)
                for b in ([1.0, 2.0], x[:, 1], x[:, 1] .+ [1.0, 2.0])
                    expr_vec = @expression(model, A * b)
                    expr_adj = @expression(model, b' * A)
                    @test size(expr_vec) == (2,)
                    @test size(expr_adj) == (1, 2)
                    if expr_vec isa AbstractArray{Float64}
                        @test A * b ≈ expr_vec
                        @test b' * A ≈ expr_adj
                    else
                        for i in 1:2
                            @test isequal_canonical(
                                expr_vec[i],
                                sum(A[i, j] * b[j] for j in 1:2),
                            )
                            @test isequal_canonical(
                                expr_adj[i],
                                sum(b[j] * A[j, i] for j in 1:2),
                            )
                        end
                    end
                end
            end
            return
        end
    end
end

function test_matmul_lower_triangular()
    model = Model()
    @variable(model, x[1:2, 1:2])
    b = [1.0, 2.0]
    @test_broken LinearAlgebra.LowerTriangular(x) * b isa Vector{AffExpr}
    return
end

function test_matmul_upper_triangular()
    model = Model()
    @variable(model, x[1:2, 1:2])
    b = [1.0, 2.0]
    @test_broken LinearAlgebra.UpperTriangular(x) * b isa Vector{AffExpr}
    return
end

end  # module

TestLinearAlgebra.runtests()
