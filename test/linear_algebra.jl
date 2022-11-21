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
    # TODO(odow): wrapper type that needs oneunit(VariableRef).
    # :UnitLowerTriangular,
    # :UnitUpperTriangular,
)
    f = getfield(LinearAlgebra, T)
    macro_matvec = Symbol("test_$(T)_macro_matvec")
    @eval begin
        function $macro_matvec()
            model = Model()
            @variable(model, x[1:2, 1:2])
            for X in (x, [1.0 2.0; 3.0 4.0])
                if X isa Matrix{VariableRef} && $f == LinearAlgebra.Hermitian
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
    if T != :Hermitian
        matvec = Symbol("test_$(T)_matvec")
        test_broken = T in (:LowerTriangular, :UpperTriangular)
        @eval begin
            function $matvec()
                model = Model()
                @variable(model, x[1:2, 1:2])
                if $test_broken
                    @test_broken $f(x) * [1.0, 2.0] isa Vector{AffExpr}
                else
                    @test $f(x) * [1.0, 2.0] isa Vector{AffExpr}
                end
                return
            end
        end
    end
end

end  # module

TestLinearAlgebra.runtests()
