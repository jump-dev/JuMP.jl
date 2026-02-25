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
what works, what doesn't, and make sure that we don't regress over time.
"""
module TestLinearAlgebra

using Test

using JuMP
import LinearAlgebra

function _test_matvec(model, A, b)
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
    matvec = Symbol("test_$(T)_matvec")
    matvec_alloc = Symbol("test_$(T)_matvec_alloc")
    @eval begin
        function $macro_matvec()
            model = Model()
            @variable(model, x[1:2, 1:2])
            for X in (x, [1.0 2.0; 3.0 4.0])
                if X isa Matrix{VariableRef} && $f == LinearAlgebra.Hermitian
                    continue  # Cannot call f(X) for this wrapper type.
                end
                _test_matvec(model, $f(X), [1.0, 2.0])
                _test_matvec(model, $f(X), x[:, 1])
                _test_matvec(model, $f(X), x[:, 1] .+ [1.0, 2.0])
            end
            return
        end
        function $matvec()
            if $f == LinearAlgebra.Hermitian
                return  # Cannot call f(X) for this wrapper type.
            end
            model = Model()
            @variable(model, x[1:2, 1:2])
            @test $f(x) * [1.0, 2.0] isa Vector{AffExpr}
            return
        end
        function $matvec_alloc()
            if $f == LinearAlgebra.Hermitian
                return
            end
            model = Model()
            @variable(model, x[1:2, 1:2])
            A, b = $f(x), [1.0, 2.0]
            alloc_macro = @allocated(@expression(model, A * b))
            # This is a bit of a weird one: we want to test that A * b in
            # the REPL is efficient. We can do that by checking that there
            # are the same number (or fewer) of allocations as the macro, in
            # the hope that the macro is efficient.
            @test @allocated(A * b) <= alloc_macro
            return
        end
    end
end

end  # module
