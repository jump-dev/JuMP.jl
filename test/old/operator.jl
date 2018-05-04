#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/operator.jl
# Testing operator overloading is correct
#############################################################################
using JuMP
using Base.Test
using OffsetArrays

# To ensure the tests work on both Windows and Linux/OSX, we need
# to use the correct comparison operators in the strings
const  leq = JuMP.repl[:leq]
const  geq = JuMP.repl[:geq]
const   eq = JuMP.repl[:eq]
const Vert = JuMP.repl[:Vert]
const sub2 = JuMP.repl[:sub2]

# For "DimensionMismatch when performing vector-matrix multiplication with custom types #988"
import Base: +, *
struct MyType{T}
    a::T
end
struct MySumType{T}
    a::T
end
Base.one(::Type{MyType{T}}) where {T} = MyType(one(T))
Base.zero(::Type{MySumType{T}}) where {T} = MySumType(zero(T))
Base.zero(::MySumType{T}) where {T} = MySumType(zero(T))
Base.transpose(t::MyType) = MyType(t.a)
Base.transpose(t::MySumType) = MySumType(t.a)
+(t1::MyT, t2::MyS) where {MyT<:Union{MyType, MySumType}, MyS<:Union{MyType, MySumType}} = MySumType(t1.a+t2.a)
*(t1::MyType{S}, t2::T) where {S, T} = MyType(t1.a*t2)
*(t1::S, t2::MyType{T}) where {S, T} = MyType(t1*t2.a)
*(t1::MyType{S}, t2::MyType{T}) where {S, T} = MyType(t1.a*t2.a)


@testset "Operator overloads" begin

    _lt(x,y) = (x.col < y.col)
    function sort_expr!(x::AffExpr)
        idx = sortperm(x.vars, lt=_lt)
        x.vars = x.vars[idx]
        x.coeffs = x.coeffs[idx]
        return x
    end

    vec_eq(x,y) = vec_eq([x;], [y;])

    function vec_eq(x::AbstractArray, y::AbstractArray)
        size(x) == size(y) || return false
        for i in 1:length(x)
            v, w = convert(AffExpr,x[i]), convert(AffExpr,y[i])
            sort_expr!(v)
            sort_expr!(w)
            string(v) == string(w) || return false
        end
        return true
    end

    function vec_eq(x::Array{QuadExpr}, y::Array{QuadExpr})
        size(x) == size(y) || return false
        for i in 1:length(x)
            string(x[i]) == string(y[i]) || return false
        end
        return true
    end

    @testset "Higher-level operators" begin
    @testset "sum" begin
        sum_m = Model()
        @variable(sum_m, 0 ≤ matrix[1:3,1:3] ≤ 1, start = 1)
        # sum(j::JuMPArray{Variable})
        @test string(sum(matrix)) == "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"
        # sum(j::JuMPArray{Variable}) in a macro
        @objective(sum_m, Max, sum(matrix))
        @test string(sum_m.obj) == "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"

        # sum{T<:Real}(j::JuMPArray{T})
        @test isapprox(sum(getvalue(matrix)), 9)
        # sum(j::Array{Variable})
        @test string(sum(matrix[1:3,1:3])) == string(sum(matrix))
        # sum(affs::Array{AffExpr})
        @test string(sum([2*matrix[i,j] for i in 1:3, j in 1:3])) == "2 matrix[1,1] + 2 matrix[2,1] + 2 matrix[3,1] + 2 matrix[1,2] + 2 matrix[2,2] + 2 matrix[3,2] + 2 matrix[1,3] + 2 matrix[2,3] + 2 matrix[3,3]"

        S = [1,3]
        @variable(sum_m, x[S], start=1)
        # sum(j::JuMPDict{Variable})
        @test length(string(sum(x))) == 11 # order depends on hashing
        @test contains(string(sum(x)),"x[1]") == true
        @test contains(string(sum(x)),"x[3]") == true
        # sum{T<:Real}(j::JuMPDict{T})
        @test sum(getvalue(x)) == 2
    end

    @testset "dot" begin
        dot_m = Model()
        @variable(dot_m, 0 ≤ x[1:3] ≤ 1)
        c = vcat(1:3)
        @test string(dot(c,x)) == "x[1] + 2 x[2] + 3 x[3]"
        @test string(dot(x,c)) == "x[1] + 2 x[2] + 3 x[3]"

        A = [1 3 ; 2 4]
        @variable(dot_m, 1 ≤ y[1:2,1:2] ≤ 1)
        @test string(vecdot(A,y)) == "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"
        @test string(vecdot(y,A)) == "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"

        B = ones(2,2,2)
        @variable(dot_m, 0 ≤ z[1:2,1:2,1:2] ≤ 1)
        @test string(vecdot(B,z)) == "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"
        @test string(vecdot(z,B)) == "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"

        @objective(dot_m, Max, dot(x, ones(3)) - vecdot(y, ones(2,2)) )
        #solve(dot_m)
        for i in 1:3
            setvalue(x[i], 1)
        end
        for i in 1:2, j in 1:2
            setvalue(y[i,j], 1)
        end
        @test isapprox(dot(c, getvalue(x)), 6.0)
        @test isapprox(vecdot(A, getvalue(y)), 10.0)

        # https://github.com/JuliaOpt/JuMP.jl/issues/656
        issue656 = Model()
        @variable(issue656, x)
        floats = Float64[i for i in 1:2]
        anys   = Array{Any}(2)
        anys[1] = 10
        anys[2] = 20 + x
        @test dot(floats, anys) == 10 + 40 + 2x

        # https://github.com/JuliaOpt/JuMP.jl/pull/943
        pull943 = Model()
        @variable(pull943, x[1 : 10^6]);
        setvalue.(x, 1 : 10^6)
        @expression(pull943, testsum, sum(x[i] * i for i = 1 : 10^6))
        @expression(pull943, testdot1, dot(x, 1 : 10^6))
        @expression(pull943, testdot2, dot(1 : 10^6, x))
        @test isapprox(getvalue(testsum), getvalue(testdot1))
        @test isapprox(getvalue(testsum), getvalue(testdot2))
    end
    end

    @testset "Vectorized operations" begin

    @testset "Transpose" begin
        m = Model()
        @variable(m, x[1:3])
        @variable(m, y[1:2,1:3])
        @variable(m, z[2:5])
        @test vec_eq(x', [x[1] x[2] x[3]])
        @test vec_eq(transpose(x), [x[1] x[2] x[3]])
        @test vec_eq(y', [y[1,1] y[2,1]
                          y[1,2] y[2,2]
                          y[1,3] y[2,3]])
        @test vec_eq(transpose(y),
                    [y[1,1] y[2,1]
                     y[1,2] y[2,2]
                     y[1,3] y[2,3]])
        @test_throws ErrorException z'
        @test_throws ErrorException transpose(z)
    end

    @testset "Vectorized arithmetic" begin
        m = Model()
        @variable(m, x[1:3])
        A = [2 1 0
             1 2 1
             0 1 2]
        B = sparse(A)
        @variable(m, X11)
        @variable(m, X23)
        X = sparse([1, 2], [1, 3], [X11, X23], 3, 3) # for testing Variable
        @variable(m, Xd[1:3, 1:3])
        Y = sparse([1, 2], [1, 3], [2X11, 4X23], 3, 3) # for testing GenericAffExpr
        Yd = [2X11 0    0
              0    0 4X23
              0    0    0]
        Z = sparse([1, 2], [1, 3], [X11^2, 2X23^2], 3, 3) # for testing GenericQuadExpr
        Zd = [X11^2 0      0
              0     0 2X23^2
              0     0      0]
        v = [4, 5, 6]
        @test vec_eq(A*x, [2x[1] +  x[2]
                            2x[2] + x[1] + x[3]
                            x[2] + 2x[3]])
        @test vec_eq(A*x, B*x)
        @test vec_eq(A*x, @JuMP.Expression(B*x))
        @test vec_eq(@JuMP.Expression(A*x), @JuMP.Expression(B*x))
        @test vec_eq(x'*A, [2x[1]+x[2]; 2x[2]+x[1]+x[3]; x[2]+2x[3]]')
        @test vec_eq(x'*A, x'*B)
        @test vec_eq(x'*A, @JuMP.Expression(x'*B))
        @test vec_eq(@JuMP.Expression(x'*A), @JuMP.Expression(x'*B))
        @test vec_eq(x'*A*x, [2x[1]*x[1] + 2x[1]*x[2] + 2x[2]*x[2] + 2x[2]*x[3] + 2x[3]*x[3]])
        @test vec_eq(x'A*x, x'*B*x)
        @test vec_eq(x'*A*x, @JuMP.Expression(x'*B*x))
        @test vec_eq(@JuMP.Expression(x'*A*x), @JuMP.Expression(x'*B*x))

        y = A*x
        @test vec_eq(-x, [-x[1], -x[2], -x[3]])
        @test vec_eq(-y, [-2x[1] -  x[2]
                           -x[1] - 2x[2] -  x[3]
                                   -x[2] - 2x[3]])
        @test vec_eq(y + 1, [2x[1] +  x[2]         + 1
                              x[1] + 2x[2] +  x[3] + 1
                                      x[2] + 2x[3] + 1])
        @test vec_eq(y - 1, [2x[1] +  x[2]         - 1
                              x[1] + 2x[2] +  x[3] - 1
                                      x[2] + 2x[3] - 1])
        @test vec_eq(y + 2ones(3), [2x[1] +  x[2]         + 2
                                     x[1] + 2x[2] +  x[3] + 2
                                             x[2] + 2x[3] + 2])
        @test vec_eq(y - 2ones(3), [2x[1] +  x[2]         - 2
                                     x[1] + 2x[2] +  x[3] - 2
                                             x[2] + 2x[3] - 2])
        @test vec_eq(2ones(3) + y, [2x[1] +  x[2]         + 2
                                     x[1] + 2x[2] +  x[3] + 2
                                             x[2] + 2x[3] + 2])
        @test vec_eq(2ones(3) - y, [-2x[1] -  x[2]         + 2
                                     -x[1] - 2x[2] -  x[3] + 2
                                             -x[2] - 2x[3] + 2])
        @test vec_eq(y + x, [3x[1] +  x[2]
                              x[1] + 3x[2] +  x[3]
                                      x[2] + 3x[3]])
        @test vec_eq(x + y, [3x[1] +  x[2]
                              x[1] + 3x[2] +  x[3]
                                      x[2] + 3x[3]])
        @test vec_eq(2y + 2x, [6x[1] + 2x[2]
                               2x[1] + 6x[2] + 2x[3]
                                       2x[2] + 6x[3]])
        @test vec_eq(y - x, [ x[1] + x[2]
                              x[1] + x[2] + x[3]
                                     x[2] + x[3]])
        @test vec_eq(x - y, [-x[1] - x[2]
                             -x[1] - x[2] - x[3]
                                    -x[2] - x[3]])
        @test vec_eq(y + x[:], [3x[1] +  x[2]
                                 x[1] + 3x[2] +  x[3]
                                         x[2] + 3x[3]])
        @test vec_eq(x[:] + y, [3x[1] +  x[2]
                                 x[1] + 3x[2] +  x[3]
                                         x[2] + 3x[3]])

        @test vec_eq(@JuMP.Expression(A*x/2), A*x/2)
        @test vec_eq(X*v,  [4X11; 6X23; 0])
        @test vec_eq(v'*X,  [4X11  0   5X23])
        @test vec_eq(v.'*X, [4X11  0   5X23])
        @test vec_eq(X'*v,  [4X11;  0;  5X23])
        @test vec_eq(X.'*v, [4X11; 0;  5X23])
        @test vec_eq(X*A,  [2X11  X11  0
                            0     X23  2X23
                            0     0    0   ])
        @test vec_eq(A*X,  [2X11  0    X23
                            X11   0    2X23
                            0     0    X23])
        @test vec_eq(A*X', [2X11  0    0
                            X11   X23  0
                            0     2X23 0])
        @test vec_eq(X'*A, [2X11  X11  0
                            0     0    0
                            X23   2X23 X23])
        @test vec_eq(X.'*A, [2X11 X11  0
                             0    0    0
                             X23  2X23 X23])
        @test vec_eq(A'*X, [2X11  0 X23
                            X11   0 2X23
                            0     0 X23])
        @test vec_eq(X.'*A, X'*A)
        @test vec_eq(A.'*X, A'*X)
        @test vec_eq(X*A, X*B)
        @test vec_eq(Y'*A, Y.'*A)
        @test vec_eq(A*Y', A*Y.')
        @test vec_eq(Z'*A, Z.'*A)
        @test vec_eq(Xd'*Y, Xd.'*Y)
        @test vec_eq(Y'*Xd, Y.'*Xd)
        @test vec_eq(Xd'*Xd, Xd.'*Xd)
        # @test_broken vec_eq(A*X, B*X)
        # @test_broken vec_eq(A*X', B*X')
        @test vec_eq(X'*A, X'*B)
        # @test_broken(X'*X, X.'*X) # sparse quadratic known to be broken, see #912
    end

    @testset "Dot-ops" begin
        m = Model()
        @variable(m, x[1:2,1:2])
        A = [1 2;
             3 4]
        B = sparse(A)
        y = SparseMatrixCSC(2, 2, copy(B.colptr), copy(B.rowval), vec(x))
        @test vec_eq(A.+x, [1+x[1,1]  2+x[1,2];
                                       3+x[2,1]  4+x[2,2]])
        @test vec_eq(A.+x, B.+x)
        # @test vec_eq(A.+x, A.+y) == true
        # @test vec_eq(A.+y, B.+y) == true
        @test vec_eq(x.+A, [1+x[1,1]  2+x[1,2];
                                       3+x[2,1]  4+x[2,2]])
        @test vec_eq(x.+A, x.+B) == true
        @test vec_eq(x.+A, y.+A)
        @test vec_eq(x .+ x, [2x[1,1] 2x[1,2]; 2x[2,1] 2x[2,2]])
        # @test vec_eq(y.+A, y.+B) == true
        @test vec_eq(A.-x, [1-x[1,1]  2-x[1,2];
                                       3-x[2,1]  4-x[2,2]])
        @test vec_eq(A.-x, B.-x)
        @test vec_eq(A.-x, A.-y)
        @test vec_eq(x .- x, [zero(AffExpr) for _1 in 1:2, _2 in 1:2])
        # @test vec_eq(A.-y, B.-y) == true
        @test vec_eq(x.-A, [-1+x[1,1]  -2+x[1,2];
                                       -3+x[2,1]  -4+x[2,2]])
        @test vec_eq(x.-A, x.-B)
        @test vec_eq(x.-A, y.-A)
        # @test vec_eq(y.-A, y.-B) == true
        @test vec_eq(A.*x, [1*x[1,1]  2*x[1,2];
                                       3*x[2,1]  4*x[2,2]])
        @test vec_eq(A.*x, B.*x)
        @test vec_eq(A.*x, A.*y)
        # @test vec_eq(A.*y, B.*y) == true

        @test vec_eq(x.*A, [1*x[1,1]  2*x[1,2];
                                       3*x[2,1]  4*x[2,2]])
        @test vec_eq(x.*A, x.*B)
        @test vec_eq(x.*A, y.*A)
        # @test vec_eq(y.*A, y.*B) == true

        @test vec_eq(x .* x, [x[1,1]^2 x[1,2]^2; x[2,1]^2 x[2,2]^2])
        @test_throws ErrorException vec_eq(A./x, [1*x[1,1]  2*x[1,2];
                                              3*x[2,1]  4*x[2,2]])
        @test vec_eq(x./A, [1/1*x[1,1]  1/2*x[1,2];
                                       1/3*x[2,1]  1/4*x[2,2]])
        @test vec_eq(x./A, x./B)
        @test vec_eq(x./A, y./A)
        # @test vec_eq(A./y, B./y) == true

        @test vec_eq((2*x) / 3, full((2*y) / 3))
        @test vec_eq(2 * (x/3), full(2 * (y/3)))
        @test vec_eq(x[1,1] * A, full(x[1,1] * B))
    end

    @testset "Vectorized comparisons" begin
        m = Model()
        @variable(m, x[1:3])
        A = [1 2 3
             0 4 5
             6 0 7]
        B = sparse(A)
        if VERSION < v"0.6.0-dev.2074" # julia PR #19670
            @constraint(m, x'*A*x .>= 1)
        else
            # force vector output
            @constraint(m, reshape(x,(1,3))*A*x .>= 1)
        end
        @test vec_eq(m.quadconstr[1].terms, [x[1]*x[1] + 2x[1]*x[2] + 4x[2]*x[2] + 9x[1]*x[3] + 5x[2]*x[3] + 7x[3]*x[3] - 1])
        @test m.quadconstr[1].sense == :(>=)
        if VERSION < v"0.6.0-dev.2074" # julia PR #19670
            @constraint(m, x'*A*x .>= 1)
        else
            @constraint(m, x'*A*x >= 1)
        end
        @test vec_eq(m.quadconstr[1].terms, m.quadconstr[2].terms)

        mat = [ 3x[1] + 12x[3] +  4x[2]
                2x[1] + 12x[2] + 10x[3]
               15x[1] +  5x[2] + 21x[3]]

        @constraint(m, (x'A)' + 2A*x .<= 1)
        terms = map(v->v.terms, m.linconstr[1:3])
        lbs   = map(v->v.lb,    m.linconstr[1:3])
        ubs   = map(v->v.ub,    m.linconstr[1:3])
        @test vec_eq(terms, mat)
        @test lbs == fill(-Inf, 3)
        @test ubs == fill(   1, 3)
        @test vec_eq((x'A)' + 2A*x, (x'A)' + 2B*x)
        @test vec_eq((x'A)' + 2A*x, (x'B)' + 2A*x)
        @test vec_eq((x'A)' + 2A*x, (x'B)' + 2B*x)
        @test vec_eq((x'A)' + 2A*x, @JuMP.Expression((x'A)' + 2A*x))
        @test vec_eq((x'A)' + 2A*x, @JuMP.Expression((x'B)' + 2A*x))
        @test vec_eq((x'A)' + 2A*x, @JuMP.Expression((x'A)' + 2B*x))
        @test vec_eq((x'A)' + 2A*x, @JuMP.Expression((x'B)' + 2B*x))

        @constraint(m, -1 .<= (x'A)' + 2A*x .<= 1)
        terms = map(v->v.terms, m.linconstr[4:6])
        lbs   = map(v->v.lb,    m.linconstr[4:6])
        ubs   = map(v->v.ub,    m.linconstr[4:6])
        @test vec_eq(terms, mat) == true
        @test lbs == fill(-1, 3)
        @test ubs == fill( 1, 3)

        @constraint(m, -[1:3;] .<= (x'A)' + 2A*x .<= 1)
        terms = map(v->v.terms, m.linconstr[7:9])
        lbs   = map(v->v.lb,    m.linconstr[7:9])
        ubs   = map(v->v.ub,    m.linconstr[7:9])
        @test vec_eq(terms, mat) == true
        @test lbs == -[1:3;]
        @test ubs == fill( 1, 3)

        @constraint(m, -[1:3;] .<= (x'A)' + 2A*x .<= [3:-1:1;])
        terms = map(v->v.terms, m.linconstr[10:12])
        lbs   = map(v->v.lb,    m.linconstr[10:12])
        ubs   = map(v->v.ub,    m.linconstr[10:12])
        @test vec_eq(terms, mat) == true
        @test lbs == -[1:3;]
        @test ubs == [3:-1:1;]

        @constraint(m, -[1:3;] .<= (x'A)' + 2A*x .<= 3)
        terms = map(v->v.terms, m.linconstr[13:15])
        lbs   = map(v->v.lb,    m.linconstr[13:15])
        ubs   = map(v->v.ub,    m.linconstr[13:15])
        @test vec_eq(terms, mat) == true
        @test lbs == -[1:3;]
        @test ubs == fill(3,3)
    end

    end

    @testset "JuMPArray concatenation" begin
        m = Model()
        @variable(m, x[1:3])
        @variable(m, y[1:3,1:3])
        @variable(m, z[1:1])
        @variable(m, w[1:1,1:3])

        @test vec_eq([x y], [x[1] y[1,1] y[1,2] y[1,3]
                                        x[2] y[2,1] y[2,2] y[2,3]
                                        x[3] y[3,1] y[3,2] y[3,3]])

        @test vec_eq([x 2y+1], [x[1] 2y[1,1]+1 2y[1,2]+1 2y[1,3]+1
                                           x[2] 2y[2,1]+1 2y[2,2]+1 2y[2,3]+1
                                           x[3] 2y[3,1]+1 2y[3,2]+1 2y[3,3]+1])

        @test vec_eq([1 x'], [1 x[1] x[2] x[3]])
        @test vec_eq([2x;x], [2x[1],2x[2],2x[3],x[1],x[2],x[3]])
        # vcat on JuMPArray
        @test vec_eq([x;x], [x[1],x[2],x[3],x[1],x[2],x[3]])
        # hcat on JuMPArray
        @test vec_eq([x x], [x[1] x[1]
                                        x[2] x[2]
                                        x[3] x[3]])
        # hvcat on JuMPArray
        tmp1 = [z w; x y]
        tmp2 = [z[1] w[1,1] w[1,2] w[1,3]
                x[1] y[1,1] y[1,2] y[1,3]
                x[2] y[2,1] y[2,2] y[2,3]
                x[3] y[3,1] y[3,2] y[3,3]]
        @test vec_eq(tmp1, tmp2)
        tmp3 = [1 2x'
                x 2y-x*x']
        tmp4 = [1    2x[1]               2x[2]               2x[3]
                x[1] -x[1]*x[1]+2y[1,1]  -x[1]*x[2]+2y[1,2]  -x[1]*x[3] + 2y[1,3]
                x[2] -x[1]*x[2]+2y[2,1]  -x[2]*x[2]+2y[2,2]  -x[2]*x[3] + 2y[2,3]
                x[3] -x[1]*x[3]+2y[3,1]  -x[2]*x[3]+2y[3,2]  -x[3]*x[3] + 2y[3,3]]
        @test vec_eq(tmp3, tmp4)

        A = sprand(3, 3, 0.2)
        B = full(A)
        @test vec_eq([A y], [B y])
    end

    @testset "Operators for non-Array AbstractArrays" begin
        m = Model()
        @variable(m, x[1:3])

        # This is needed to compare arrays that have nonstandard indexing
        elements_equal(A::AbstractArray{T, N}, B::AbstractArray{T, N}) where {T, N} = all(a == b for (a, b) in zip(A, B))

        for x2 in (OffsetArray(x, -length(x)), view(x, :), sparse(x))
            @test elements_equal(+x, +x2)
            @test elements_equal(-x, -x2)
            @test elements_equal(x + first(x), x2 + first(x2))
            @test elements_equal(x - first(x), x2 - first(x2))
            @test elements_equal(first(x) - x, first(x2) - x2)
            @test elements_equal(first(x) + x, first(x2) + x2)
            @test elements_equal(2 * x, 2 * x2)
            @test elements_equal(first(x) + x2, first(x2) + x)
            @test sum(x) == sum(x2)
            if !JuMP.one_indexed(x2)
                @test_throws DimensionMismatch x + x2
            end
        end
    end

    @testset "Norm and diagm for non-Array AbstractArrays" begin
        m = Model()
        @variable(m, x[1:3])

        for x2 in (OffsetArray(x, -length(x)), view(x, :), sparse(x))
            if !JuMP.one_indexed(x2)
                @test_throws AssertionError diagm(x2)
            else
                @test diagm(x) == diagm(x2)
            end
            @test norm(x).terms == norm(x2).terms
        end
    end

    @testset "DimensionMismatch when performing vector-matrix multiplication with custom types #988" begin
        m = Model()
        @variable m Q[1:3, 1:3] SDP

        x = [MyType(1), MyType(2), MyType(3)]
        y = Q * x
        z = x' * Q
        ElemT = MySumType{JuMP.GenericAffExpr{Float64,JuMP.Variable}}
        @test typeof(y) == Vector{ElemT}
        @test size(y) == (3,)
        @test typeof(z) == (isdefined(Base, :RowVector) ? RowVector{ElemT, ConjArray{ElemT, 1, Vector{ElemT}}} : Matrix{ElemT})
        @test size(z) == (1, 3)
        for i in 1:3
            # Q is symmetric
            a = zero(JuMP.GenericAffExpr{Float64,JuMP.Variable})
            a += Q[1,i]
            a += 2Q[2,i]
            a += 3Q[3,i]
            # Q[1,i] + 2Q[2,i] + 3Q[3,i] is rearranged as 2 Q[2,3] + Q[1,3] + 3 Q[3,3]
            @test z[i].a == y[i].a == a
        end
    end
end
