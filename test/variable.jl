#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/variable.jl
# Testing for VariableRef
#############################################################################
using JuMP
import JuMP.repl
using Compat
using Compat.Test

function variables_test(ModelType::Type{<:JuMP.AbstractModel}, VariableRefType::Type{<:JuMP.AbstractVariableRef})
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}

    @testset "constructors" begin
        # Constructors
        mcon = ModelType()

        @testset "No bound" begin
            @variable(mcon, nobounds)
            @test !JuMP.haslowerbound(nobounds)
            @test !JuMP.hasupperbound(nobounds)
            @test !JuMP.isfixed(nobounds)
            @test JuMP.name(nobounds) == "nobounds"

            @test zero(nobounds) isa AffExprType
            @test one(nobounds) isa AffExprType
        end

        @testset "Lower bound" begin
            @variable(mcon, lbonly >= 0, Bin)
            @test JuMP.haslowerbound(lbonly)
            @test JuMP.lowerbound(lbonly) == 0.0
            @test !JuMP.hasupperbound(lbonly)
            @test !JuMP.isfixed(lbonly)
            @test JuMP.isbinary(lbonly)
            @test !JuMP.isinteger(lbonly)
            @test isequal(mcon[:lbonly],lbonly)
            @test_throws ErrorException @variable(mcon, lbonly)
        end

        @testset "Upper bound" begin
            @variable(mcon, ubonly <= 1, Int)
            @test !JuMP.haslowerbound(ubonly)
            @test JuMP.hasupperbound(ubonly)
            @test JuMP.upperbound(ubonly) == 1.0
            @test !JuMP.isfixed(ubonly)
            @test !JuMP.isbinary(ubonly)
            @test JuMP.isinteger(ubonly)
            @test isequal(mcon[:ubonly],ubonly)
        end

        @testset "Interval" begin
            function has_bounds(var, lb, ub)
                @test JuMP.haslowerbound(var)
                @test JuMP.lowerbound(var) == lb
                @test JuMP.hasupperbound(var)
                @test JuMP.upperbound(var) == ub
                @test !JuMP.isfixed(var)
            end

            @variable(mcon, 0 <= bothb1 <= 1)
            has_bounds(bothb1, 0.0, 1.0)
            @variable(mcon, 0 ≤  bothb2 ≤  1)
            has_bounds(bothb2, 0.0, 1.0)
            @variable(mcon, 1 >= bothb3 >= 0)
            has_bounds(bothb3, 0.0, 1.0)
            @variable(mcon, 1 ≥  bothb4 ≥  0)
            has_bounds(bothb4, 0.0, 1.0)
            @test_macro_throws ErrorException @variable(mcon, 1 ≥ bothb5 ≤ 0)
            @test_macro_throws ErrorException @variable(mcon, 1 ≤ bothb6 ≥ 0)
        end

        @testset "Fix" begin
            @variable(mcon, fixed == 1.0)
            @test !JuMP.haslowerbound(fixed)
            @test !JuMP.hasupperbound(fixed)
            @test JuMP.isfixed(fixed)
            @test JuMP.fixvalue(fixed) == 1.0
        end

        @testset "Custom index sets" begin
            @variable(mcon, onerangeub[-7:1] <= 10, Int)
            @variable(mcon, manyrangelb[0:1,10:20,1:1] >= 2)
            @test JuMP.haslowerbound(manyrangelb[0,15,1])
            @test JuMP.lowerbound(manyrangelb[0,15,1]) == 2
            @test !JuMP.hasupperbound(manyrangelb[0,15,1])

            s = ["Green","Blue"]
            @variable(mcon, x[i=-10:10,s] <= 5.5, Int, start=i+1)
            @test JuMP.upperbound(x[-4,"Green"]) == 5.5
            @test JuMP.name(x[-10,"Green"]) == "x[-10,Green]"
            # TODO: broken because of https://github.com/JuliaOpt/MathOptInterface.jl/issues/302
            #@test JuMP.startvalue(x[-3,"Blue"]) == -2
            @test isequal(mcon[:onerangeub][-7],onerangeub[-7])
            @test_throws KeyError mcon[:foo]
        end

        @testset "Anonymous variable" begin
            @test_throws ErrorException @variable(mcon, [(0,0)]) # #922
            x = @variable(mcon, [(0,2)])
            @test JuMP.name(x[0]) == ""
            @test JuMP.name(x[2]) == ""
        end
    end

    @testset "isvalid and delete variable" begin
        m = ModelType()
        @variable(m, x)
        @test MOI.isvalid(m, x)
        MOI.delete!(m, x)
        @test !MOI.isvalid(m, x)
    end


    @testset "get and set bounds" begin
        m = ModelType()
        @variable(m, 0 <= x <= 2)
        @test JuMP.lowerbound(x) == 0
        @test JuMP.upperbound(x) == 2
        setlowerbound(x, 1)
        @test JuMP.lowerbound(x) == 1
        setupperbound(x, 3)
        @test JuMP.upperbound(x) == 3
        @variable(m, q, Bin)
        @test !JuMP.haslowerbound(q)
        @test !JuMP.hasupperbound(q)

        @variable(m, 0 <= y <= 1, Bin)
        @test JuMP.lowerbound(y) == 0
        @test JuMP.upperbound(y) == 1

        @variable(m, fixedvar == 2)
        @test JuMP.fixvalue(fixedvar) == 2.0
        JuMP.fix(fixedvar, 5)
        @test JuMP.fixvalue(fixedvar) == 5
        @test_throws AssertionError JuMP.lowerbound(fixedvar)
        @test_throws AssertionError JuMP.upperbound(fixedvar)
    end

    @testset "get and set start" begin
        m = ModelType()
        @variable(m, x[1:3])
        x0 = collect(1:3)
        JuMP.setstartvalue.(x, x0)
        @test JuMP.startvalue.(x) == x0
        @test JuMP.startvalue.([x[1],x[2],x[3]]) == x0

        @variable(m, y[1:3,1:2])
        @test_throws DimensionMismatch JuMP.setstartvalue.(y, collect(1:6))
    end

    @testset "get and set integer/binary" begin
        m = ModelType()
        @variable(m, x[1:3])

        JuMP.setinteger(x[2])
        @test JuMP.isinteger(x[2])
        JuMP.unsetinteger(x[2])
        @test !JuMP.isinteger(x[2])

        JuMP.setbinary(x[1])
        @test JuMP.isbinary(x[1])
        @test_throws AssertionError JuMP.setinteger(x[1])
        JuMP.unsetbinary(x[1])
        @test !JuMP.isbinary(x[1])
        # TODO test binary/integer keyword arguments
    end

    @testset "repeated elements in index set (issue #199)" begin
        repeatmod = ModelType()
        s = [:x,:x,:y]
        @test_throws ErrorException @variable(repeatmod, x[s], container = JuMPArray)
        @test_throws ErrorException @variable(repeatmod, x[s], container = Dict)
        @test_throws ErrorException @variable(repeatmod, x[s,[1]], container = JuMPArray)
        @test_throws ErrorException @variable(repeatmod, x[s,[1]], container = Dict)
    end

    @testset "Base.OneTo as index set (#933)" begin
        m = ModelType()
        x = @variable(m, [Base.OneTo(3), 1:2], container=Auto)
        @test x isa Matrix{VariableRefType}
        @test size(x) == (3,2)
        x = @variable(m, [Base.OneTo(3), 1:2], container=Array)
        @test x isa Matrix{VariableRefType}
        @test size(x) == (3,2)
        x = @variable(m, [Base.OneTo(3), 1:2], container=JuMPArray)
        @test x isa JuMPArray{VariableRefType}
        @test length.(indices(x)) == (3,2)
    end

    @testset "basename= in @variable" begin
        m = ModelType()
        @variable(m, x)
        @test JuMP.name(x) == "x"

        y = @variable(m, basename="foo")
        @test JuMP.name(y) == "foo"

        z = @variable(m, z[i=2:3], basename="t")
        @test JuMP.name(z[2]) == "t[2]"
        @test JuMP.name(z[3]) == "t[3]"
    end

    # TODO reenable when printing comes back
    # @testset "condition in indexing" begin
    #    fa = repl[:for_all]
    #    inset, dots = repl[:in], repl[:dots]
    #    condmod = ModelType()
    #    @variable(condmod, x[i=1:10; iseven(i)])
    #    @variable(condmod, y[j=1:10,k=3:2:9; isodd(j+k) && k <= 8])
    #    @test length(x.tupledict) == 5
    #    @test length(y.tupledict) == 15
    #    @test string(condmod) == "Min 0\nSubject to\n x[i] $fa i $inset {1,2,$dots,9,10} s.t. iseven(i)\n y[j,k] $fa j $inset {1,2,$dots,9,10}, k $inset {3,5,7,9} s.t. isodd(j + k) and k <= 8\n"
    # end

    @testset "@variable returning Array{VariableRef}" begin
        m = ModelType()
        @variable(m, x[1:3,1:4,1:2], start = 0.0)
        @variable(m, y[1:0], start = 0.0)
        @variable(m, z[1:4], start = 0.0)

        @test typeof(x) == Array{VariableRefType,3}
        @test typeof(y) == Vector{VariableRefType}
        @test typeof(z) == Vector{VariableRefType}

        @test typeof(JuMP.startvalue.(x)) == Array{Float64,3}
        # No type to infer for an empty collection.
        @test typeof(JuMP.startvalue.(y)) == Array{Any,1}
        @test typeof(JuMP.startvalue.(z)) == Array{Float64,1}
    end

    @testset "startvalue on empty things" begin
        m = ModelType()
        @variable(m, x[1:4,  1:0,1:3], start = 0) # Array{VariableRef}
        @variable(m, y[1:4,  2:1,1:3], start = 0) # JuMPArray
        @variable(m, z[1:4,Set(),1:3], start = 0) # Dict

        @test JuMP.startvalue.(x) == Array{Float64}(4, 0, 3)
        # TODO: Decide what to do here. I don't know if we still need to test this given broadcast syntax.
        # @test typeof(JuMP.startvalue(y)) <: JuMP.JuMPArray{Float64}
        # @test JuMP.size(JuMP.startvalue(y)) == (4,0,3)
        # @test typeof(JuMP.startvalue(z)) == JuMP.JuMPArray{Float64,3,Tuple{UnitRange{Int},Set{Any},UnitRange{Int}}}
        # @test length(JuMP.startvalue(z)) == 0
    end

# Slices three-dimensional JuMPArray x[I,J,K]
# I,J,K can be singletons, ranges, colons, etc.
function sliceof(x, I, J, K)
    y = Array{VariableRefType}(length(I), length(J), length(K))

    ii = 1
    jj = 1
    kk = 1
    for i in I
        for j in J
            for k in K
                y[ii,jj,kk] = x[i,j,k]
                kk += 1
            end
            jj += 1
            kk = 1
        end
        ii += 1
        jj = 1
    end
    idx = [length(I)==1, length(J)==1, length(K)==1]
    squeeze(y, tuple(find(idx)...))
end

    @testset "Slices of JuMPArray (#684)" begin
        m = ModelType()
        @variable(m, x[1:3, 1:4,1:2], container=JuMPArray)
        @variable(m, y[1:3,-1:2,3:4])
        @variable(m, z[1:3,-1:2:4,3:4])
        @variable(m, w[1:3,-1:2,[:red,"blue"]])

        #@test x[:] == vec(sliceof(x, 1:3, 1:4, 1:2))
        @test x isa JuMPArray
        @test x[:,:,:].data == sliceof(x, 1:3, 1:4, 1:2)
        @test x[1,:,:].data == sliceof(x, 1, 1:4, 1:2)
        @test x[1,:,2].data == sliceof(x, 1, 1:4, 2)
        @test_throws KeyError x[1,:,3]
        #@test x[1:2,:,:].data == sliceof(x, 1:2, 1:4, 1:2)
        #@test x[1:2,:,2].data == sliceof(x, 1:2, 1:4, 2)
        #@test x[1:2,:,1:2].data == sliceof(x, 1:2, 1:4, 1:2)
        @test_throws KeyError x[1:2,:,1:3]

        #@test y[:] == vec(sliceof(y, 1:3, -1:2, 3:4))
        @test y[:,:,:].data == sliceof(y, 1:3, -1:2, 3:4)
        @test y[1,:,:].data == sliceof(y, 1, -1:2, 3:4)
        @test y[1,:,4].data == sliceof(y, 1, -1:2, 4)
        @test_throws KeyError y[1,:,5]
        # @test y[1:2,:,:] == sliceof(y, 1:2, -1:2, 3:4)
        # @test y[1:2,:,4] == sliceof(y, 1:2, -1:2, 4)
        # @test y[1:2,:,3:4] == sliceof(y, 1:2, -1:2, 3:4)
        # @test_throws BoundsError y[1:2,:,1:3]

        #@test z[:] == vec(sliceof(z, 1:3, -1:2:4, 3:4))
        @test z[:,1,:].data == sliceof(z, 1:3, 1, 3:4)
        @test z[1,1,:].data == sliceof(z, 1, 1, 3:4)
        @test_throws KeyError z[:,5,3]
        # @test z[1:2,1,:] == sliceof(z, 1:2, 1, 3:4)
        # @test z[1:2,1,4] == sliceof(z, 1:2, 1, 4)
        # @test z[1:2,1,3:4] == sliceof(z, 1:2, 1, 3:4)
        # @test_throws BoundsError z[1:2,1,1:3]

        #@test w[:] == vec(sliceof(w, 1:3, -1:2, [:red,"blue"]))
        @test w[:,:,:] == w
        @test w[1,:,"blue"].data == sliceof(w, 1, -1:2, ["blue"])
        @test w[1,:,:red].data == sliceof(w, 1, -1:2, [:red])
        @test_throws KeyError w[1,:,"green"]
        # @test w[1:2,:,"blue"] == sliceof(w, 1:2, -1:2, ["blue"])
        # @test_throws ErrorException w[1:2,:,[:red,"blue"]]
    end

    @testset "Can't use end for indexing a JuMPArray or Dict" begin
        m = ModelType()
        @variable(m, x[0:2,1:4])
        @variable(m, y[i=1:4,j=1:4;true])
        @variable(m, z[0:2])
        @test_throws ErrorException x[end,1]
        @test_throws ErrorException x[end-1]
        @test_throws ErrorException x[0,end-1]
        @test_throws MethodError y[end,end-1]
        @test_throws MethodError y[end,1]
        @test_throws ErrorException z[end]
    end

    @testset "Unsigned dimension lengths" begin
        m = ModelType()
        t = UInt(4)
        @variable(m, x[1:t])
        #@constraintref(y[1:t])
        @test JuMP.num_variables(m) == 4
    end

    # TODO decide what to do here
    # @testset "getstart on sparse array (#889)" begin
    #     m = ModelType()
    #     @variable(m, x)
    #     @variable(m, y)
    #     X = sparse([1, 3], [2, 3], [x, y])
    #
    #     @test typeof(X) == SparseMatrixCSC{VariableRefType, Int}
    #
    #     setstart(x, 1)
    #     setstart(y, 2)
    #
    #     Y = getstart.(X)
    #     @test typeof(Y) == SparseMatrixCSC{Float64, Int}
    #     @test Y == sparse([1, 3], [2, 3], [1, 2])
    # end

    @testset "Symmetric variable" begin
        m = ModelType()

        @variable m x[1:2, 1:2] Symmetric
        @test x isa Symmetric
        @test x[1, 2] === x[2, 1]

        y = @variable m [1:2, 1:2] Symmetric
        @test y isa Symmetric
        @test y[1, 2] === y[2, 1]
    end
end

@testset "Variables for JuMP.Model" begin
    variables_test(Model, VariableRef{Model{JuMP.NonDirectBackendType}})
end

@testset "Variables for JuMPExtension.MyModel" begin
    variables_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end
