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
using Compat
using Compat.Test

function variables_test(ModelType::Type{<:JuMP.AbstractModel}, VariableRefType::Type{<:JuMP.AbstractVariableRef})
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}

    @testset "Constructors" begin
        mcon = ModelType()

        @testset "No bound" begin
            @variable(mcon, nobounds)
            @test !JuMP.has_lower_bound(nobounds)
            @test !JuMP.has_upper_bound(nobounds)
            @test !JuMP.is_fixed(nobounds)
            @test JuMP.name(nobounds) == "nobounds"

            @test zero(nobounds) isa AffExprType
            @test one(nobounds) isa AffExprType
        end

        @testset "Lower bound" begin
            @variable(mcon, lbonly >= 0, Bin)
            @test JuMP.has_lower_bound(lbonly)
            @test JuMP.lower_bound(lbonly) == 0.0
            @test !JuMP.has_upper_bound(lbonly)
            @test !JuMP.is_fixed(lbonly)
            @test JuMP.is_binary(lbonly)
            @test !JuMP.is_integer(lbonly)
            @test isequal(mcon[:lbonly], lbonly)
            JuMP.delete_lower_bound(lbonly)
            @test !JuMP.has_lower_bound(lbonly)
            # Name already used
            @test_throws ErrorException @variable(mcon, lbonly)
        end

        @testset "Lower bound with constant on lhs" begin
            @variable(mcon, 0 <= lblhs, Bin)
            @test JuMP.has_lower_bound(lblhs)
            @test JuMP.lower_bound(lblhs) == 0.0
            @test !JuMP.has_upper_bound(lblhs)
            @test !JuMP.is_fixed(lblhs)
            @test JuMP.is_binary(lblhs)
            @test !JuMP.is_integer(lblhs)
            @test isequal(mcon[:lblhs],lblhs)
        end

        @testset "Upper bound" begin
            @variable(mcon, ubonly <= 1, Int)
            @test !JuMP.has_lower_bound(ubonly)
            @test JuMP.has_upper_bound(ubonly)
            @test JuMP.upper_bound(ubonly) == 1.0
            @test !JuMP.is_fixed(ubonly)
            @test !JuMP.is_binary(ubonly)
            @test JuMP.is_integer(ubonly)
            @test isequal(mcon[:ubonly], ubonly)
            JuMP.delete_upper_bound(ubonly)
            @test !JuMP.has_upper_bound(ubonly)
        end

        @testset "Upper bound" begin
            @variable(mcon, 1 >= ublhs, Int)
            @test !JuMP.has_lower_bound(ublhs)
            @test JuMP.has_upper_bound(ublhs)
            @test JuMP.upper_bound(ublhs) == 1.0
            @test !JuMP.is_fixed(ublhs)
            @test !JuMP.is_binary(ublhs)
            @test JuMP.is_integer(ublhs)
            @test isequal(mcon[:ublhs],ublhs)
        end

        @testset "Interval" begin
            function has_bounds(var, lb, ub)
                @test JuMP.has_lower_bound(var)
                @test JuMP.lower_bound(var) == lb
                @test JuMP.has_upper_bound(var)
                @test JuMP.upper_bound(var) == ub
                @test !JuMP.is_fixed(var)
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
            @test !JuMP.has_lower_bound(fixed)
            @test !JuMP.has_upper_bound(fixed)
            @test JuMP.is_fixed(fixed)
            @test JuMP.fix_value(fixed) == 1.0
            JuMP.unfix(fixed)
            @test !JuMP.is_fixed(fixed)
        end

        @testset "Custom index sets" begin
            @variable(mcon, onerangeub[-7:1] <= 10, Int)
            @variable(mcon, manyrangelb[0:1,10:20,1:1] >= 2)
            @test JuMP.has_lower_bound(manyrangelb[0,15,1])
            @test JuMP.lower_bound(manyrangelb[0,15,1]) == 2
            @test !JuMP.has_upper_bound(manyrangelb[0,15,1])

            s = ["Green","Blue"]
            @variable(mcon, x[i=-10:10,s] <= 5.5, Int, start=i+1)
            @test JuMP.upper_bound(x[-4,"Green"]) == 5.5
            @test JuMP.name(x[-10,"Green"]) == "x[-10,Green]"
            # TODO: broken because of https://github.com/JuliaOpt/MathOptInterface.jl/issues/302
            #@test JuMP.start_value(x[-3,"Blue"]) == -2
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
        @test JuMP.is_valid(m, x)
        JuMP.delete(m, x)
        @test !JuMP.is_valid(m, x)
    end


    @testset "get and set bounds" begin
        m = ModelType()
        @variable(m, 0 <= x <= 2)
        @test JuMP.lower_bound(x) == 0
        @test JuMP.upper_bound(x) == 2
        set_lower_bound(x, 1)
        @test JuMP.lower_bound(x) == 1
        set_upper_bound(x, 3)
        @test JuMP.upper_bound(x) == 3
        @variable(m, q, Bin)
        @test !JuMP.has_lower_bound(q)
        @test !JuMP.has_upper_bound(q)

        @variable(m, 0 <= y <= 1, Bin)
        @test JuMP.lower_bound(y) == 0
        @test JuMP.upper_bound(y) == 1

        @variable(m, fixedvar == 2)
        @test JuMP.fix_value(fixedvar) == 2.0
        JuMP.fix(fixedvar, 5)
        @test JuMP.fix_value(fixedvar) == 5
        @test_throws AssertionError JuMP.lower_bound(fixedvar)
        @test_throws AssertionError JuMP.upper_bound(fixedvar)
    end

    @testset "get and set start" begin
        m = ModelType()
        @variable(m, x[1:3])
        x0 = collect(1:3)
        JuMP.set_start_value.(x, x0)
        @test JuMP.start_value.(x) == x0
        @test JuMP.start_value.([x[1],x[2],x[3]]) == x0

        @variable(m, y[1:3,1:2])
        @test_throws DimensionMismatch JuMP.set_start_value.(y, collect(1:6))
    end

    @testset "get and set integer/binary" begin
        m = ModelType()
        @variable(m, x[1:3])

        JuMP.set_integer(x[2])
        @test JuMP.is_integer(x[2])
        JuMP.unset_integer(x[2])
        @test !JuMP.is_integer(x[2])

        JuMP.set_binary(x[1])
        @test JuMP.is_binary(x[1])
        @test_throws Exception JuMP.set_integer(x[1])
        JuMP.unset_binary(x[1])
        @test !JuMP.is_binary(x[1])

        @variable(m, y, binary = true)
        @test JuMP.is_binary(y)
        @test_throws Exception JuMP.set_integer(y)
        JuMP.unset_binary(y)
        @test !JuMP.is_binary(y)

        @variable(m, z, integer = true)
        @test JuMP.is_integer(z)
        @test_throws Exception JuMP.set_binary(z)
        JuMP.unset_integer(z)
        @test !JuMP.is_integer(z)
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
        @test length.(Compat.axes(x)) == (3,2)
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

    @testset "condition in indexing" begin
        function test_one_dim(x)
            @test length(x) == 5
            for i in 1:10
                if iseven(i)
                    @test haskey(x, i)
                else
                    @test !haskey(x, i)
                end
            end
        end

        function test_two_dim(y)
            @test length(y) == 15
            for j in 1:10, k in 3:2:9
                if isodd(j+k) && k <= 8
                    @test haskey(y, (j,k))
                else
                    @test !haskey(y, (j,k))
                end
            end
        end


        model = ModelType()
        # Parses as ref on 0.7.
        @variable(model, named_one_dim[i = 1:10; iseven(i)])
        # Parses as vcat on 0.7.
        anon_one_dim = @variable(model, [i = 1:10; iseven(i)])
        # Parses as typed_vcat on 0.7.
        @variable(model, named_two_dim[j = 1:10, k=3:2:9;
                                       isodd(j + k) && k <= 8])
        # Parses as vect on 0.7.
        anon_two_dim = @variable(model, [j = 1:10, k=3:2:9;
                                         isodd(j + k) && k <= 8])

        test_one_dim(named_one_dim)
        test_one_dim(anon_one_dim)
        test_two_dim(named_two_dim)
        test_two_dim(anon_two_dim)
    end

    @testset "@variable returning Array{VariableRef}" begin
        m = ModelType()
        @variable(m, x[1:3,1:4,1:2], start = 0.0)
        @variable(m, y[1:0], start = 0.0)
        @variable(m, z[1:4], start = 0.0)

        @test typeof(x) == Array{VariableRefType,3}
        @test typeof(y) == Vector{VariableRefType}
        @test typeof(z) == Vector{VariableRefType}

        @test typeof(JuMP.start_value.(x)) == Array{Float64,3}
        # No type to infer for an empty collection.
        @test typeof(JuMP.start_value.(y)) == Array{Any,1}
        @test typeof(JuMP.start_value.(z)) == Array{Float64,1}
    end

    @testset "start_value on empty things" begin
        m = ModelType()
        @variable(m, x[1:4,  1:0,1:3], start = 0) # Array{VariableRef}
        @variable(m, y[1:4,  2:1,1:3], start = 0) # JuMPArray
        @variable(m, z[1:4,Set(),1:3], start = 0) # Dict

        @test JuMP.start_value.(x) == Array{Float64}(undef, 4, 0, 3)
        # TODO: Decide what to do here. I don't know if we still need to test this given broadcast syntax.
        # @test typeof(JuMP.start_value(y)) <: JuMP.JuMPArray{Float64}
        # @test JuMP.size(JuMP.start_value(y)) == (4,0,3)
        # @test typeof(JuMP.start_value(z)) == JuMP.JuMPArray{Float64,3,Tuple{UnitRange{Int},Set{Any},UnitRange{Int}}}
        # @test length(JuMP.start_value(z)) == 0
    end

# Slices three-dimensional JuMPArray x[I,J,K]
# I,J,K can be singletons, ranges, colons, etc.
function sliceof(x, I, J, K)
    y = Array{VariableRefType}(undef, length(I), length(J), length(K))

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
    @static if VERSION < v"0.7.0-"
        squeeze(y, tuple(findall(idx)...))
    else
        dropdims(y, dims=tuple(findall(idx)...))
    end
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

    @testset "end for indexing a JuMPArray" begin
        model = ModelType()
        @variable(model, x[0:2,1:4])
        @variable(model, z[0:2])

        if VERSION >= v"0.7-"
            @test x[end,1] == x[2,1]
            @test x[0,end-1] == x[0,3]
            @test z[end] == z[2]
            # TODO: This probably isn't hard to make work.
            @test_throws ErrorException x[end-1]
        else
            @test_throws ErrorException x[end,1]
            @test_throws ErrorException x[end-1]
            @test_throws ErrorException x[0,end-1]
            @test_throws ErrorException z[end]
        end
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
    variables_test(Model, VariableRef)
end

@testset "Variables for JuMPExtension.MyModel" begin
    variables_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end
