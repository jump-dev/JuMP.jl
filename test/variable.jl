#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/variable.jl
# Testing for VariableRef
#############################################################################

module VariableTests

using JuMP

import LinearAlgebra: Symmetric
using Pukeko

include("utilities.jl")
include("JuMPExtension.jl")
const BOTH_MODEL_AND_VAR = [
    (Model, VariableRef),
    (JuMPExtension.MyModel, JuMPExtension.MyVariableRef)]
const BOTH_MODEL = [Model, JuMPExtension.MyModel]

function test_name(variable, name)
    @test JuMP.name(variable) == name
    @test variable == JuMP.variable_by_name(JuMP.owner_model(variable), name)
end

# Slices three-dimensional DenseAxisArray x[I,J,K]
# I,J,K can be singletons, ranges, colons, etc.
function sliceof(VariableRefType, x, I, J, K)
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
    dropdims(y, dims=tuple(findall(idx)...))
end

function variable_no_bound(ModelType, VariableRefType)
    model = ModelType()
    @variable(model, nobounds)
    @test !JuMP.has_lower_bound(nobounds)
    @test !JuMP.has_upper_bound(nobounds)
    @test !JuMP.is_fixed(nobounds)
    test_name(nobounds, "nobounds")
    @test zero(nobounds) isa JuMP.GenericAffExpr{Float64, VariableRefType}
    @test one(nobounds) isa JuMP.GenericAffExpr{Float64, VariableRefType}
end
@parametric variable_no_bound BOTH_MODEL_AND_VAR

function variable_lower_bound_rhs(ModelType)
    model = ModelType()
    @variable(model, lbonly >= 0, Bin)
    @test JuMP.has_lower_bound(lbonly)
    @test JuMP.lower_bound(lbonly) == 0.0
    @test !JuMP.has_upper_bound(lbonly)
    @test !JuMP.is_fixed(lbonly)
    @test JuMP.is_binary(lbonly)
    @test !JuMP.is_integer(lbonly)
    @test isequal(model[:lbonly], lbonly)
    JuMP.delete_lower_bound(lbonly)
    @test !JuMP.has_lower_bound(lbonly)
    # Name already used
    @test_throws ErrorException @variable(model, lbonly)
end
@parametric variable_lower_bound_rhs BOTH_MODEL

function variable_lower_bound_lhs(ModelType)
    model = ModelType()
    @variable(model, 0 <= lblhs, Bin)
    @test JuMP.has_lower_bound(lblhs)
    @test JuMP.lower_bound(lblhs) == 0.0
    @test !JuMP.has_upper_bound(lblhs)
    @test !JuMP.is_fixed(lblhs)
    @test JuMP.is_binary(lblhs)
    @test !JuMP.is_integer(lblhs)
    @test isequal(model[:lblhs], lblhs)
end
@parametric variable_lower_bound_lhs BOTH_MODEL

function variable_upper_bound_rhs(ModelType)
    model = ModelType()
    @variable(model, ubonly <= 1, Int)
    @test !JuMP.has_lower_bound(ubonly)
    @test JuMP.has_upper_bound(ubonly)
    @test JuMP.upper_bound(ubonly) == 1.0
    @test !JuMP.is_fixed(ubonly)
    @test !JuMP.is_binary(ubonly)
    @test JuMP.is_integer(ubonly)
    @test isequal(model[:ubonly], ubonly)
    JuMP.delete_upper_bound(ubonly)
    @test !JuMP.has_upper_bound(ubonly)
end
@parametric variable_upper_bound_rhs BOTH_MODEL

function variable_upper_bound_lhs(ModelType)
    model = ModelType()
    @variable(model, 1 >= ublhs, Int)
    @test !JuMP.has_lower_bound(ublhs)
    @test JuMP.has_upper_bound(ublhs)
    @test JuMP.upper_bound(ublhs) == 1.0
    @test !JuMP.is_fixed(ublhs)
    @test !JuMP.is_binary(ublhs)
    @test JuMP.is_integer(ublhs)
    @test isequal(model[:ublhs],ublhs)
end
@parametric variable_upper_bound_lhs BOTH_MODEL

function variable_interval(ModelType)
    function has_bounds(var, lb, ub)
        @test JuMP.has_lower_bound(var)
        @test JuMP.lower_bound(var) == lb
        @test JuMP.has_upper_bound(var)
        @test JuMP.upper_bound(var) == ub
        @test !JuMP.is_fixed(var)
    end
    model = ModelType()
    @variable(model, 0 <= bothb1 <= 1)
    has_bounds(bothb1, 0.0, 1.0)
    @variable(model, 0 ≤  bothb2 ≤  1)
    has_bounds(bothb2, 0.0, 1.0)
    @variable(model, 1 >= bothb3 >= 0)
    has_bounds(bothb3, 0.0, 1.0)
    @variable(model, 1 ≥  bothb4 ≥  0)
    has_bounds(bothb4, 0.0, 1.0)
    @test_macro_throws ErrorException @variable(model, 1 ≥ bothb5 ≤ 0)
    @test_macro_throws ErrorException @variable(model, 1 ≤ bothb6 ≥ 0)
end
@parametric variable_interval BOTH_MODEL

function variable_fix(ModelType)
    model = ModelType()
    @variable(model, fixed == 1.0)
    @test !JuMP.has_lower_bound(fixed)
    @test !JuMP.has_upper_bound(fixed)
    @test JuMP.is_fixed(fixed)
    @test JuMP.fix_value(fixed) == 1.0
    JuMP.unfix(fixed)
    @test !JuMP.is_fixed(fixed)
    JuMP.set_lower_bound(fixed, 0.0)
    @test_throws Exception JuMP.fix(fixed, 1.0)
    JuMP.fix(fixed, 1.0; force = true)
    @test !JuMP.has_lower_bound(fixed)
    @test !JuMP.has_upper_bound(fixed)
    @test JuMP.is_fixed(fixed)
    @test JuMP.fix_value(fixed) == 1.0
end
@parametric variable_fix BOTH_MODEL

function variable_custom_index_sets(ModelType)
    model = ModelType()
    @variable(model, onerangeub[-7:1] <= 10, Int)
    @variable(model, manyrangelb[0:1, 10:20, 1:1] >= 2)
    @test JuMP.has_lower_bound(manyrangelb[0, 15, 1])
    @test JuMP.lower_bound(manyrangelb[0, 15, 1]) == 2
    @test !JuMP.has_upper_bound(manyrangelb[0, 15, 1])

    s = ["Green","Blue"]
    @variable(model, x[i=-10:10, s] <= 5.5, Int, start=i+1)
    @test JuMP.upper_bound(x[-4, "Green"]) == 5.5
    test_name(x[-10, "Green"], "x[-10,Green]")
    # TODO: broken because of
    #       https://github.com/JuliaOpt/MathOptInterface.jl/issues/302
    # @test JuMP.start_value(x[-3, "Blue"]) == -2
    @test isequal(model[:onerangeub][-7], onerangeub[-7])
    @test_throws KeyError model[:foo]
end
@parametric variable_custom_index_sets BOTH_MODEL

function variable_anonymous(ModelType)
    model = ModelType()
    @test_throws ErrorException @variable(model, [(0, 0)])  # #922
    x = @variable(model, [(0, 2)])
    @test JuMP.name(x[0]) == ""
    @test JuMP.name(x[2]) == ""
end
@parametric variable_anonymous BOTH_MODEL

function variable_is_valid_delete(ModelType)
    model = ModelType()
    @variable(model, x)
    @test JuMP.is_valid(model, x)
    JuMP.delete(model, x)
    @test !JuMP.is_valid(model, x)
    second_model = ModelType()
    @test_throws Exception JuMP.delete(second_model, x)
end
@parametric variable_is_valid_delete BOTH_MODEL

function variable_bounds_set_get(ModelType)
    model = ModelType()
    @variable(model, 0 <= x <= 2)
    @test JuMP.lower_bound(x) == 0
    @test JuMP.upper_bound(x) == 2
    set_lower_bound(x, 1)
    @test JuMP.lower_bound(x) == 1
    set_upper_bound(x, 3)
    @test JuMP.upper_bound(x) == 3
    @variable(model, q, Bin)
    @test !JuMP.has_lower_bound(q)
    @test !JuMP.has_upper_bound(q)

    @variable(model, 0 <= y <= 1, Bin)
    @test JuMP.lower_bound(y) == 0
    @test JuMP.upper_bound(y) == 1

    @variable(model, fixedvar == 2)
    @test JuMP.fix_value(fixedvar) == 2.0
    JuMP.fix(fixedvar, 5)
    @test JuMP.fix_value(fixedvar) == 5
    @test_throws Exception JuMP.lower_bound(fixedvar)
    @test_throws Exception JuMP.upper_bound(fixedvar)
end
@parametric variable_bounds_set_get BOTH_MODEL

function variable_starts_set_get(ModelType)
    model = ModelType()
    @variable(model, x[1:3])
    x0 = collect(1:3)
    JuMP.set_start_value.(x, x0)
    @test JuMP.start_value.(x) == x0
    @test JuMP.start_value.([x[1],x[2],x[3]]) == x0

    @variable(model, y[1:3,1:2])
    @test_throws DimensionMismatch JuMP.set_start_value.(y, collect(1:6))
end
@parametric variable_starts_set_get BOTH_MODEL

function variable_integrality_set_get(ModelType)
    model = ModelType()
    @variable(model, x[1:3])

    JuMP.set_integer(x[2])
    JuMP.set_integer(x[2])  # test duplicated call
    @test JuMP.is_integer(x[2])
    JuMP.unset_integer(x[2])
    @test !JuMP.is_integer(x[2])

    JuMP.set_binary(x[1])
    JuMP.set_binary(x[1])  # test duplicated call
    @test JuMP.is_binary(x[1])
    @test_throws Exception JuMP.set_integer(x[1])
    JuMP.unset_binary(x[1])
    @test !JuMP.is_binary(x[1])

    @variable(model, y, binary = true)
    @test JuMP.is_binary(y)
    @test_throws Exception JuMP.set_integer(y)
    JuMP.unset_binary(y)
    @test !JuMP.is_binary(y)

    @variable(model, z, integer = true)
    @test JuMP.is_integer(z)
    @test_throws Exception JuMP.set_binary(z)
    JuMP.unset_integer(z)
    @test !JuMP.is_integer(z)
end
@parametric variable_integrality_set_get BOTH_MODEL

function variable_repeated_elements(ModelType)
    # Tests repeated elements in index set throw error (JuMP issue #199).
    model = ModelType()
    index_set = [:x,:x,:y]
    @test_throws ErrorException (
        @variable(model, unused_variable[index_set], container=DenseAxisArray))
    @test_throws ErrorException (
        @variable(model, unused_variable[index_set], container=SparseAxisArray))
    @test_throws ErrorException (
        @variable(model, unused_variable[index_set, [1]], container=DenseAxisArray))
    @test_throws ErrorException (
        @variable(model, unused_variable[index_set, [1]], container=SparseAxisArray))
end
@parametric variable_repeated_elements BOTH_MODEL

function variable_oneto_index_set(ModelType, VariableRefType)
    # Tests that Base.OneTo can be used in index set (JuMP issue #933).
    model = ModelType()
    auto_var = @variable(model, [Base.OneTo(3), 1:2], container=Auto)
    @test auto_var isa Matrix{VariableRefType}
    @test size(auto_var) == (3, 2)
    array_var = @variable(model, [Base.OneTo(3), 1:2], container=Array)
    @test array_var isa Matrix{VariableRefType}
    @test size(array_var) == (3, 2)
    denseaxisarray_var = @variable(model, [Base.OneTo(3), 1:2], container=DenseAxisArray)
    @test denseaxisarray_var isa JuMP.Containers.DenseAxisArray{VariableRefType}
    @test length.(axes(denseaxisarray_var)) == (3, 2)
end
@parametric variable_oneto_index_set BOTH_MODEL_AND_VAR

function variable_base_name_in_macro(ModelType)
    model = ModelType()
    @variable(model, normal_var)
    test_name(normal_var, "normal_var")
    no_indices = @variable(model, base_name="foo")
    test_name(no_indices, "foo")
    # Note that `z` will be ignored in name.
    indices = @variable(model, z[i=2:3], base_name="t")
    test_name(indices[2], "t[2]")
    test_name(indices[3], "t[3]")
end
@parametric variable_base_name_in_macro BOTH_MODEL

function variable_name(ModelType)
    model = ModelType()
    @variable(model, x)
    test_name(x, "x")
    JuMP.set_name(x, "y")
    @test JuMP.variable_by_name(model, "x") isa Nothing
    test_name(x, "y")
    y = @variable(model, base_name="y")
    err(name) = ErrorException("Multiple variables have the name $name.")
    @test_throws err("y") JuMP.variable_by_name(model, "y")
    JuMP.set_name(y, "x")
    test_name(x, "y")
    test_name(y, "x")
    JuMP.set_name(x, "x")
    @test_throws err("x") JuMP.variable_by_name(model, "x")
    @test JuMP.variable_by_name(model, "y") isa Nothing
    JuMP.set_name(y, "y")
    test_name(x, "x")
    test_name(y, "y")
end
@parametric variable_name BOTH_MODEL

function variable_condition_in_indexing(ModelType)
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
    @variable(model, named_one_dim[i=1:10; iseven(i)])
    test_one_dim(named_one_dim)
    # Parses as vcat on 0.7.
    anon_one_dim = @variable(model, [i=1:10; iseven(i)])
    test_one_dim(anon_one_dim)
    # Parses as typed_vcat on 0.7.
    @variable(model, named_two_dim[j=1:10, k=3:2:9; isodd(j + k) && k <= 8])
    test_two_dim(named_two_dim)
    # Parses as vect on 0.7.
    anon_two_dim = @variable(model, [j=1:10, k=3:2:9; isodd(j + k) && k <= 8])
    test_two_dim(anon_two_dim)
end
@parametric variable_condition_in_indexing BOTH_MODEL

function variable_macro_return_type(ModelType, VariableRefType)
    model = ModelType()
    @variable(model, x[1:3, 1:4, 1:2], start=0.0)
    @test typeof(x) == Array{VariableRefType,3}
    @test typeof(JuMP.start_value.(x)) == Array{Float64,3}

    @variable(model, y[1:0], start=0.0)
    @test typeof(y) == Vector{VariableRefType}
    # No type to infer for an empty collection.
    @test typeof(JuMP.start_value.(y)) == Array{Any,1}

    @variable(model, z[1:4], start = 0.0)
    @test typeof(z) == Vector{VariableRefType}
    @test typeof(JuMP.start_value.(z)) == Array{Float64,1}
end
@parametric variable_macro_return_type BOTH_MODEL_AND_VAR

function variable_start_value_on_empty(ModelType)
    model = ModelType()
    @variable(model, x[1:4,  1:0,1:3], start = 0)  # Array{VariableRef}
    @variable(model, y[1:4,  2:1,1:3], start = 0)  # DenseAxisArray
    @variable(model, z[1:4,Set(),1:3], start = 0)  # SparseAxisArray

    @test JuMP.start_value.(x) == Array{Float64}(undef, 4, 0, 3)
    # TODO: Decide what to do here. I don't know if we still need to test this
    #       given broadcast syntax.
    # @test typeof(JuMP.start_value(y)) <: JuMP.DenseAxisArray{Float64}
    # @test JuMP.size(JuMP.start_value(y)) == (4,0,3)
    # @test typeof(JuMP.start_value(z)) ==
    #   JuMP.DenseAxisArray{Float64,3,Tuple{UnitRange{Int},Set{Any},UnitRange{Int}}}
    # @test length(JuMP.start_value(z)) == 0
end
@parametric variable_start_value_on_empty BOTH_MODEL

function variable_denseaxisarray_slices(ModelType, VariableRefType)
    # Test slicing DenseAxisArrays (JuMP issue #684).
    model = ModelType()
    @variable(model, x[1:3, 1:4, 1:2], container=DenseAxisArray)
    @variable(model, y[1:3, -1:2, 3:4])
    @variable(model, z[1:3, -1:2:4, 3:4])
    @variable(model, w[1:3, -1:2,[:red, "blue"]])

    #@test x[:] == vec(sliceof(VariableRefType, x, 1:3, 1:4, 1:2))
    @test x isa JuMP.Containers.DenseAxisArray
    @test x[:, :, :].data == sliceof(VariableRefType, x, 1:3, 1:4, 1:2)
    @test x[1, :, :].data == sliceof(VariableRefType, x, 1, 1:4, 1:2)
    @test x[1, :, 2].data == sliceof(VariableRefType, x, 1, 1:4, 2)
    @test_throws KeyError x[1, :, 3]
    #@test x[1:2,:,:].data == sliceof(VariableRefType, x, 1:2, 1:4, 1:2)
    #@test x[1:2,:,2].data == sliceof(VariableRefType, x, 1:2, 1:4, 2)
    #@test x[1:2,:,1:2].data == sliceof(VariableRefType, x, 1:2, 1:4, 1:2)
    @test_throws KeyError x[1:2, :, 1:3]

    #@test y[:] == vec(sliceof(VariableRefType, y, 1:3, -1:2, 3:4))
    @test y[:, :, :].data == sliceof(VariableRefType, y, 1:3, -1:2, 3:4)
    @test y[1, :, :].data == sliceof(VariableRefType, y, 1, -1:2, 3:4)
    @test y[1, :, 4].data == sliceof(VariableRefType, y, 1, -1:2, 4)
    @test_throws KeyError y[1, :, 5]
    # @test y[1:2,:,:] == sliceof(VariableRefType, y, 1:2, -1:2, 3:4)
    # @test y[1:2,:,4] == sliceof(VariableRefType, y, 1:2, -1:2, 4)
    # @test y[1:2,:,3:4] == sliceof(VariableRefType, y, 1:2, -1:2, 3:4)
    # @test_throws BoundsError y[1:2,:,1:3]

    #@test z[:] == vec(sliceof(VariableRefType, z, 1:3, -1:2:4, 3:4))
    @test z[:, 1, :].data == sliceof(VariableRefType, z, 1:3, 1, 3:4)
    @test z[1, 1, :].data == sliceof(VariableRefType, z, 1, 1, 3:4)
    @test_throws KeyError z[:, 5, 3]
    # @test z[1:2,1,:] == sliceof(VariableRefType, z, 1:2, 1, 3:4)
    # @test z[1:2,1,4] == sliceof(VariableRefType, z, 1:2, 1, 4)
    # @test z[1:2,1,3:4] == sliceof(VariableRefType, z, 1:2, 1, 3:4)
    # @test_throws BoundsError z[1:2,1,1:3]

    #@test w[:] == vec(sliceof(VariableRefType, w, 1:3, -1:2, [:red,"blue"]))
    @test w[:, :, :] == w
    @test w[1, :, "blue"].data == sliceof(VariableRefType, w, 1, -1:2, ["blue"])
    @test w[1, :, :red].data == sliceof(VariableRefType, w, 1, -1:2, [:red])
    @test_throws KeyError w[1, :, "green"]
    # @test w[1:2,:,"blue"] == sliceof(VariableRefType, w, 1:2, -1:2, ["blue"])
    # @test_throws ErrorException w[1:2,:,[:red,"blue"]]
end
@parametric variable_denseaxisarray_slices BOTH_MODEL_AND_VAR

function variable_end_indexing(ModelType)
    model = ModelType()
    @variable(model, x[0:2, 1:4])
    @variable(model, z[0:2])
    @test x[end,1] == x[2, 1]
    @test x[0, end-1] == x[0, 3]
    @test z[end] == z[2]
    # TODO: It is redirected to x[11] as it is the 11th element but linear
    #       indexing is not supported
    @test_throws KeyError x[end-1]
end
@parametric variable_end_indexing BOTH_MODEL

function variable_unsigned_index(ModelType)
    # Tests unsigned int can be used to construct index set (JuMP issue #857).
    model = ModelType()
    t = UInt(4)
    @variable(model, x[1:t])
    @test JuMP.num_variables(model) == 4
end
@parametric variable_unsigned_index BOTH_MODEL

function variable_symmetric(ModelType)
    model = ModelType()

    @variable model x[1:2, 1:2] Symmetric
    @test x isa Symmetric
    @test x[1, 2] === x[2, 1]

    y = @variable model [1:2, 1:2] Symmetric
    @test y isa Symmetric
    @test y[1, 2] === y[2, 1]
end
@parametric variable_symmetric BOTH_MODEL

end  # module VariableTests

import Pukeko
Pukeko.run_tests(VariableTests, exclude_name="test_name")