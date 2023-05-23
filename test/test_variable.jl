#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# test/variable.jl
# Testing for VariableRef
#############################################################################

module TestVariable

using JuMP
using Test

import LinearAlgebra

include(joinpath(@__DIR__, "utilities.jl"))

function _test_variable_name_util(variable, s_name)
    @test s_name == @inferred name(variable)
    @test variable == variable_by_name(owner_model(variable), s_name)
    return
end

# Slices three-dimensional DenseAxisArray x[I,J,K]
# I,J,K can be singletons, ranges, colons, etc.
function _sliceof_util(VariableRefType, x, I, J, K)
    y = Array{VariableRefType}(undef, length(I), length(J), length(K))
    ii = 1
    for i in I
        jj = 1
        for j in J
            kk = 1
            for k in K
                y[ii, jj, kk] = x[i, j, k]
                kk += 1
            end
            jj += 1
        end
        ii += 1
    end
    idx = [length(I) == 1, length(J) == 1, length(K) == 1]
    return dropdims(y; dims = tuple(findall(idx)...))
end

function test_extension_variable_no_bound(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, nobounds)
    @test !has_lower_bound(nobounds)
    @test !has_upper_bound(nobounds)
    @test !is_fixed(nobounds)
    _test_variable_name_util(nobounds, "nobounds")
    @test zero(nobounds) isa GenericAffExpr{T,VariableRefType}
    @test one(nobounds) isa GenericAffExpr{T,VariableRefType}
    return
end

function test_extension_variable_lower_bound_rhs(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, lbonly >= 0, Bin)
    @test has_lower_bound(lbonly)
    @test 0.0 == @inferred lower_bound(lbonly)
    @test !has_upper_bound(lbonly)
    @test !is_fixed(lbonly)
    @test is_binary(lbonly)
    @test !is_integer(lbonly)
    @test isequal(model[:lbonly], lbonly)
    delete_lower_bound(lbonly)
    @test !has_lower_bound(lbonly)
    # Name already used
    @test_throws ErrorException @variable(model, lbonly)
    return
end

function test_extension_variable_lower_bound_lhs(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, 0 <= lblhs, Bin)
    @test has_lower_bound(lblhs)
    @test 0.0 == @inferred lower_bound(lblhs)
    @test !has_upper_bound(lblhs)
    @test !is_fixed(lblhs)
    @test is_binary(lblhs)
    @test !is_integer(lblhs)
    @test isequal(model[:lblhs], lblhs)
    return
end

function test_extension_variable_upper_bound_rhs(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, ubonly <= 1, Int)
    @test !has_lower_bound(ubonly)
    @test has_upper_bound(ubonly)
    @test 1.0 == @inferred upper_bound(ubonly)
    @test !is_fixed(ubonly)
    @test !is_binary(ubonly)
    @test is_integer(ubonly)
    @test isequal(model[:ubonly], ubonly)
    delete_upper_bound(ubonly)
    @test !has_upper_bound(ubonly)
    return
end

function test_extension_variable_upper_bound_lhs(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, 1 >= ublhs, Int)
    @test !has_lower_bound(ublhs)
    @test has_upper_bound(ublhs)
    @test 1.0 == @inferred upper_bound(ublhs)
    @test !is_fixed(ublhs)
    @test !is_binary(ublhs)
    @test is_integer(ublhs)
    @test isequal(model[:ublhs], ublhs)
    return
end

function _has_bounds(var, lb, ub)
    @test has_lower_bound(var)
    @test lb == @inferred lower_bound(var)
    @test has_upper_bound(var)
    @test ub == @inferred upper_bound(var)
    @test !is_fixed(var)
    return
end

function test_extension_variable_interval(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, 0 <= bothb1 <= 1)
    _has_bounds(bothb1, 0.0, 1.0)
    @variable(model, 0 ≤ bothb2 ≤ 1)
    _has_bounds(bothb2, 0.0, 1.0)
    @variable(model, 1 >= bothb3 >= 0)
    _has_bounds(bothb3, 0.0, 1.0)
    @variable(model, 1 ≥ bothb4 ≥ 0)
    _has_bounds(bothb4, 0.0, 1.0)
    @test_macro_throws ErrorException @variable(model, 1 ≥ bothb5 ≤ 0)
    @test_macro_throws ErrorException @variable(model, 1 ≤ bothb6 ≥ 0)
    return
end

function test_extension_variable_fix(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, fixed == 1.0)
    @test !has_lower_bound(fixed)
    @test !has_upper_bound(fixed)
    @test is_fixed(fixed)
    @test 1.0 == @inferred fix_value(fixed)
    unfix(fixed)
    @test !is_fixed(fixed)
    set_lower_bound(fixed, 0.0)
    @test_throws Exception fix(fixed, 1.0)
    fix(fixed, 1.0; force = true)
    @test !has_lower_bound(fixed)
    @test !has_upper_bound(fixed)
    @test is_fixed(fixed)
    @test 1.0 == @inferred fix_value(fixed)
    return
end

function test_extension_variable_custom_index_sets(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, onerangeub[-7:1] <= 10, Int)
    @variable(model, manyrangelb[0:1, 10:20, 1:1] >= 2)
    @test has_lower_bound(manyrangelb[0, 15, 1])
    @test 2 == @inferred lower_bound(manyrangelb[0, 15, 1])
    @test !has_upper_bound(manyrangelb[0, 15, 1])

    s = ["Green", "Blue"]
    @variable(model, x[i = -10:10, s] <= 5.5, Int, start = i + 1)
    @test 5.5 == @inferred upper_bound(x[-4, "Green"])
    _test_variable_name_util(x[-10, "Green"], "x[-10,Green]")
    @test start_value(x[-3, "Blue"]) == -2
    @test isequal(model[:onerangeub][-7], onerangeub[-7])
    @test_throws KeyError model[:foo]
    return
end

function test_extension_variable_anonymous(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @test_throws ErrorException @variable(model, [(0, 0)])  # #922
    x = @variable(model, [(0, 2)])
    @test "" == @inferred name(x[0])
    @test "" == @inferred name(x[2])
    return
end

function test_extension_variable_is_valid_delete(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test is_valid(model, x)
    delete(model, x)
    @test !is_valid(model, x)
    second_model = ModelType()
    @test_throws Exception delete(second_model, x)
    return
end

function test_extension_variable_bounds_set_get(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, 0 <= x <= 2)
    @test 0 == @inferred lower_bound(x)
    @test 2 == @inferred upper_bound(x)
    set_lower_bound(x, 1)
    @test 1 == @inferred lower_bound(x)
    set_upper_bound(x, 3)
    @test 3 == @inferred upper_bound(x)
    @variable(model, q, Bin)
    @test !has_lower_bound(q)
    @test !has_upper_bound(q)
    @variable(model, 0 <= y <= 1, Bin)
    @test 0 == @inferred lower_bound(y)
    @test 1 == @inferred upper_bound(y)
    @variable(model, fixedvar == 2)
    @test 2.0 == @inferred fix_value(fixedvar)
    fix(fixedvar, 5)
    @test 5 == @inferred fix_value(fixedvar)
    @test_throws Exception lower_bound(fixedvar)
    @test_throws Exception upper_bound(fixedvar)
    return
end

function test_extension_variable_starts_set_get(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:3])
    x0 = collect(1:3)
    set_start_value.(x, x0)
    @test start_value.(x) == x0
    @test start_value.([x[1], x[2], x[3]]) == x0
    @test has_start_value(x[1]) == true
    @variable(model, y[1:3, 1:2])
    @test has_start_value(y[1, 1]) == false
    @test_throws DimensionMismatch set_start_value.(y, collect(1:6))
    return
end

function test_extension_variable_integrality_set_get(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:3])
    set_integer(x[2])
    set_integer(x[2])  # test duplicated call
    @test is_integer(x[2])
    unset_integer(x[2])
    @test !is_integer(x[2])
    set_binary(x[1])
    set_binary(x[1])  # test duplicated call
    @test is_binary(x[1])
    @test_throws Exception set_integer(x[1])
    unset_binary(x[1])
    @test !is_binary(x[1])
    @variable(model, y, binary = true)
    @test is_binary(y)
    @test_throws Exception set_integer(y)
    unset_binary(y)
    @test !is_binary(y)
    @variable(model, z, integer = true)
    @test is_integer(z)
    @test_throws Exception set_binary(z)
    unset_integer(z)
    @test !is_integer(z)
    return
end

function test_extension_variable_repeated_elements(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    # Tests repeated elements in index set throw error (JuMP issue #199).
    model = ModelType()
    index_set = [:x, :x, :y]
    @test_throws ErrorException @variable(
        model,
        unused_variable[index_set],
        container = DenseAxisArray
    )
    @test_throws ErrorException @variable(
        model,
        unused_variable[index_set],
        container = SparseAxisArray
    )
    @test_throws ErrorException @variable(
        model,
        unused_variable[index_set, [1]],
        container = DenseAxisArray
    )
    @test_throws ErrorException @variable(
        model,
        unused_variable[index_set, [1]],
        container = SparseAxisArray
    )
    return
end

function test_extension_variable_oneto_index_set(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    # Tests that Base.OneTo can be used in index set (JuMP issue #933).
    model = ModelType()
    auto_var = @variable(model, [Base.OneTo(3), 1:2], container = Auto)
    @test auto_var isa Matrix{VariableRefType}
    @test (3, 2) == @inferred size(auto_var)
    array_var = @variable(model, [Base.OneTo(3), 1:2], container = Array)
    @test array_var isa Matrix{VariableRefType}
    @test (3, 2) == @inferred size(array_var)
    denseaxisarray_var =
        @variable(model, [Base.OneTo(3), 1:2], container = DenseAxisArray)
    @test denseaxisarray_var isa Containers.DenseAxisArray{VariableRefType}
    @test length.(axes(denseaxisarray_var)) == (3, 2)
    return
end

function test_extension_variable_base_name_in_macro(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, normal_var)
    _test_variable_name_util(normal_var, "normal_var")
    no_indices = @variable(model, base_name = "foo")
    _test_variable_name_util(no_indices, "foo")
    # Note that `z` will be ignored in name.
    indices = @variable(model, z[i = 2:3], base_name = "t")
    _test_variable_name_util(indices[2], "t[2]")
    _test_variable_name_util(indices[3], "t[3]")
    return
end

function test_extension_variable_name(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    _test_variable_name_util(x, "x")
    set_name(x, "y")
    @test variable_by_name(model, "x") isa Nothing
    _test_variable_name_util(x, "y")
    y = @variable(model, base_name = "y")
    err(name) = ErrorException("Multiple variables have the name $name.")
    @test_throws err("y") variable_by_name(model, "y")
    set_name(y, "x")
    _test_variable_name_util(x, "y")
    _test_variable_name_util(y, "x")
    set_name(x, "x")
    @test_throws err("x") variable_by_name(model, "x")
    @test variable_by_name(model, "y") isa Nothing
    set_name(y, "y")
    _test_variable_name_util(x, "x")
    _test_variable_name_util(y, "y")
    return
end

function test_extension_variable_condition_in_indexing(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    function test_one_dim(x)
        @test 5 == @inferred length(x)
        for i in 1:10
            @test haskey(x, i) == iseven(i)
        end
        return
    end
    function test_two_dim(y)
        @test 15 == @inferred length(y)
        for j in 1:10, k in 3:2:9
            @test haskey(y, (j, k)) == (isodd(j + k) && k <= 8)
        end
        return
    end
    model = ModelType()
    # Parses as ref on 0.7.
    @variable(model, named_one_dim[i = 1:10; iseven(i)])
    test_one_dim(named_one_dim)
    # Parses as vcat on 0.7.
    anon_one_dim = @variable(model, [i = 1:10; iseven(i)])
    test_one_dim(anon_one_dim)
    # Parses as typed_vcat on 0.7.
    @variable(model, named_two_dim[j = 1:10, k = 3:2:9; isodd(j + k) && k <= 8])
    test_two_dim(named_two_dim)
    # Parses as vect on 0.7.
    anon_two_dim =
        @variable(model, [j = 1:10, k = 3:2:9; isodd(j + k) && k <= 8])
    test_two_dim(anon_two_dim)
    return
end

function test_extension_variable_macro_return_type(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x[1:3, 1:4, 1:2], start = T(0))
    @test typeof(x) == Array{VariableRefType,3}
    @test typeof(start_value.(x)) == Array{T,3}
    @variable(model, y[1:0], start = T(0))
    @test typeof(y) == Vector{VariableRefType}
    # No type to infer for an empty collection.
    @test typeof(start_value.(y)) == Vector{Union{Nothing,T}}
    @variable(model, z[1:4], start = T(0))
    @test typeof(z) == Vector{VariableRefType}
    @test typeof(start_value.(z)) == Vector{T}
    return
end

function test_extension_variable_start_value_on_empty(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x[1:4, 1:0, 1:3], start = 0)  # Array{VariableRef}
    @variable(model, y[1:4, 2:1, 1:3], start = 0)  # DenseAxisArray
    @variable(model, z[1:4, Set(), 1:3], start = 0)  # SparseAxisArray

    @test start_value.(x) == Array{T}(undef, 4, 0, 3)
    # TODO: Decide what to do here. I don't know if we still need to test this
    #       given broadcast syntax.
    # @test typeof(start_value(y)) <: DenseAxisArray{Float64}
    # @test size(start_value(y)) == (4,0,3)
    # @test typeof(start_value(z)) ==
    #   DenseAxisArray{Float64,3,Tuple{UnitRange{Int},Set{Any},UnitRange{Int}}}
    # @test length(start_value(z)) == 0
    return
end

function test_extension_variable_denseaxisarray_slices(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    # Test slicing DenseAxisArrays (JuMP issue #684).
    model = ModelType()
    @variable(model, x[1:3, 1:4, 1:2], container = DenseAxisArray)
    @variable(model, y[1:3, -1:2, 3:4])
    @variable(model, z[1:3, -1:2:4, 3:4])
    @variable(model, w[1:3, -1:2, [:red, "blue"]])

    #@test x[:] == vec(_sliceof_util(VariableRefType, x, 1:3, 1:4, 1:2))
    @test x isa Containers.DenseAxisArray
    @test x[:, :, :].data == _sliceof_util(VariableRefType, x, 1:3, 1:4, 1:2)
    @test x[1, :, :].data == _sliceof_util(VariableRefType, x, 1, 1:4, 1:2)
    @test x[1, :, 2].data == _sliceof_util(VariableRefType, x, 1, 1:4, 2)
    @test_throws KeyError x[1, :, 3]
    #@test x[1:2,:,:].data == _sliceof_util(VariableRefType, x, 1:2, 1:4, 1:2)
    #@test x[1:2,:,2].data == _sliceof_util(VariableRefType, x, 1:2, 1:4, 2)
    #@test x[1:2,:,1:2].data == _sliceof_util(VariableRefType, x, 1:2, 1:4, 1:2)
    @test_throws KeyError x[1:2, :, 1:3]

    #@test y[:] == vec(_sliceof_util(VariableRefType, y, 1:3, -1:2, 3:4))
    @test y[:, :, :].data == _sliceof_util(VariableRefType, y, 1:3, -1:2, 3:4)
    @test y[1, :, :].data == _sliceof_util(VariableRefType, y, 1, -1:2, 3:4)
    @test y[1, :, 4].data == _sliceof_util(VariableRefType, y, 1, -1:2, 4)
    @test_throws KeyError y[1, :, 5]
    # @test y[1:2,:,:] == _sliceof_util(VariableRefType, y, 1:2, -1:2, 3:4)
    # @test y[1:2,:,4] == _sliceof_util(VariableRefType, y, 1:2, -1:2, 4)
    # @test y[1:2,:,3:4] == _sliceof_util(VariableRefType, y, 1:2, -1:2, 3:4)
    # @test_throws BoundsError y[1:2,:,1:3]

    #@test z[:] == vec(_sliceof_util(VariableRefType, z, 1:3, -1:2:4, 3:4))
    @test z[:, 1, :].data == _sliceof_util(VariableRefType, z, 1:3, 1, 3:4)
    @test z[1, 1, :].data == _sliceof_util(VariableRefType, z, 1, 1, 3:4)
    @test_throws KeyError z[:, 5, 3]
    # @test z[1:2,1,:] == _sliceof_util(VariableRefType, z, 1:2, 1, 3:4)
    # @test z[1:2,1,4] == _sliceof_util(VariableRefType, z, 1:2, 1, 4)
    # @test z[1:2,1,3:4] == _sliceof_util(VariableRefType, z, 1:2, 1, 3:4)
    # @test_throws BoundsError z[1:2,1,1:3]

    #@test w[:] == vec(_sliceof_util(VariableRefType, w, 1:3, -1:2, [:red,"blue"]))
    @test w[:, :, :] == w
    @test w[1, :, "blue"].data ==
          _sliceof_util(VariableRefType, w, 1, -1:2, ["blue"])
    @test w[1, :, :red].data ==
          _sliceof_util(VariableRefType, w, 1, -1:2, [:red])
    @test_throws KeyError w[1, :, "green"]
    # @test w[1:2,:,"blue"] == _sliceof_util(VariableRefType, w, 1:2, -1:2, ["blue"])
    # @test_throws ErrorException w[1:2,:,[:red,"blue"]]
    return
end

function test_extension_variable_end_indexing(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[0:2, 1:4])
    @variable(model, z[0:2])
    @test x[end, 1] == x[2, 1]
    @test x[0, end-1] == x[0, 3]
    @test z[end] == z[2]
    # TODO: It is redirected to x[11] as it is the 11th element but linear
    #       indexing is not supported
    @test_throws BoundsError x[end-1]
    return
end

function test_extension_variable_unsigned_index(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    # Tests unsigned int can be used to construct index set (JuMP issue #857).
    model = ModelType()
    t = UInt(4)
    @variable(model, x[1:t])
    @test 4 == @inferred num_variables(model)
    return
end

function test_extension_variable_symmetric(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2, 1:2], Symmetric)
    @test x isa LinearAlgebra.Symmetric
    @test x[1, 2] === x[2, 1]
    @test model[:x] === x
    y = @variable(model, [1:2, 1:2], Symmetric)
    @test y isa LinearAlgebra.Symmetric
    @test y[1, 2] === y[2, 1]
    return
end

function test_extension_variable_skewsymmetric(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2, 1:2] in SkewSymmetricMatrixSpace())
    @test x[1, 2] == -x[2, 1]
    @test iszero(x[1, 1])
    @test iszero(x[2, 2])
    @test model[:x] === x
    @test num_variables(model) == 1
    @variable(model, z[1:3, 1:3] in SkewSymmetricMatrixSpace())
    @test z[1, 2] == -z[2, 1]
    @test z[1, 3] == -z[3, 1]
    @test z[2, 3] == -z[3, 2]
    @test iszero(z[1, 1])
    @test iszero(z[2, 2])
    @test iszero(z[3, 3])
    @test model[:z] === z
    @test num_variables(model) == 4
    y = @variable(model, [1:3, 1:3] in SkewSymmetricMatrixSpace())
    @test y[1, 2] == -y[2, 1]
    @test y[2, 3] == -y[3, 2]
    @test iszero(y[3, 3])
    @test num_variables(model) == 7
    return
end

function test_extension_variables_constrained_on_creation_errors(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, x[1:2] in SecondOrderCone(), set = PSDCone())`: " *
            "Cannot use set keyword because the variable is already " *
            "constrained to `$(Expr(:escape, :(SecondOrderCone())))`.",
        ),
        @variable(model, x[1:2] in SecondOrderCone(), set = PSDCone()),
    )
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, x[1:2] in SecondOrderCone(), PSD)`: " *
            "Cannot pass `PSD` as a positional argument because the variable " *
            "is already constrained to `$(Expr(:escape, :(SecondOrderCone())))`.",
        ),
        @variable(model, x[1:2] in SecondOrderCone(), PSD),
    )
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, x[1:2, 1:2], PSD, Symmetric)`: " *
            "Cannot pass `Symmetric` as a positional argument because the " *
            "variable is already constrained to `$(PSDCone())`.",
        ),
        @variable(model, x[1:2, 1:2], PSD, Symmetric),
    )
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, x[1:2], set = SecondOrderCone(), set = PSDCone())`: " *
            "`set` keyword argument was given 2 times.",
        ),
        @variable(model, x[1:2], set = SecondOrderCone(), set = PSDCone()),
    )
    return
end

function test_extension_variables_constrained_on_creation(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2] in SecondOrderCone())
    @test num_constraints(model, typeof(x), MOI.SecondOrderCone) == 1
    @test name(x[1]) == "x[1]"
    @test name(x[2]) == "x[2]"

    @variable(model, [1:2] in SecondOrderCone())
    @test num_constraints(model, typeof(x), MOI.SecondOrderCone) == 2

    @variable(model, [1:3] ∈ MOI.SecondOrderCone(3))
    @test num_constraints(model, typeof(x), MOI.SecondOrderCone) == 3

    z = @variable(model, z ∈ MOI.Semiinteger(1.0, 2.0))
    @test num_constraints(model, typeof(z), MOI.Semiinteger{Float64}) == 1

    @variable(model, set = MOI.Semiinteger(1.0, 2.0))
    @test num_constraints(model, typeof(z), MOI.Semiinteger{Float64}) == 2

    X = @variable(model, [1:3, 1:3] in PSDCone())
    @test X isa LinearAlgebra.Symmetric
    @test num_constraints(
        model,
        typeof(x),
        MOI.PositiveSemidefiniteConeTriangle,
    ) == 1
    return
end

function test_extension_batch_delete_variables(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:3] >= 1)
    @objective(model, Min, sum([1, 2, 3] .* x))
    @test all(is_valid.(model, x))
    delete(model, x[[1, 3]])
    @test all((!is_valid).(model, x[[1, 3]]))
    @test is_valid(model, x[2])
    second_model = ModelType()
    @test_throws Exception delete(second_model, x[2])
    @test_throws Exception delete(second_model, x[[1, 3]])
    return
end

function test_all_variable()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    @test [x, y] == @inferred all_variables(model)
    return
end

function test_macro_variables()
    model = Model()
    @variables model begin
        0 ≤ x[i = 1:2] ≤ i
        y ≥ 2, Int, (start = 0.7)
        z ≤ 3, (start = 10)
        q, (Bin, start = 0.5)
    end

    @test "x[1]" == @inferred name(x[1])
    @test 0 == @inferred lower_bound(x[1])
    @test 1 == @inferred upper_bound(x[1])
    @test !is_binary(x[1])
    @test !is_integer(x[1])
    @test start_value(x[1]) === nothing

    @test "x[2]" == @inferred name(x[2])
    @test 0 == @inferred lower_bound(x[2])
    @test 2 == @inferred upper_bound(x[2])
    @test !is_binary(x[2])
    @test !is_integer(x[2])
    @test start_value(x[2]) === nothing

    @test "y" == @inferred name(y)
    @test 2 == @inferred lower_bound(y)
    @test !has_upper_bound(y)
    @test !is_binary(y)
    @test is_integer(y)
    @test start_value(y) === 0.7

    @test "z" == @inferred name(z)
    @test !has_lower_bound(z)
    @test 3 == @inferred upper_bound(z)
    @test !is_binary(z)
    @test !is_integer(z)
    @test start_value(z) === 10.0

    @test "q" == @inferred name(q)
    @test !has_lower_bound(q)
    @test !has_upper_bound(q)
    @test is_binary(q)
    @test !is_integer(q)
    @test start_value(q) === 0.5
    return
end

function test_dual_variable()
    model = Model()
    @variable(model, x == 0)
    exception = ErrorException(
        "To query the dual variables associated with a variable bound, first " *
        "obtain a constraint reference using one of `UpperBoundRef`, `LowerBoundRef`, " *
        "or `FixRef`, and then call `dual` on the returned constraint reference.\nFor " *
        "example, if `x <= 1`, instead of `dual(x)`, call `dual(UpperBoundRef(x))`.",
    )
    @test_throws exception dual(x)
    return
end

function test_value_containers()
    model = Model()
    @variable(model, x[1:2])
    exception = ErrorException(
        "`JuMP.value` is not defined for collections of JuMP types. Use " *
        "Julia's broadcast syntax instead: `JuMP.value.(x)`.",
    )
    @test_throws exception value(x)
    return
end

function test_get_variable_coefficient()
    m = Model()
    x = @variable(m, x)
    y = @variable(m, y)
    @test coefficient(x, x) == 1.0
    @test coefficient(x, y) == 0.0
    @test coefficient(x, x, x) == 0.0
    @test coefficient(x, y, x) == coefficient(x, x, y) == 0.0
    @test coefficient(x, y, y) == 0.0
    return
end

function _mock_reduced_cost_util(
    obj_sense::OptimizationSense,
    #obj_value::Float64,
    var_obj_coeff::Float64,
    #var_value,
    var_bound_type::Symbol,
    var_bounds_dual = nothing,
    has_duals::Bool = var_bounds_dual !== nothing,
)
    mockoptimizer =
        MOIU.MockOptimizer(MOIU.Model{Float64}(); eval_objective_value = false)
    m = direct_model(mockoptimizer)
    if var_bound_type === :lower
        @variable(m, x >= 0)
        if has_duals
            @assert isa(var_bounds_dual, Float64)
            has_duals && MOI.set(
                mockoptimizer,
                MOI.ConstraintDual(),
                optimizer_index(LowerBoundRef(x)),
                var_bounds_dual,
            )
        end
    elseif var_bound_type === :upper
        @variable(m, x <= 10)
        if has_duals
            @assert isa(var_bounds_dual, Float64)
            has_duals && MOI.set(
                mockoptimizer,
                MOI.ConstraintDual(),
                optimizer_index(UpperBoundRef(x)),
                var_bounds_dual,
            )
        end
    elseif var_bound_type === :fixed
        @variable(m, x == 10)
        if has_duals
            @assert isa(var_bounds_dual, Float64)
            MOI.set(
                mockoptimizer,
                MOI.ConstraintDual(),
                optimizer_index(FixRef(x)),
                var_bounds_dual,
            )
        end
    elseif var_bound_type === :both
        @variable(m, 0 <= x <= 10)
        if has_duals
            @assert length(var_bounds_dual) == 2
            @assert eltype(var_bounds_dual) == Float64
            lb_dual, ub_dual = var_bounds_dual
            MOI.set(
                mockoptimizer,
                MOI.ConstraintDual(),
                optimizer_index(LowerBoundRef(x)),
                lb_dual,
            )
            MOI.set(
                mockoptimizer,
                MOI.ConstraintDual(),
                optimizer_index(UpperBoundRef(x)),
                ub_dual,
            )
        end
    elseif var_bound_type === :none
        @variable(m, x)
        @assert var_bounds_dual === nothing
    else
        error("unrecognized bound type")
    end
    @objective(m, obj_sense, var_obj_coeff * x)
    if has_duals
        MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
        MOI.set(mockoptimizer, MOI.ResultCount(), 1)
        MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    end
    optimize!(m)
    return x
end

function test_reduced_cost()
    Min = MIN_SENSE
    Max = MAX_SENSE
    # The method should always fail if duals are not available.
    x = _mock_reduced_cost_util(Min, 1.0, :none)
    @test_throws ErrorException reduced_cost(x)
    x = _mock_reduced_cost_util(Min, 1.0, :fixed)
    @test_throws ErrorException reduced_cost(x)
    x = _mock_reduced_cost_util(Min, 1.0, :lower)
    @test_throws ErrorException reduced_cost(x)
    x = _mock_reduced_cost_util(Min, 1.0, :upper)
    @test_throws ErrorException reduced_cost(x)
    x = _mock_reduced_cost_util(Min, 1.0, :both)
    @test_throws ErrorException reduced_cost(x)
    # My reimplementation of the tests suggested by @odow.
    # Note that the floating point values are compared by equality because
    # there is no risk of the solver messing this up (mocks are being used).
    # First the fixed variable tests.
    x = _mock_reduced_cost_util(Min, 1.0, :none, nothing, true) # free var
    @test reduced_cost(x) == 0.0
    x = _mock_reduced_cost_util(Min, 1.0, :fixed, 1.0) # min x, x == 10
    @test reduced_cost(x) == 1.0
    x = _mock_reduced_cost_util(Max, 1.0, :fixed, -1.0) # max x, x == 10
    @test reduced_cost(x) == 1.0
    x = _mock_reduced_cost_util(Min, -1.0, :fixed, -1.0) # min -x, x == 10
    @test reduced_cost(x) == -1.0
    x = _mock_reduced_cost_util(Max, -1.0, :fixed, 1.0) # max -x, x == 10
    @test reduced_cost(x) == -1.0
    # Then the double bounded variables.
    #x = _mock_reduced_cost_util(Min, 1.0, :both, (0.0, 1.0)) # min x, 0 <= x <= 10
    #@test reduced_cost(x) == 1.0
    #x = _mock_reduced_cost_util(Max, 1.0, :both, (-1.0, 0.0)) # max x, 0 <= x <= 10
    #@test reduced_cost(x) == 1.0
    #x = _mock_reduced_cost_util(Min, -1.0, :both, (-1.0, 0.0)) # min -x, 0 <= x <= 10
    #@test reduced_cost(x) == -1.0
    #x = _mock_reduced_cost_util(Max, -1.0, :both, (0.0, 1.0)) # max -x, 0 <= x <= 10
    #@test reduced_cost(x) == -1.0
    # Test for a single upper bound and a single lower bound.
    x = _mock_reduced_cost_util(Min, 1.0, :lower, 1.0) # min x, 0 <= x
    @test reduced_cost(x) == 1.0
    x = _mock_reduced_cost_util(Max, 1.0, :upper, 1.0) # max x, x <= 10
    @test reduced_cost(x) == 1.0
    return
end

function test_value()
    @test value(1) === 1
    @test value(1.0) === 1.0
    @test value(JuMP._MA.Zero()) === 0.0
    return
end

function test_value_var()
    model = Model()
    @variable(model, x[1:2])
    vals = Dict(x[1] => 1.0, x[2] => 2.0)
    f = vidx -> vals[vidx]
    @test value(f, x[1]) === 1.0
    @test value(f, x[2]) === 2.0
    return
end

function test_relax_integrality()
    model = Model()
    @variable(model, x, Bin)
    @variable(model, -1 <= y <= 2, Bin)
    @variable(model, 0.1 <= z <= 0.6, Bin)
    @variable(model, a, Int)
    @variable(model, -1 <= b <= 2, Int)
    unrelax = relax_integrality(model)

    @test !is_binary(x)
    @test lower_bound(x) == 0.0
    @test upper_bound(x) == 1.0

    @test !is_binary(y)
    @test lower_bound(y) == 0.0
    @test upper_bound(y) == 1.0

    @test !is_binary(z)
    @test lower_bound(z) == 0.1
    @test upper_bound(z) == 0.6

    @test !is_integer(a)
    @test !has_lower_bound(a)
    @test !has_upper_bound(a)

    @test !is_integer(b)
    @test lower_bound(b) == -1.0
    @test upper_bound(b) == 2.0

    unrelax()

    @test is_binary(x)
    @test !has_lower_bound(x)
    @test !has_upper_bound(x)

    @test is_binary(y)
    @test lower_bound(y) == -1.0
    @test upper_bound(y) == 2.0

    @test is_binary(z)
    @test lower_bound(z) == 0.1
    @test upper_bound(z) == 0.6

    @test is_integer(a)
    @test !has_lower_bound(a)
    @test !has_upper_bound(a)

    @test is_integer(b)
    @test lower_bound(b) == -1.0
    @test upper_bound(b) == 2.0

    fix(x, 1)
    unrelax = relax_integrality(model)
    @test !is_binary(x)
    @test is_fixed(x)
    unrelax()
    @test is_binary(x)
    @test is_fixed(x)

    fix(x, 0)
    unrelax = relax_integrality(model)
    @test !is_binary(x)
    @test is_fixed(x)
    unrelax()
    @test is_binary(x)
    @test is_fixed(x)

    fix(a, 1)
    unrelax = relax_integrality(model)
    @test !is_integer(a)
    @test is_fixed(a)
    unrelax()
    @test is_integer(a)
    @test is_fixed(a)
    return
end

function test_relax_integrality_error_cases()
    model = Model()
    @variable(model, x)
    @constraint(model, x in MOI.Semicontinuous(1.0, 2.0))
    err = ErrorException(
        "Support for relaxing semicontinuous constraints " *
        "is not yet implemented.",
    )
    @test_throws err relax_integrality(model)

    model = Model()
    @variable(model, x)
    @constraint(model, x in MOI.Semiinteger(1.0, 2.0))
    err = ErrorException(
        "Support for relaxing semi-integer constraints " *
        "is not yet implemented.",
    )
    @test_throws err relax_integrality(model)

    model = Model()
    @variable(model, x, Bin)
    fix(x, 2)
    err = ErrorException(
        "The model has no valid relaxation: binary variable " *
        "fixed out of bounds.",
    )
    @test_throws err relax_integrality(model)

    model = Model()
    @variable(model, x, Bin)
    fix(x, -1)
    err = ErrorException(
        "The model has no valid relaxation: binary variable " *
        "fixed out of bounds.",
    )
    @test_throws err relax_integrality(model)
    return
end

function test_unknown_size_dense()
    model = Model()
    f = Iterators.filter(k -> isodd(k), 1:10)
    @variable(model, x[f])
    @test length(x) == 5
    return
end

function test_unknown_size_sparse()
    model = Model()
    f = Iterators.filter(k -> isodd(k), 1:10)
    @variable(model, x[i = f; i < 5])
    @test length(x) == 2
    return
end

function test_start_value()
    model = Model()
    @variable(model, x)
    @test start_value(x) === nothing
    set_start_value(x, 1.0)
    @test start_value(x) == 1.0
    set_start_value(x, nothing)
    @test start_value(x) === nothing
    set_start_value(x, 1)
    @test start_value(x) == 1.0
    return
end

function test_inf_lower_bound()
    for y in [-Inf, Inf, NaN]
        model = Model()
        @variable(model, x >= y)
        @test !has_lower_bound(x)
        @test_throws(
            ErrorException(
                "Unable to set lower bound to $y. To remove the bound, use " *
                "`delete_lower_bound`.",
            ),
            set_lower_bound(x, y),
        )
    end
    return
end

function test_inf_upper_bound()
    for y in [-Inf, Inf, NaN]
        model = Model()
        @variable(model, x <= y)
        @test !has_upper_bound(x)
        @test_throws(
            ErrorException(
                "Unable to set upper bound to $y. To remove the bound, use " *
                "`delete_upper_bound`.",
            ),
            set_upper_bound(x, y),
        )
    end
    return
end

function test_inf_fixed()
    for y in [-Inf, Inf, NaN]
        model = Model()
        @test_throws(
            ErrorException("Unable to fix variable to $y"),
            @variable(model, x == y),
        )
        @variable(model, x)
        @test_throws(ErrorException("Unable to fix variable to $y"), fix(x, y))
    end
    return
end

struct _UnsupportedVariableName <: MOI.AbstractOptimizer end
MOI.add_variable(::_UnsupportedVariableName) = MOI.VariableIndex(1)
MOI.is_empty(::_UnsupportedVariableName) = true

function test_unsupported_VariableName()
    model = direct_model(_UnsupportedVariableName())
    @variable(model, x)
    @test x isa VariableRef
    @test name(x) == ""
    return
end

function test_error_messages()
    model = Model()
    @variable(model, x)
    err = try
        x >= 1
    catch err
        err
    end
    function f(s)
        return ErrorException(
            replace(replace(err.msg, ">= 1" => "$(s) 1"), "`>=`" => "`$(s)`"),
        )
    end
    @test_throws err 1 >= x
    @test_throws f("<=") x <= 1
    @test_throws f("<=") 1 <= x
    @test_throws f(">") x > 1
    @test_throws f(">") 1 > x
    @test_throws f("<") x < 1
    @test_throws f("<") 1 < x
    return
end

function test_extension_rational_inf_bounds(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    u = Rational{Int}(Inf)
    @variable(model, -u <= x <= u)
    @test has_lower_bound(x) == false
    @test has_upper_bound(x) == false
    return
end

function test_ConstraintRef()
    model = Model()
    @variable(model, x >= 0)
    c = LowerBoundRef(x)
    @constraint(model, c2, x >= 0)
    @test VariableRef(c) == x
    @test_throws(MethodError, VariableRef(c2))
    return
end

function test_extension_start_value_nothing(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x, start = nothing)
    @test start_value(x) === nothing
    return
end

function test_VariableRef()
    model = Model()
    x = VariableRef(model)
    @test x isa VariableRef
    @test name(x) == ""
    @test num_variables(model) == 1
    return
end

function test_VariableIndex_VariableRef()
    model = Model()
    @variable(model, x)
    @test MOI.VariableIndex(x) === index(x)
    return
end

function test_VariableIndex_VariableRef_fix_with_upper_bound()
    model = Model()
    @variable(model, x <= 2)
    fix(x, 1.0; force = true)
    @test is_fixed(x)
    @test fix_value(x) == 1.0
    @test !has_upper_bound(x)
    return
end

function test_complex_variable()
    model = Model()
    @variable(
        model,
        x in ComplexPlane(),
        start = 5 + 6im,
        lower_bound = 1 + 2im,
        upper_bound = 3 + 4im
    )
    xr = first(x.terms).first
    xi = collect(x.terms)[2].first
    @test lower_bound(xr) == 1
    @test upper_bound(xr) == 3
    @test start_value(xr) == 5
    @test lower_bound(xi) == 2
    @test upper_bound(xi) == 4
    @test start_value(xi) == 6
    @test num_variables(model) == 2
    v = all_variables(model)
    @test xr == v[1]
    @test name(v[1]) == "real(x)"
    @test xi == v[2]
    @test name(v[2]) == "imag(x)"
    @test x == v[1] + v[2] * im
    @test conj(x) == v[1] - v[2] * im
    return
end

function test_complex_variable()
    model = Model()
    @variable(model, x[1:2] in ComplexPlane())
    @test x[1] isa GenericAffExpr{ComplexF64,VariableRef}
    @test x[2] isa GenericAffExpr{ComplexF64,VariableRef}
    @test num_variables(model) == 4
    return
end

function test_extension_complex_variable_errors(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @test_throws ErrorException @variable(model, x in ComplexPlane(), Int)
    @test_throws ErrorException @variable(model, x in ComplexPlane(), Bin)
    return
end

function _test_Hermitian(model, H)
    @test H isa LinearAlgebra.Hermitian
    Q = parent(H)
    @test num_variables(model) == 4
    v = all_variables(model)
    _test_variable_name_util(v[1], "real(H[1,1])")
    @test Q[1, 1] == 1v[1]
    _test_variable_name_util(v[2], "real(H[1,2])")
    _test_variable_name_util(v[4], "imag(H[1,2])")
    @test Q[1, 2] == v[2] + v[4] * im
    _test_variable_name_util(v[3], "real(H[2,2])")
    @test Q[2, 2] == 1v[3]
    @test Q[2, 1] == conj(Q[1, 2])
    @test H[2, 1] == conj(H[1, 2])
end

function test_extension_variable_Hermitian_tag(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, H[1:2, 1:2], Hermitian)
    _test_Hermitian(model, H)
    @test num_constraints(model; count_variable_in_set_constraints = true) == 0
    return
end

function test_extension_variable_Hermitian(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, H[1:2, 1:2] in HermitianMatrixSpace())
    _test_Hermitian(model, H)
    @test num_constraints(model; count_variable_in_set_constraints = true) == 0
    return
end

function test_extension_Hermitian_PSD(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, H[1:2, 1:2] in HermitianPSDCone())
    _test_Hermitian(model, H)
    con_refs = all_constraints(
        model,
        Vector{VariableRefType},
        MOI.HermitianPositiveSemidefiniteConeTriangle,
    )
    @test length(con_refs) == 1
    con = constraint_object(con_refs[1])
    @test jump_function(con) == all_variables(model)
    @test moi_set(con) == MOI.HermitianPositiveSemidefiniteConeTriangle(2)
    return
end

function _test_Hermitian_errors(model, set)
    @test_throws ErrorException @variable(
        model,
        H[i = 1:2, j = 1:2] in set,
        lower_bound = (i + j) * im
    )
    @test_throws ErrorException @variable(
        model,
        H[i = 1:2, j = 1:2] in set,
        upper_bound = (i + j) * im
    )
    @test_throws ErrorException @variable(
        model,
        H[i = 1:2, j = 1:2] in set,
        start = (i + j) * im
    )
    @test_throws ErrorException @variable(
        model,
        H[i = 1:2, j = 1:2] in set,
        integer = i > j
    )
    @test_throws ErrorException @variable(
        model,
        H[i = 1:2, j = 1:2] in set,
        Bin
    )
    @test_throws ErrorException @variable(
        model,
        H[i = 1:2, j = 1:2] in set,
        Int,
    )
    @test_throws ErrorException @variable(
        model,
        H[i = 1:2, j = 1:2] in set,
        integer = i != j,
    )
    return
end

function test_extension_Hermitian_errors(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    _test_Hermitian_errors(ModelType(), HermitianMatrixSpace())
    return
end

function test_extension_Hermitian_PSD_errors(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    _test_Hermitian_errors(ModelType(), HermitianPSDCone())
    return
end

function test_Hermitian_PSD_keyword()
    model = Model()
    @test_throws ErrorException @variable(
        model,
        H[i = 1:2, j = 1:2] in HermitianPSDCone(),
        integer = i != j,
    )
    @variable(
        model,
        H[i = 1:2, j = 1:2] in HermitianPSDCone(),
        lower_bound = (i + j) + (i - j) * im,
        upper_bound = i * j + (j - i) * im
    )
    v = all_variables(model)
    for i in 1:2, j in 1:2
        @test value(lower_bound, H[i, j]) == (i + j) + (i - j) * im
    end
    for i in 1:2, j in 1:2
        @test value(upper_bound, H[i, j]) == (i * j) + (j - i) * im
    end
    @variable(
        model,
        T[i = 1:2, j = 1:2] in HermitianPSDCone(),
        start = (i + j) + (j - i) * im
    )
    for i in 1:2, j in 1:2
        @test value(start_value, T[i, j]) == (i + j) + (j - i) * im
    end
    return
end

function test_Hermitian_PSD_anon()
    model = Model()
    y = @variable(model, [1:2, 1:2] in HermitianPSDCone())
    @test y isa LinearAlgebra.Hermitian
    x = parent(y)
    @test sprint(show, x[1, 1]) == "_[1]"
    @test sprint(show, y[1, 1]) == "_[1]"
    @test sprint(show, x[1, 2]) == "_[2] + _[4] im"
    @test sprint(show, y[1, 2]) == "_[2] + _[4] im"
    @test sprint(show, x[2, 1]) == "_[2] - _[4] im"
    @test sprint(show, y[2, 1]) == "_[2] - _[4] im"
    @test sprint(show, x[2, 2]) == "_[3]"
    @test sprint(show, y[2, 2]) == "_[3]"
    return
end

function test_fix_discrete_variables()
    le = JuMP._math_symbol(MIME("text/plain"), :leq)
    ge = JuMP._math_symbol(MIME("text/plain"), :geq)
    eq = JuMP._math_symbol(MIME("text/plain"), :eq)
    model = Model()
    @variable(model, x, Bin, start = 1)
    @variable(model, 1 <= y <= 10, Int, start = 2)
    @objective(model, Min, x + y)
    @test sprint(print, model) ==
          "Min x + y\nSubject to\n y $ge 1\n y $le 10\n y integer\n x binary\n"
    undo_relax = fix_discrete_variables(start_value, model)
    @test sprint(print, model) == "Min x + y\nSubject to\n x $eq 1\n y $eq 2\n"
    undo_relax()
    @test sprint(print, model) ==
          "Min x + y\nSubject to\n y $ge 1\n y $le 10\n y integer\n x binary\n"
    return
end

function test_fix_discrete_variables_value()
    model = Model()
    @variable(model, x, Bin, start = 1)
    @variable(model, 1 <= y <= 10, Int, start = 2)
    @objective(model, Min, x + y)
    @test_throws OptimizeNotCalled fix_discrete_variables(model)
    return
end

function test_string_names_set_on_creation_constrained_on_creation()
    model = Model()
    set_string_names_on_creation(model, false)
    @variable(model, x in MOI.ZeroOne())
    @test x isa VariableRef
    @test isempty(name(x))
    @variable(model, y[1:2] in MOI.ZeroOne())
    @test y[1] isa VariableRef
    @test all(isempty, name.(y))
    return
end

function test_dependent_set_variable_macro()
    for C in (Array, Containers.DenseAxisArray, Containers.SparseAxisArray)
        model = Model()
        @variable(model, x[i = 1:3] in MOI.EqualTo(1.0 * i), container = C)
        @test all(is_fixed.(x))
        @test [fix_value(x[i]) for i in 1:3] == [1.0, 2.0, 3.0]
    end
    for C in (Array, Containers.DenseAxisArray, Containers.SparseAxisArray)
        model = Model()
        @variable(model, x[i = 1:3] in MOI.ZeroOne(), container = C)
        @test all(is_binary.(x))
    end
    for C in (Array, Containers.DenseAxisArray, Containers.SparseAxisArray)
        model = Model()
        @variable(model, x[i = 1:3] in MOI.GreaterThan(sqrt(i)), container = C)
        @test all(has_lower_bound.(x))
        @test [lower_bound(x[i]) for i in 1:3] == sqrt.([1.0, 2.0, 3.0])
    end
    return
end

end  # module TestVariable
