#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

import JuMP
import Test

# It is important to test this _before_ calling `include_modules_to_test`
# because some of the tests introduce new ambiguities.
if VERSION < v"1.8"
    # In Julia v1.6.x is one ambiguity with a method in StaticArrays that we
    # can't easily work-around without importing StaticArrays.
    Test.@test length(Test.detect_ambiguities(JuMP; recursive = true)) == 1
else
    Test.@test isempty(Test.detect_ambiguities(JuMP; recursive = true))
end

include("Kokako.jl")

const MODULES_TO_TEST = Kokako.include_modules_to_test(JuMP)

include(joinpath(@__DIR__, "JuMPExtension.jl"))

if isempty(ARGS)
    # JuMPExtension.jl also contains some tests.
    push!(MODULES_TO_TEST, "JuMPExtension.jl" => JuMPExtension)
end

Kokako.run_tests(MODULES_TO_TEST)

# `Float32` is cheaper to test than `BigFloat` and is enough
# for the purpose of testing that the types don't get promoted to `Float64`
Kokako.run_tests(
    MODULES_TO_TEST,
    JuMP.GenericModel{Float32},
    JuMP.GenericVariableRef{Float32};
    test_prefix = "test_extension_",
    include_names = Dict(
        "test_mutable_arithmetics.jl" => ["test_extension_promote_operation"],
    ),
)

Kokako.run_tests(
    MODULES_TO_TEST,
    JuMPExtension.MyModel,
    JuMPExtension.MyVariableRef;
    test_prefix = "test_extension_",
    include_names = Dict(
        "test_mutable_arithmetics.jl" => ["test_extension_promote_operation"],
    ),
)
