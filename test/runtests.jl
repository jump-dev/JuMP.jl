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
Test.@test isempty(Test.detect_ambiguities(JuMP))

# TODO(odow): there are still some ambiguities in Containers
# @test isempty(Test.detect_ambiguities(JuMP.Containers))

include("Kokako.jl")

const MODULES_TO_TEST = Kokako.include_modules_to_test(@__DIR__)

include(joinpath(@__DIR__, "JuMPExtension.jl"))

if isempty(ARGS)
    # JuMPExtension.jl also contains some tests.
    push!(MODULES_TO_TEST, "JuMPExtension.jl" => JuMPExtension)
end

Kokako.run_tests(MODULES_TO_TEST)

Kokako.run_tests(
    MODULES_TO_TEST,
    JuMPExtension.MyModel,
    JuMPExtension.MyVariableRef;
    test_prefix = "test_extension_",
    include_names = Dict(
        "test_mutable_arithmetics.jl" => ["test_extension_promote_operation"],
    ),
)
