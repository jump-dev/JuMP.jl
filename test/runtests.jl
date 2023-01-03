#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

using Test

# Run these tests only if the user does not specify a particular file to run.
if isempty(ARGS)
    # It is important to run these test _before_ any other packages are loaded.
    include("hygiene.jl")
    import JuMP
    @test isempty(Test.detect_ambiguities(JuMP))
    # TODO(odow): there are still some ambiguities in Containers
    # @test isempty(Test.detect_ambiguities(JuMP.Containers))
end

function run_tests(
    mod::Module,
    args...;
    test_prefix::String = "test_",
    include_names::Vector{String} = String[""],
    exclude_names::Vector{String} = String[],
)
    for name in names(mod; all = true)
        test_function, test_name = getfield(mod, name), string(name)
        if !(test_function isa Function)
            continue
        elseif !startswith(test_name, test_prefix)
            continue
        elseif !any(needle -> occursin(needle, test_name), include_names)
            continue
        elseif any(needle -> occursin(needle, test_name), exclude_names)
            continue
        end
        @testset "$name" begin
            test_function(args...)
        end
    end
    return
end

const FILES_TO_TEST = String[]

if isempty(ARGS)
    for (root, dirs, files) in walkdir(@__DIR__)
        for file in files
            if startswith(file, "test_") && endswith(file, ".jl")
                push!(FILES_TO_TEST, joinpath(root, file))
            end
        end
    end
else
    append!(FILES_TO_TEST, ARGS)
end

include(joinpath(@__DIR__, "JuMPExtension.jl"))

const EXTENSION_INCLUDES = Dict(
    "test_mutable_arithmetics.jl" => ["test_extension_promote_operation"],
)

for filename in FILES_TO_TEST
    file = last(split(filename, "/"))
    mod = include(filename)
    @testset "$file" begin
        run_tests(mod)  # Test default tests
        run_tests(  # Test JuMPExtension
            mod,
            JuMPExtension.MyModel,
            JuMPExtension.MyVariableRef;
            test_prefix = "test_extension_",
            include_names = get(EXTENSION_INCLUDES, file, [""]),
        )
    end
end
