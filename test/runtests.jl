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
import ParallelTestRunner
import Test

# It is important to test this _before_ calling `include_modules_to_test`
# because some of the tests introduce new ambiguities.
Test.@test isempty(Test.detect_ambiguities(JuMP; recursive = true))

is_test_file(f) = startswith(f, "test_") && endswith(f, ".jl")

testsuite = Dict{String,Expr}()
for (root, dirs, files) in walkdir(@__DIR__)
    for file in joinpath.(root, filter(is_test_file, files))
        testsuite[file] = :(run_tests(include($file)))
    end
end

const init_code = quote
    include(joinpath(@__DIR__, "parallel_test_setup.jl"))
end
ParallelTestRunner.runtests(JuMP, ARGS; testsuite, init_code)
