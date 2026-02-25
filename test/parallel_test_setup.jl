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

include(joinpath(@__DIR__, "JuMPExtension.jl"))

function run_tests(m::Module)
    _startswith(k) = f -> startswith("$f", k)
    # Default tests
    for f in filter(_startswith("test_"), names(m; all = true))
        Test.@testset "$f" begin
            getfield(m, f)()
        end
    end
    # Test with {Float32}
    for f in filter(_startswith("test_extension_"), names(m; all = true))
        Test.@testset "$f::Float32" begin
            getfield(m, f)(
                JuMP.GenericModel{Float32},
                JuMP.GenericVariableRef{Float32},
            )
        end
    end
    # Test JuMP extension
    for f in filter(_startswith("test_extension_"), names(m; all = true))
        Test.@testset "$f::JuMPExtension" begin
            getfield(m, f)(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
        end
    end
    return
end
