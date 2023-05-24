#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestContainersDimensionalData

using Test

using DimensionalData
using JuMP
    
function test_dimension_data_integration()
    if !isdefined(Base, :get_extension)
        return
    end
    model = Model()
    @variable(model, x[i = 2:4, j = ["a", "b"]], container = DimArray)
    @test x isa DimArray
    @test_throws(
        MethodError,
        @variable(model, z[i = 1:3, j = i:2], container = DimArray),
    )
    @test_throws(
        ArgumentError,
        @variable(model, w[1:3, j = 1:2], container = DimArray),
    )
    return
end

end
