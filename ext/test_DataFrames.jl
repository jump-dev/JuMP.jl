#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestContainersDataFrames

using Test

using DataFrames
using JuMP

function test_dimension_data_variable()
    model = Model()
    @variable(model, x[i = 2:4, j = 1:2], container = DataFrame)
    @variable(model, y[i = 2:4, j = 1:2; isodd(i+j)], container = DataFrame)

    return
end

end
