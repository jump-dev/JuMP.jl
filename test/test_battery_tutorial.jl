#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestBatteryTutorial

import Test

function test_battery_tutorial()
    Base.include(
        Module(:BatteryTutorial),
        joinpath(
            @__DIR__,
            "..",
            "docs",
            "src",
            "tutorials",
            "applications",
            "battery.jl",
        ),
    )
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    Test.@testset "Battery tutorial" begin
        test_battery_tutorial()
    end
end

end # module
