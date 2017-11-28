@testset "Solvehook" begin
    m = Model()
    @test m.solvehook === nothing
    called = false
    function hook(m)
        called = true
    end
    JuMP.setsolvehook(m, hook)
    @test !called
    solve(m)
    @test called
end
