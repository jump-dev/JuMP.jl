let
    # Intially
    # max 1.1x + 1.0y
    # st     x +    y <= 3
    #     0 <= x <= 3
    #     1 <= y <= 3
    # x* = 2, y* = 1
    m = Model()
    @defVar(m, 0 <= x <= 3)
    @defVar(m, 1 <= y <= 3)
    @setObjective(m, :Max, 1.1x + 1.0y)
    maincon = @addConstraint(m, x + y <= 3)
    solve(m)
    @test_approx_eq_eps getValue(x) 2.0 1e-6
    @test_approx_eq_eps getValue(y) 1.0 1e-6

    # Test adding a variable
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 3
    #     0 <= x <= 3
    #     1 <= y <= 3
    #     0 <= z <= 5
    # x* = 0, y* = 1, z* = 2
    @defVar(m, 0 <= z <= 5, 100.0, [maincon], [1.0])
    @test m.internalModelLoaded
    solve(m)
    @test_approx_eq_eps getValue(x) 0.0 1e-6
    @test_approx_eq_eps getValue(y) 1.0 1e-6
    @test_approx_eq_eps getValue(z) 2.0 1e-6
    @test_approx_eq_eps getObjectiveValue(m) 201.0 1e-6

    # Test changing bounds
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 3
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 2
    # x* = 1, y* = 0, z* = 2
    setLower(y, 0.0)
    setUpper(z, 2.0)
    solve(m)
    @test_approx_eq_eps getValue(x) 1.0 1e-6
    @test_approx_eq_eps getValue(y) 0.0 1e-6
    @test_approx_eq_eps getValue(z) 2.0 1e-6
    m.internalModelLoaded = false

    # Test changing problem type
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 3
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 1.5, Integer
    # x* = 2, y* = 0, z* = 1
    setUpper(z, 1.5)
    m.colCat[3] = JuMP.INTEGER
    solve(m)
    @test_approx_eq_eps getValue(x) 2.0 1e-6
    @test_approx_eq_eps getValue(y) 0.0 1e-6
    @test_approx_eq_eps getValue(z) 1.0 1e-6

    # Test changing constraint bound (<= constraint)
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 2
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 1.5, Integer
    # x* = 1, y* = 0, z* = 1
    chgConstrRHS(maincon, 2.0)
    solve(m)
    @test_approx_eq_eps getValue(x) 1.0 1e-6
    @test_approx_eq_eps getValue(y) 0.0 1e-6
    @test_approx_eq_eps getValue(z) 1.0 1e-6

    # Test adding a constraint
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 2
    #        x        +      z <= 0
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 1.5, Integer
    # x* = 0, y* = 2, z* = 0
    xz0ref = @addConstraint(m, x + z <=0)
    solve(m)
    @test_approx_eq_eps getValue(x) 0.0 1e-6
    @test_approx_eq_eps getValue(y) 2.0 1e-6
    @test_approx_eq_eps getValue(z) 0.0 1e-6

    # Test adding a range constraint and modifying it
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 2
    #        x        +      z <= 0
    #                 -10 <= x <= 10
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 1.5, Integer
    rangeref = @addConstraint(m, -10 <= x <= 10)
    @test_throws chgConstrRHS(rangeref, 11)

    # Test changing constraint bound (>= constraint)
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 2
    #        x        +      z <= 1
    #        x +    y          >= 1
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 1.5, Integer
    # x* = 1, y* = 0, z* = 1
    chgConstrRHS(xz0ref, 2.0)
    xyg0ref = @addConstraint(m, x + y >= 0)
    solve(m)
    chgConstrRHS(xyg0ref, 1)
    solve(m)
    @test_approx_eq_eps getValue(x) 1.0 1e-6
    @test_approx_eq_eps getValue(y) 0.0 1e-6
    @test_approx_eq_eps getValue(z) 1.0 1e-6
end

let
    m = Model()
    @defVar(m, x >= 0)
    @defVar(m, y >= 0)
    @addConstraint(m, x + y == 1)
    @setObjective(m, Max, y)
    solve(m; load_model_only=true)
    @assert getInternalModel(m) != nothing
    @assert m.internalModelLoaded == true
    stat = solve(m)
    @test stat == :Optimal  
    @test getValue(x) == 0.
    @test getValue(y) == 1.
    @test getObjectiveValue(m) == 1.
    @test getDual(x) == -1.
    @test getDual(y) == 0.
end