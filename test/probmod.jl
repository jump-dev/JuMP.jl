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

# Test adding a "decoupled" variable (#205)
let
    m = Model()
    @defVar(m, x >= 0)
    @setObjective(m, Min, x)
    solve(m)
    @defVar(m, y >= 0)
    @addConstraint(m, x + y == 1)
    @setObjective(m, Min, 2x+y)
    solve(m)
    @test_approx_eq getValue(x) 0.0
    @test_approx_eq getValue(y) 1.0
end

function methods_test(solvername, solverobj, supp)
    mod = Model(solver=solverobj)
    @defVar(mod, x >= 0)
    @addConstraint(mod, 2x == 2)
    solve(mod)
    internal_mod = getInternalModel(mod)
    for (it,(meth, args)) in enumerate(mpb_methods)
        if supp[it]
            @test applicable(meth, internal_mod, args...)
            @test method_exists(meth, map(typeof, tuple(internal_mod, args...)))
        end
    end
end

# test there were no regressions in applicable
const mpb_methods = [(MathProgBase.addquadconstr!, (Cint[1],Float64[1.0],Cint[1],Cint[1],Float64[1],'>',1.0)),
                     (MathProgBase.setquadobjterms!, (Cint[1], Cint[1], Float64[1.0])),
                     (MathProgBase.addconstr!,   ([1],[1.0],1.0,1.0)),
                     (MathProgBase.addsos1!,     ([1],[1.0])),
                     (MathProgBase.addsos2!,     ([1],[1.0])),
                     (MathProgBase.addvar!,      ([1],[1.0],1.0,1.0,1.0)),
                     (MathProgBase.setvarLB!,    ([1.0],)),
                     (MathProgBase.setvarUB!,    ([1.0],)),
                     (MathProgBase.setconstrLB!, ([1.0],)),
                     (MathProgBase.setconstrUB!, ([1.0],)),
                     (MathProgBase.setobj!,      ([1.0],)),
                     (MathProgBase.setsense!,    (:Min,)),
                     (MathProgBase.setvartype!,  (['C'],)),
                     (MathProgBase.getinfeasibilityray, ()),
                     (MathProgBase.getunboundedray, ()),
                     (MathProgBase.getreducedcosts, ()),
                     (MathProgBase.getconstrduals, ()),
                     (MathProgBase.setwarmstart!, ([1.0]))]

if Pkg.installed("Gurobi") != nothing
    using Gurobi
    supp = (true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true)
    methods_test("Gurobi", GurobiSolver(), supp)
end

if Pkg.installed("CPLEX") != nothing
    using CPLEX
    supp = (true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true)
    methods_test("CPLEX", CplexSolver(), supp)
end
if Pkg.installed("Clp") != nothing
    using Clp
    supp = (false,false,true,false,false,true,true,true,true,true,true,true,true,true,true,true,true,false)
    methods_test("Clp", ClpSolver(), supp)
end
if Pkg.installed("Cbc") != nothing
    # no-op, since Cbc doesn't support any
end
if Pkg.installed("GLPK") != nothing
    using GLPKMathProgInterface
    supp = (false,false,true,false,false,true,true,true,true,true,true,true,false,true,true,true,true,false)
    methods_test("GLPK", GLPKSolverLP(), supp)
    supp = (false,false,true,false,false, true,true,true,true,true,true,true,true,false,false,false,false,false)
    methods_test("GLPK", GLPKSolverMIP(), supp)
end
if Pkg.installed("Mosek") != nothing
    using Mosek
    supp = (true,true,true,false,false,true,true,true,true,true,true,true,true,true,true,true,true,false)
    methods_test("Mosek", MosekSolver(), supp)
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
