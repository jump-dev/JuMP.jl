# Pull packages in that are available
if Pkg.installed("Gurobi") != nothing  
    using Gurobi
end
if Pkg.installed("CPLEXLink") != nothing
    using CPLEXLink
end

# Only write test code once, use macros to copy-paste it
macro callback_test(solvername)
    if solvername == :Gurobi
        pkg_install  = :(Pkg.installed("Gurobi"))
        mod_creation = :(modQ = Model(solver=GurobiSolver()))
    elseif solvername == :CPLEXLink
        pkg_install  = :(Pkg.installed("CPLEXLink"))
        mod_creation = :(mod = Model(solver=CplexSolver()))
    end

    quote
        if $pkg_install != nothing  
            let
                $mod_creation
                @defVar(modQ, -2 <= x <= 2 )
                @defVar(modQ, -2 <= y <= 2 )

                @setObjective(modQ, Min, x - y )
                addConstraint(modQ, x + x*x + x*y + y*y <= 1 )

                status = solve(modQ)
                @test status == :Optimal
                @test_approx_eq_eps modQ.objVal -1-4/sqrt(3) 1e-6
                @test_approx_eq_eps (getValue(x) + getValue(y)) -1/3 1e-3
            end

            let
                $mod_creation

                @defVar(modQ, -2 <= x <= 2, Int )
                @defVar(modQ, -2 <= y <= 2, Int )

                @setObjective(modQ, Min, x - y )
                addConstraint(modQ, x + x*x + x*y + y*y <= 1 )

                status = solve(modQ)
                @test status == :Optimal
                @test_approx_eq_eps modQ.objVal -3 1e-6
                @test_approx_eq_eps (getValue(x) + getValue(y)) -1 1e-6
            end
        end  # pkg test
    end  # quote
end  # macro
