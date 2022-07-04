# This precompile statements are a mix of hand-generated ones (e.g., to ensure
# coverage across the breadth of common constraint types, as well as some
# SnoopCompile ones that are expensive.)

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    T = Float64
    scalar_sets = (
        MOI.LessThan{T},
        MOI.GreaterThan{T},
        MOI.EqualTo{T},
        MOI.Interval{T},
        MOI.Integer,
        MOI.ZeroOne,
        MOI.Semiinteger{T},
        MOI.Semicontinuous{T},
    )
    scalar_functions = (
        MOI.VariableIndex,
        MOI.ScalarAffineFunction{T},
        MOI.ScalarQuadraticFunction{T},
    )
    vector_sets = (
        MOI.Nonnegatives,
        MOI.Nonpositives,
        MOI.Zeros,
        MOI.Reals,
        MOI.SecondOrderCone,
        MOI.RotatedSecondOrderCone,
        MOI.ExponentialCone,
        MOI.DualExponentialCone,
        MOI.PositiveSemidefiniteConeSquare,
        MOI.PositiveSemidefiniteConeTriangle,
    )
    vector_functions = (
        MOI.VectorOfVariables,
        MOI.VectorAffineFunction{T},
        MOI.VectorQuadraticFunction{T},
    )
    constraints = vcat(
        [(F, S) for F in scalar_functions for S in scalar_sets],
        [(F, S) for F in vector_functions for S in vector_sets],
    )
    for (F, S) in constraints
        Base.precompile(
            _moi_add_constraint,
            (
                MOIU.CachingOptimizer{
                    MOI.AbstractOptimizer,
                    MOIU.UniversalFallback{MOIU.Model{Float64}},
                },
                F,
                S,
            ),
        )
    end

    Base.precompile(Tuple{typeof(model_string),Type,Model})   # time: 0.25714603
    return
end
