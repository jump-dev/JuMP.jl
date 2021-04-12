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
        MOI.SingleVariable,
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
            moi_add_constraint,
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

    # Nonlinear interface
    Base.precompile(
        Tuple{
            typeof(MOI.eval_hessian_lagrangian),
            NLPEvaluator,
            SubArray{Float64,1,Vector{Float64},Tuple{UnitRange{Int64}},true},
            Vector{Float64},
            Float64,
            SubArray{Float64,1,Vector{Float64},Tuple{UnitRange{Int64}},true},
        },
    )   # time: 0.83568186
    Base.precompile(
        Tuple{
            typeof(MOI.eval_objective_gradient),
            NLPEvaluator,
            Vector{Float64},
            Vector{Float64},
        },
    )   # time: 0.50052494
    Base.precompile(Tuple{typeof(MOI.initialize),NLPEvaluator,Vector{Symbol}})   # time: 0.34357876
    Base.precompile(Tuple{typeof(MOI.features_available),NLPEvaluator})   # time: 0.014501549
    Base.precompile(
        Tuple{
            typeof(MOI.eval_constraint_jacobian),
            NLPEvaluator,
            SubArray{Float64,1,Vector{Float64},Tuple{UnitRange{Int64}},true},
            Vector{Float64},
        },
    )   # time: 0.00857948
    Base.precompile(
        Tuple{
            typeof(MOI.eval_constraint),
            NLPEvaluator,
            SubArray{Float64,1,Vector{Float64},Tuple{UnitRange{Int64}},true},
            Vector{Float64},
        },
    )   # time: 0.003288632
    Base.precompile(Tuple{typeof(MOI.jacobian_structure),NLPEvaluator})   # time: 0.002370995
    Base.precompile(
        Tuple{typeof(MOI.eval_objective),NLPEvaluator,Vector{Float64}},
    )   # time: 0.002276394
    return
end
