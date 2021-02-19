const __bodyfunction__ = Dict{Method,Any}()

# Find keyword "body functions" (the function that contains the body
# as written by the developer, called after all missing keyword-arguments
# have been assigned values), in a manner that doesn't depend on
# gensymmed names.
# `mnokw` is the method that gets called when you invoke it without
# supplying any keywords.
function __lookup_kwbody__(mnokw::Method)
    function getsym(arg)
        isa(arg, Symbol) && return arg
        @assert isa(arg, GlobalRef)
        return arg.name
    end

    f = get(__bodyfunction__, mnokw, nothing)
    if f === nothing
        fmod = mnokw.module
        # The lowered code for `mnokw` should look like
        #   %1 = mkw(kwvalues..., #self#, args...)
        #        return %1
        # where `mkw` is the name of the "active" keyword body-function.
        ast = Base.uncompressed_ast(mnokw)
        if isa(ast, Core.CodeInfo) && length(ast.code) >= 2
            callexpr = ast.code[end-1]
            if isa(callexpr, Expr) && callexpr.head == :call
                fsym = callexpr.args[1]
                if isa(fsym, Symbol)
                    f = getfield(fmod, fsym)
                elseif isa(fsym, GlobalRef)
                    if fsym.mod === Core && fsym.name === :_apply
                        f = getfield(mnokw.module, getsym(callexpr.args[2]))
                    elseif fsym.mod === Core && fsym.name === :_apply_iterate
                        f = getfield(mnokw.module, getsym(callexpr.args[3]))
                    else
                        f = getfield(fsym.mod, fsym.name)
                    end
                else
                    f = missing
                end
            else
                f = missing
            end
        else
            f = missing
        end
        __bodyfunction__[mnokw] = f
    end
    return f
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(MathOptInterface.eval_hessian_lagrangian),NLPEvaluator,SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true},Vector{Float64},Float64,SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}})   # time: 0.83568186
    Base.precompile(Tuple{typeof(MathOptInterface.eval_objective_gradient),NLPEvaluator,Vector{Float64},Vector{Float64}})   # time: 0.50052494
    Base.precompile(Tuple{typeof(MathOptInterface.initialize),NLPEvaluator,Vector{Symbol}})   # time: 0.34357876
    let fbody = try __lookup_kwbody__(which(Model, (Type,))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Bool,Base.Iterators.Pairs{Union{}, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},Type{Model},Type,))
        end
    end   # time: 0.26403505
    Base.precompile(Tuple{typeof(model_string),Type,Model})   # time: 0.25714603
    Base.precompile(Tuple{typeof(_moi_add_variable),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},Model,ScalarVariable{Int64, Float64, Float64, Float64},String})   # time: 0.21578874
    Base.precompile(Tuple{typeof(moi_add_constraint),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}})   # time: 0.19299112
    Base.precompile(Tuple{typeof(moi_add_constraint),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},MathOptInterface.ScalarQuadraticFunction{Float64},MathOptInterface.LessThan{Float64}})   # time: 0.12382254
    Base.precompile(Tuple{typeof(_hessian_slice),NLPEvaluator,_FunctionStorage,Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true},Float64,Int64,Vector{Float64},Type{Val{2}}})   # time: 0.123344906
    Base.precompile(Tuple{typeof(_process_NL_expr),Symbol,Expr})   # time: 0.08418905
    Base.precompile(Tuple{typeof(show),IOContext{IOBuffer},AffExpr})   # time: 0.058707304
    Base.precompile(Tuple{typeof(_moi_add_constrained_variables),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},Vector{ScalarVariable{Float64, Float64, Float64, Float64}},MathOptInterface.PositiveSemidefiniteConeTriangle,Vector{String}})   # time: 0.056060467
    Base.precompile(Tuple{typeof(_constraint_macro),Tuple{Symbol, Expr},Symbol,typeof(parse_constraint_expr),LineNumberNode})   # time: 0.037720144
    Base.precompile(Tuple{typeof(dual),ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.GreaterThan{Float64}}, ScalarShape}})   # time: 0.027809434
    isdefined(JuMP, Symbol("#108#109")) && Base.precompile(Tuple{getfield(JuMP, Symbol("#108#109")),Expr})   # time: 0.02720738
    Base.precompile(Tuple{Core.kwftype(typeof(register)),NamedTuple{(:autodiff,), Tuple{Bool}},typeof(register),Model,Symbol,Int64,Function})   # time: 0.026307618
    Base.precompile(Tuple{Core.kwftype(typeof(constraint_string)),NamedTuple{(:in_math_mode,), Tuple{Bool}},typeof(constraint_string),Type,ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarQuadraticFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}})   # time: 0.024550961
    Base.precompile(Tuple{typeof(moi_add_constraint),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.PositiveSemidefiniteConeSquare})   # time: 0.022663178
    Base.precompile(Tuple{typeof(_parse_NL_expr),Symbol,Expr,Symbol,Symbol,Symbol})   # time: 0.02125863
    Base.precompile(Tuple{Type{MathOptInterface.VectorOfVariables},Vector{VariableRef}})   # time: 0.01940772
    Base.precompile(Tuple{typeof(+),AffExpr,AffExpr})   # time: 0.017755583
    Base.precompile(Tuple{typeof(moi_add_constraint),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.SecondOrderCone})   # time: 0.016445328
    Base.precompile(Tuple{typeof(moi_add_constraint),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.EqualTo{Float64}})   # time: 0.016153162
    Base.precompile(Tuple{typeof(moi_add_constraint),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}})   # time: 0.015447589
    Base.precompile(Tuple{Core.kwftype(typeof(constraint_string)),NamedTuple{(:in_math_mode,), Tuple{Bool}},typeof(constraint_string),Type,ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape}})   # time: 0.015426918
    Base.precompile(Tuple{typeof(add_constraint),Model,VectorConstraint{AffExpr, MathOptInterface.PositiveSemidefiniteConeSquare, SquareMatrixShape},String})   # time: 0.014703301
    Base.precompile(Tuple{typeof(MathOptInterface.features_available),NLPEvaluator})   # time: 0.014501549
    Base.precompile(Tuple{Core.kwftype(typeof(constraint_string)),NamedTuple{(:in_math_mode,), Tuple{Bool}},typeof(constraint_string),Type,ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable, MathOptInterface.GreaterThan{Float64}}, ScalarShape}})   # time: 0.01252703
    Base.precompile(Tuple{Type{_FunctionStorage},Vector{NodeData},Vector{Float64},Int64,JuMP._Derivatives.Coloring.IndexedSet,Bool,Vector{_SubexpressionStorage},Vector{Int64},Vector{Linearity},Vector{Set{Tuple{Int64, Int64}}},Vector{Vector{Int64}},Dict{MathOptInterface.VariableIndex, Int64}})   # time: 0.012121839
    Base.precompile(Tuple{typeof(setindex!),Model,JuMP.Containers.DenseAxisArray{VariableRef, 3, Tuple{Vector{String}, Vector{Int64}, Base.OneTo{Int64}}, Tuple{Dict{String, Int64}, Dict{Int64, Int64}, Dict{Int64, Int64}}},Symbol})   # time: 0.009960893
    Base.precompile(Tuple{typeof(function_string),Type,AffExpr,Bool})   # time: 0.009605979
    Base.precompile(Tuple{typeof(_moi_fix),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},VariableRef,Float64,Bool})   # time: 0.009475288
    Base.precompile(Tuple{typeof(_constraint_macro),Tuple{Symbol, Expr},Symbol,typeof(parse_SD_constraint_expr),LineNumberNode})   # time: 0.009388473
    Base.precompile(Tuple{typeof(reshape_vector),Vector{VariableRef},SymmetricMatrixShape})   # time: 0.009321982
    Base.precompile(Tuple{typeof(parse_constraint_expr),Function,Expr})   # time: 0.009186855
    Base.precompile(Tuple{Core.kwftype(typeof(constraint_string)),NamedTuple{(:in_math_mode,), Tuple{Bool}},typeof(constraint_string),Type,ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable, MathOptInterface.LessThan{Float64}}, ScalarShape}})   # time: 0.009008045
    Base.precompile(Tuple{typeof(MathOptInterface.eval_constraint_jacobian),NLPEvaluator,SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true},Vector{Float64}})   # time: 0.00857948
    Base.precompile(Tuple{Core.kwftype(typeof(constraint_string)),NamedTuple{(:in_math_mode,), Tuple{Bool}},typeof(constraint_string),Type,ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}})   # time: 0.007398971
    Base.precompile(Tuple{typeof(_is_info_keyword),Expr})   # time: 0.005932119
    Base.precompile(Tuple{typeof(_constraint_macro),Tuple{Symbol, Expr, Expr},Symbol,typeof(parse_constraint_expr),LineNumberNode})   # time: 0.005548552
    isdefined(JuMP, Symbol("#90#99")) && Base.precompile(Tuple{getfield(JuMP, Symbol("#90#99")),Expr})   # time: 0.00545616
    Base.precompile(Tuple{typeof(all_constraints),Model,Type{QuadExpr},Type{MathOptInterface.LessThan{Float64}}})   # time: 0.005269228
    isdefined(JuMP, Symbol("#90#99")) && Base.precompile(Tuple{getfield(JuMP, Symbol("#90#99")),Symbol})   # time: 0.005126045
    Base.precompile(Tuple{typeof(setindex!),Model,Symmetric{VariableRef, Matrix{VariableRef}},Symbol})   # time: 0.004744124
    Base.precompile(Tuple{typeof(set_integer),VariableRef})   # time: 0.004611087
    isdefined(JuMP, Symbol("#90#99")) && Base.precompile(Tuple{getfield(JuMP, Symbol("#90#99")),Any})   # time: 0.004526849
    Base.precompile(Tuple{typeof(setindex!),Model,JuMP.Containers.DenseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape}, 3, Tuple{StepRange{Int64, Int64}, StepRange{Int64, Int64}, Base.OneTo{Int64}}, Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}, Dict{Int64, Int64}}},Symbol})   # time: 0.004247034
    Base.precompile(Tuple{typeof(all_constraints),Model,Type{VariableRef},Type{MathOptInterface.GreaterThan{Float64}}})   # time: 0.004184557
    Base.precompile(Tuple{typeof(all_constraints),Model,Type{VariableRef},Type{MathOptInterface.LessThan{Float64}}})   # time: 0.003693318
    Base.precompile(Tuple{typeof(MathOptInterface.eval_constraint),NLPEvaluator,SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true},Vector{Float64}})   # time: 0.003288632
    Base.precompile(Tuple{typeof(all_constraints),Model,Type{AffExpr},Type{MathOptInterface.EqualTo{Float64}}})   # time: 0.003248726
    Base.precompile(Tuple{typeof(parse_variable),Function,_VariableInfoExpr,Int64,Symbol,Symbol,Symbol,Float64})   # time: 0.003094496
    Base.precompile(Tuple{typeof(_moi_add_variable),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},Model,ScalarVariable{Int64, Int64, Float64, Float64},String})   # time: 0.003087982
    Base.precompile(Tuple{typeof(set_objective_coefficient),Model,VariableRef,Float64})   # time: 0.003028207
    Base.precompile(Tuple{typeof(setindex!),Model,Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.GreaterThan{Float64}}, ScalarShape}},Symbol})   # time: 0.002985051
    isdefined(JuMP, Symbol("#_error#101")) && Base.precompile(Tuple{getfield(JuMP, Symbol("#_error#101")),String})   # time: 0.002983753
    isdefined(JuMP, Symbol("#_error#101")) && Base.precompile(Tuple{getfield(JuMP, Symbol("#_error#101")),String})   # time: 0.002792419
    Base.precompile(Tuple{typeof(parse_variable),Function,_VariableInfoExpr,Expr,Symbol,Expr,Symbol,Expr})   # time: 0.002732017
    Base.precompile(Tuple{typeof(setindex!),Model,JuMP.Containers.SparseAxisArray{VariableRef, 3, Tuple{String, String, Int64}},Symbol})   # time: 0.002709768
    Base.precompile(Tuple{typeof(_moi_fix),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},VariableRef,Int64,Bool})   # time: 0.002553926
    Base.precompile(Tuple{typeof(setindex!),Model,Matrix{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape}},Symbol})   # time: 0.002535749
    Base.precompile(Tuple{typeof(setindex!),Model,Matrix{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}},Symbol})   # time: 0.002429434
    Base.precompile(Tuple{typeof(all_constraints),Model,Type{AffExpr},Type{MathOptInterface.LessThan{Float64}}})   # time: 0.002425839
    Base.precompile(Tuple{typeof(setindex!),Model,Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}},Symbol})   # time: 0.002409321
    Base.precompile(Tuple{typeof(_name_call),Expr,Vector{Any}})   # time: 0.002405974
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:lower_bound,), Tuple{Int64}},Type{_VariableInfoExpr}})   # time: 0.002381082
    Base.precompile(Tuple{typeof(MathOptInterface.jacobian_structure),NLPEvaluator})   # time: 0.002370995
    Base.precompile(Tuple{typeof(_name_call),String,Vector{Any}})   # time: 0.002361506
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:start,), Tuple{Float64}},Type{_VariableInfoExpr}})   # time: 0.002346172
    Base.precompile(Tuple{typeof(MathOptInterface.eval_objective),NLPEvaluator,Vector{Float64}})   # time: 0.002276394
    Base.precompile(Tuple{typeof(MutableArithmetics.mutable_operate!),typeof(MutableArithmetics.sub_mul),AffExpr,VariableRef,UInt16,Float64})   # time: 0.002249377
    Base.precompile(Tuple{typeof(build_constraint),Function,Matrix{AffExpr},PSDCone})   # time: 0.002245585
    Base.precompile(Tuple{typeof(parse_one_operator_constraint),Function,Bool,Val{:(==)},Expr,Expr})   # time: 0.00220451
    Base.precompile(Tuple{typeof(parse_one_operator_constraint),Function,Bool,Val{:>=},Expr,Expr})   # time: 0.001952178
    Base.precompile(Tuple{typeof(_moi_set_integer),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},VariableRef})   # time: 0.001938481
    Base.precompile(Tuple{typeof(parse_one_operator_constraint),Function,Bool,Val{:<=},Expr,Int64})   # time: 0.001845112
    Base.precompile(Tuple{typeof(set_name),ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.GreaterThan{Float64}}, ScalarShape},String})   # time: 0.001767766
    Base.precompile(Tuple{typeof(set_name),ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape},String})   # time: 0.001764125
    Base.precompile(Tuple{typeof(_add_kw_args),Expr,Vector{Any}})   # time: 0.001574977
    Base.precompile(Tuple{typeof(_moi_add_variable),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},Model,ScalarVariable{Float64, Float64, Float64, Float64},String})   # time: 0.001559715
    Base.precompile(Tuple{typeof(MathOptInterface.eval_objective),_UserFunctionEvaluator,SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}})   # time: 0.001388102
    Base.precompile(Tuple{typeof(_moi_get_result),MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}},MathOptInterface.ObjectiveValue})   # time: 0.001357983
    Base.precompile(Tuple{typeof(constraint_ref_with_index),Model,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable, MathOptInterface.LessThan{Float64}}})   # time: 0.001327407
    Base.precompile(Tuple{typeof(MutableArithmetics.sub_mul),VariableRef,Float64,AffExpr})   # time: 0.001269834
    Base.precompile(Tuple{typeof(set_name),ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape},String})   # time: 0.001234877
    Base.precompile(Tuple{typeof(convert),Type{AffExpr},Union{Number, UniformScaling}})   # time: 0.00121725
    Base.precompile(Tuple{Core.kwftype(typeof(_macro_assign_and_return)),NamedTuple{(:model_for_registering,), Tuple{Expr}},typeof(_macro_assign_and_return),Expr,Symbol,Symbol})   # time: 0.001182317
    isdefined(JuMP, Symbol("#110#111")) && Base.precompile(Tuple{getfield(JuMP, Symbol("#110#111")),Tuple{Int64, MathOptInterface.VariableIndex}})   # time: 0.00111902
    Base.precompile(Tuple{typeof(set_normalized_coefficient),ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.GreaterThan{Float64}}, ScalarShape},VariableRef,Int64})   # time: 0.00108381
    Base.precompile(Tuple{typeof(parse_constraint_head),Function,Val{:call},Symbol,Symbol,Int64})   # time: 0.001021371
end
