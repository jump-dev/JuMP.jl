function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    # TODO: Base.precompile(Tuple{typeof(eval_and_check_return_type),Function,Type,JuMP._UserFunctionEvaluator,Vararg{Any, N} where N})   # time: 0.003078309
end
