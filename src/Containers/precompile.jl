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
    Base.precompile(Tuple{typeof(show),Base.TTY,DenseAxisArray{Float64, 2, Tuple{UnitRange{Int64}, Vector{Symbol}}, Tuple{Dict{Int64, Int64}, Dict{Symbol, Int64}}}})   # time: 0.06925241
    Base.precompile(Tuple{typeof(_extract_kw_args),Tuple{Symbol, Expr}})   # time: 0.02388525
    Base.precompile(Tuple{typeof(container),Function,VectorizedProductIterator{Tuple{Int64}}})   # time: 0.019638063
    Base.precompile(Tuple{typeof(_extract_kw_args),Tuple{Expr, Symbol}})   # time: 0.01300304
    Base.precompile(Tuple{typeof(container),Function,VectorizedProductIterator{Tuple{UnitRange{Int64}, Vector{Symbol}}},Type{DenseAxisArray}})   # time: 0.011042849
    Base.precompile(Tuple{typeof(_extract_kw_args),Tuple{}})   # time: 0.010988676
    Base.precompile(Tuple{typeof(container),Function,VectorizedProductIterator{Tuple{StepRange{Int64, Int64}, StepRange{Int64, Int64}, Base.OneTo{Int64}}},Type{DenseAxisArray}})   # time: 0.010785832
    Base.precompile(Tuple{typeof(container),Function,VectorizedProductIterator{Tuple{Base.OneTo{Int64}, UnitRange{Int64}, UnitRange{Int64}}},Type{DenseAxisArray}})   # time: 0.010772673
    Base.precompile(Tuple{typeof(_extract_kw_args),Tuple{Expr}})   # time: 0.010699514
    Base.precompile(Tuple{typeof(container),Function,VectorizedProductIterator{Tuple{Vector{String}}},Type{DenseAxisArray}})   # time: 0.009808394
    Base.precompile(Tuple{typeof(_extract_kw_args),Tuple{Symbol}})   # time: 0.009082955
    Base.precompile(Tuple{typeof(container),Function,VectorizedProductIterator{Tuple{UnitRange{Int64}}},Type{DenseAxisArray}})   # time: 0.008933888
    Base.precompile(Tuple{typeof(_extract_kw_args),Tuple{Symbol, Expr, Expr}})   # time: 0.007755594
    Base.precompile(Tuple{typeof(_extract_kw_args),Tuple{Expr, Expr}})   # time: 0.007691665
    Base.precompile(Tuple{typeof(_extract_kw_args),Tuple{Expr, Expr, Expr}})   # time: 0.006882744
    Base.precompile(Tuple{typeof(container_code),Vector{Any},Expr,Expr,Symbol})   # time: 0.005621142
    # TODO: isdefined(JuMP.Containers, Symbol("#37#38")) && Base.precompile(Tuple{getfield(JuMP.Containers, Symbol("#37#38")),Int64})   # time: 0.004817874
    Base.precompile(Tuple{typeof(to_index),DenseAxisArray,Any})   # time: 0.004183273
    Base.precompile(Tuple{typeof(vectorized_product),Base.OneTo{Int64},Vararg{Any, N} where N})   # time: 0.00327413
    Base.precompile(Tuple{typeof(_to_index_tuple),Tuple{Int64, Vararg{Any, N} where N},Tuple{Dict{Int64, Int64}, Dict{Symbol, Int64}}})   # time: 0.003025034
    let fbody = try __lookup_kwbody__(which(nested, (Function,Vararg{Function, N} where N,))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Function,typeof(nested),Function,Vararg{Function, N} where N,))
        end
    end   # time: 0.002819959
    Base.precompile(Tuple{typeof(haskey),NoDuplicateDict{_A, _B} where {_A, _B},Tuple{Int64, Int64, Int64}})   # time: 0.002722087
    Base.precompile(Tuple{typeof(vectorized_product),Base.OneTo{Int64},Vararg{Base.OneTo{Int64}, N} where N})   # time: 0.002584661
    Base.precompile(Tuple{typeof(_to_index_tuple),Tuple{Int64, Function},Tuple{Dict{Int64, Int64}, Dict{Symbol, Int64}}})   # time: 0.002516737
    Base.precompile(Tuple{typeof(_to_index_tuple),Tuple{String, Int64, Int64},Tuple{Dict{String, Int64}, Dict{Int64, Int64}, Dict{Int64, Int64}}})   # time: 0.002253078
    Base.precompile(Tuple{typeof(has_dependent_sets),Vector{Any},Vector{Any}})   # time: 0.00217896
    # TODO: Base.precompile(Tuple{Type{DenseAxisArray},Core.Array{T, N},Any,Tuple{Vararg{Dict, N}} where N})   # time: 0.001977609
    Base.precompile(Tuple{typeof(_explicit_oneto),Expr})   # time: 0.001850751
    Base.precompile(Tuple{typeof(haskey),NoDuplicateDict{_A, _B} where {_A, _B},Tuple{String, String, Int64}})   # time: 0.001788105
    Base.precompile(Tuple{typeof(vectorized_product),Int64})   # time: 0.001391668
    Base.precompile(Tuple{typeof(_explicit_oneto),Any})   # time: 0.001387689
    Base.precompile(Tuple{typeof(build_lookup),Vector{Int64}})   # time: 0.001371672
    Base.precompile(Tuple{typeof(depends_on),Expr,Symbol})   # time: 0.001150738
    Base.precompile(Tuple{typeof(has_colon),Tuple{Int64, Function}})   # time: 0.001108028
    Base.precompile(Tuple{typeof(container),Function,VectorizedProductIterator{Tuple{Vector{String}, Base.OneTo{Int64}}},Type{DenseAxisArray}})   # time: 0.001101885
end
