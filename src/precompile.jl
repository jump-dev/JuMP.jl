#=
BEFORE
start 0.000290968
using 0.39690596
model 0.113112627
defva 0.325597146
addcn 0.989031469
=#

# JuMP.jl
Base.precompile(Variable, (Model,Float64,Float64,Symbol,String,Float64))
Base.precompile(Variable, (Model,Float64,Int,Symbol,String,Float64))
Base.precompile(Variable, (Model,Int,Float64,Symbol,String,Float64))
Base.precompile(Variable, (Model,Int,Int,Symbol,String,Float64))
Base.precompile(addConstraint, (Model,LinearConstraint))
# macros.jl
Base.precompile(buildrefsets, (Expr,))
# parseExpr_staged.jl
Base.precompile(addToExpression_reorder, (AffExpr,(Int,Variable)))
Base.precompile(addToExpression_reorder, (AffExpr,(Float64,Int)))
Base.precompile(parseCurly, (Expr,Symbol,Vector{Any},Vector{Any},Symbol))
Base.precompile(parseSum,   (Expr,Symbol,Vector{Any},Vector{Any},Symbol))
Base.precompile(parseExpr,  (Expr,Symbol,Vector{Any},Vector{Any},Symbol))

#=
AFTER
start 0.00034573
using 0.410585875
model 0.115268569
defva 0.33364511
addcn 0.559814242
=#