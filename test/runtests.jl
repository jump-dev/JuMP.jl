using ReverseDiffSparse
using Base.Test


ex = :(sin(x[1]^2) + cos(x[2]*4)/5-2.0)

nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)

#@show nd

storage = zeros(length(nd))
partials_storage = zeros(length(nd))
reverse_storage = zeros(length(nd))

x = [2.0,3.0]
#@show x
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],x,[])
true_val = sin(x[1]^2) + cos(x[2]*4)/5 -2.0
@test isapprox(fval,true_val)

grad = zeros(2)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)

true_grad = [2*x[1]*cos(x[1]^2), -4*sin(x[2]*4)/5]
@test isapprox(grad,true_grad)

# subexpressions

nd_outer = [NodeData(SUBEXPRESSION,1,-1)]
li = list_subexpressions(nd_outer)
@test li == [1]
li,li_individual = order_subexpressions(Vector{NodeData}[nd_outer], Vector{NodeData}[nd])
@test li == [1]
@test li_individual[1] == [1]

nd_outer2 = [NodeData(CALL,operator_to_id[:+],-1),NodeData(SUBEXPRESSION,2,1),NodeData(SUBEXPRESSION,1,1)]
li = list_subexpressions(nd_outer2)
@test li == [1,2]
li,li_individual = order_subexpressions(Vector{NodeData}[nd_outer2], Vector{NodeData}[nd,nd_outer])
@test li == [1,2]
@test li_individual[1] == [1,2]

adj_outer = adjmat(nd_outer)
outer_storage = zeros(1)
outer_storage_partials = zeros(1)
fval = forward_eval(outer_storage,outer_storage_partials,nd_outer,adj_outer,[],[],x,[fval])
@test isapprox(fval,true_val)

outer_reverse_storage = zeros(1)
fill!(grad,0.0)
subexpr_output = zeros(1)
reverse_eval(outer_reverse_storage,outer_storage_partials,nd_outer,adj_outer)
reverse_extract(grad,outer_reverse_storage,nd_outer,adj_outer,subexpr_output,1.0)
@assert subexpr_output[1] == 1.0
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],subexpr_output[1])
@test isapprox(grad,true_grad)


ex = :((1/x[1])^x[2]-x[3])

nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)

storage = zeros(length(nd))
partials_storage = zeros(length(nd))
reverse_storage = zeros(length(nd))

x = [2.5,3.5,1.0]
#@show x
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],x,[])
true_val = (1/x[1])^x[2]-x[3]
@test isapprox(fval,true_val)

grad = zeros(3)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)

true_grad = [-x[2]*x[1]^(-x[2]-1), -((1/x[1])^x[2])*log(x[1]),-1]
@test isapprox(grad,true_grad)

# logical expressions
ex = :(x[1] > 0.5 && x[1] < 0.9)
nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)
storage = zeros(length(nd))
partials_storage = zeros(length(nd))
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[1.5],[])
@test fval == 0
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[0.6],[])
@test fval == 1

ex = :(ifelse(x[1] >= 0.5 || x[1] <= 0.1,x[1],5))
nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)
storage = zeros(length(nd))
partials_storage = zeros(length(nd))
reverse_storage = zeros(length(nd))
grad = zeros(1)
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[1.5],[])
@test fval == 1.5
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)
@test grad[1] == 1

fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[-0.1],[])
@test fval == -0.1
fill!(grad,0)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)
@test grad[1] == 1

fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[0.2],[])
@test fval == 5
fill!(grad,0)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)
@test grad[1] == 0

ex = :(sqrt(x[1]))
nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)
storage = zeros(length(nd))
partials_storage = zeros(length(nd))
reverse_storage = zeros(length(nd))
grad = zeros(1)
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[-1.5],[])
@test isnan(fval)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)
@test isnan(grad[1])


# parameters
nd = [NodeData(PARAMETER,1,-1)]
adj = adjmat(nd)
storage = zeros(length(nd))
partials_storage = zeros(length(nd))
fval = forward_eval(storage,partials_storage,nd,adj,[],[105.2],[-0.1],[])
@test fval == 105.2

# abs
ex = :(abs(x[1]))

nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)

storage = zeros(length(nd))
partials_storage = zeros(length(nd))
reverse_storage = zeros(length(nd))
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[2.0],[])
@test fval == 2.0
grad = zeros(1)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)
@test grad[1] == 1.0


fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[-2.0],[])
@test fval == 2.0
grad = zeros(1)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)
@test grad[1] == -1.0

# https://github.com/JuliaOpt/JuMP.jl/issues/855
ex = :(ifelse(x[1]<=3.0, (x[1]-2.0)^2, 2*log(x[1]-2.0)+1.0))

nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)

storage = zeros(length(nd))
partials_storage = zeros(length(nd))
reverse_storage = zeros(length(nd))
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[-1.0],[])
@test fval == 9.0
grad = zeros(1)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)
@test grad[1] == -6.0

fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[2.0],[])
@test fval == 0.0
grad = zeros(1)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)
@test grad[1] == 0.0


function test_linearity(ex,testval,IJ = [],indices=[])
    nd,const_values = expr_to_nodedata(ex)
    adj = adjmat(nd)
    linearity = classify_linearity(nd,adj,[])
    @test linearity[1] == testval
    idxset = Coloring.IndexedSet(100)
    edgelist = compute_hessian_sparsity(nd,adj,linearity,idxset,Array{Set{Tuple{Int,Int}}}(0), Array{Vector{Int}}(0))
    if linearity[1] != NONLINEAR
        @test length(edgelist) == 0
    elseif length(IJ) > 0
        @test IJ == edgelist
    end
    if length(indices) > 0
        ix = sort(collect(compute_gradient_sparsity(nd)))
        compute_gradient_sparsity!(idxset, nd)
        ix2 = sort(collect(idxset))
        empty!(idxset)
        @test ix == indices
        @test ix == ix2
    end
end

test_linearity(:(sin(x[1]^2) + cos(x[2]*4)-2.0), NONLINEAR, Set([(2,2),(1,1)]), [1,2])
test_linearity(:(3*4*(x[1]+x[2])), LINEAR)
test_linearity(:(x[3]*x[2]), NONLINEAR, Set([(3,2),(3,3),(2,2)]),[2,3])
test_linearity(:(3+4), CONSTANT)
test_linearity(:(sin(3)+x[1]), LINEAR)
test_linearity(:(cos(x[3])*sin(3)+x[1]), NONLINEAR, Set([(3,3)]),[1,3])
test_linearity(:(x[1]-3x[2]), LINEAR)
test_linearity(:(-x[1]), LINEAR)
test_linearity(:(+x[1]), LINEAR)
test_linearity(:(x[1]^x[2]), NONLINEAR, Set([(2,2),(1,1),(2,1)]))
test_linearity(:(x[1]/3+x[2]), LINEAR)
test_linearity(:(3/(x[1]*x[2])), NONLINEAR, Set([(2,2),(1,1),(2,1)]))
test_linearity(:(1/(x[1]+3)), NONLINEAR)
test_linearity(:(ifelse(x[1] <= 1,x[1],x[2])), PIECEWISE_LINEAR, Set([]))
test_linearity(:(ifelse(x[1] <= 1,x[1]^2,x[2])), NONLINEAR, Set([(1,1)]))
test_linearity(:(ifelse(1 <= 1,2,3)), CONSTANT)
test_linearity(:(1/ifelse(x[1] < 1, x[1],0)), NONLINEAR, Set([(1,1)]))

# eliminating fixed variables and constants
ex = :(sin(x[1]^2) + cos(x[2]*(2*2))/5-2.0)
nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)
storage = zeros(length(nd))
partials_storage = zeros(length(nd))
x = [2.0,3.0]
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],x,[])
linearity = classify_linearity(nd,adj,[],[false,true])
new_nd = simplify_constants(storage,nd,adj,const_values,linearity)
@test length(new_nd) < length(nd)
new_adj = adjmat(new_nd)
true_val = sin(x[1]^2) + cos(x[2]*4)/5 -2.0
x[2] = 100 # this shouldn't affect the answer because we said x[2] is fixed
fval = forward_eval(storage,partials_storage,new_nd,new_adj,const_values,[],x,[])
@test isapprox(fval,true_val)
# all variables fixed
linearity = classify_linearity(nd,adj,[],[true,true])
new_nd = simplify_constants(storage,nd,adj,const_values,linearity)
@test length(new_nd) == 1

# user-defined functions
using Distributions
#Φ(x,y) = Φ(y) - Φ(x)
# c(x) = cos(x)
type CDFEvaluator <: MathProgBase.AbstractNLPEvaluator
end
function MathProgBase.eval_f(::CDFEvaluator,x)
    @assert length(x) == 2
    return cdf(Normal(0,1),x[2])-cdf(Normal(0,1),x[1])
end
function MathProgBase.eval_grad_f(::CDFEvaluator,grad,x)
    grad[1] = -pdf(Normal(0,1),x[1])
    grad[2] = pdf(Normal(0,1),x[2])
end
r = ReverseDiffSparse.UserOperatorRegistry()
register_multivariate_operator!(r,:Φ,CDFEvaluator())
register_univariate_operator!(r,:c,cos,x->-sin(x),x->-cos(x))
Φ(x,y) = MathProgBase.eval_f(CDFEvaluator(),[x,y])
ex = :(Φ(x[2],x[1]-1)*c(x[3]))
nd,const_values = expr_to_nodedata(ex,r)
@test ReverseDiffSparse.has_user_multivariate_operators(nd)
adj = adjmat(nd)
storage = zeros(length(nd))
partials_storage = zeros(length(nd))
reverse_storage = zeros(length(nd))
x = [2.0,3.0,4.0]
fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],x,[],zeros(2),zeros(2),user_operators=r)
true_val = Φ(x[2],x[1]-1)*cos(x[3])
@test isapprox(fval,true_val)
grad = zeros(3)
reverse_eval(reverse_storage,partials_storage,nd,adj)
reverse_extract(grad,reverse_storage,nd,adj,[],1.0)

true_grad = [cos(x[3])*pdf(Normal(0,1),x[1]-1), -cos(x[3])*pdf(Normal(0,1),x[2]), -sin(x[3])*Φ(x[2],x[1]-1)]
@test isapprox(grad,true_grad)



using DualNumbers
using ForwardDiff

# dual forward test
function dualforward(ex, x; ignore_nan=false)
    nd,const_values = expr_to_nodedata(ex)
    adj = adjmat(nd)
    forward_storage = zeros(length(nd))
    partials_storage = zeros(length(nd))
    reverse_storage = zeros(length(nd))

    fval = forward_eval(forward_storage,partials_storage,nd,adj,const_values,[],x,[])
    grad = zeros(length(x))
    reverse_eval(reverse_storage,partials_storage,nd,adj)
    reverse_extract(grad,reverse_storage,nd,adj,[],1.0)

    zero_ϵ = zero(ForwardDiff.Partials{1,Float64})
    forward_storage_ϵ = fill(zero_ϵ,length(nd))
    partials_storage_ϵ = fill(zero_ϵ,length(nd))
    x_values_ϵ = fill(ForwardDiff.Partials((1.0,)),length(x))
    reverse_storage_ϵ = fill(zero_ϵ,length(nd))
    output_ϵ = fill(zero_ϵ,length(x))
    fval_ϵ = forward_eval_ϵ(forward_storage,forward_storage_ϵ,partials_storage,partials_storage_ϵ,nd,adj,x_values_ϵ,[])
    reverse_eval_ϵ(output_ϵ,reverse_storage,reverse_storage_ϵ,partials_storage,partials_storage_ϵ,nd,adj,[],[],2.0,zero_ϵ)
    @test isapprox(fval_ϵ[1], dot(grad,ones(length(x))))

    # compare with running dual numbers
    forward_dual_storage = zeros(DualNumbers.Dual{Float64},length(nd))
    partials_dual_storage = zeros(DualNumbers.Dual{Float64},length(nd))
    output_dual_storage = zeros(DualNumbers.Dual{Float64},length(x))
    reverse_dual_storage = zeros(DualNumbers.Dual{Float64},length(nd))
    x_dual = [DualNumbers.Dual(x[i],1.0) for i in 1:length(x)]
    fval = forward_eval(forward_dual_storage,partials_dual_storage,nd,adj,const_values,[],x_dual,[])
    reverse_eval(reverse_dual_storage,partials_dual_storage,nd,adj)
    reverse_extract(output_dual_storage,reverse_dual_storage,nd,adj,[],DualNumbers.Dual(2.0))
    for k in 1:length(nd)
        @test isapprox(epsilon(forward_dual_storage[k]), forward_storage_ϵ[k][1])
        if !(isnan(epsilon(partials_dual_storage[k])) && ignore_nan)
            @test isapprox(epsilon(partials_dual_storage[k]), partials_storage_ϵ[k][1])
        else
            @test !isnan(forward_storage_ϵ[k][1])
        end
        if !(isnan(epsilon(reverse_dual_storage[k])) && ignore_nan)
            @test isapprox(epsilon(reverse_dual_storage[k]), reverse_storage_ϵ[k][1]/2)
        else
            @test !isnan(reverse_storage_ϵ[k][1])
        end
    end
    for k in 1:length(x)
        if !(isnan(epsilon(output_dual_storage[k])) && ignore_nan)
            @test isapprox(epsilon(output_dual_storage[k]), output_ϵ[k][1])
        else
            @test !isnan(output_ϵ[k][1])
        end
    end
end

dualforward(:(sin(x[1]^2) + cos(x[2]*4)/5-2.0),[1.0,2.0])
dualforward(:(sin(x[1]^x[2]) + cos(x[2]*4)/5-2.0),[1.0,2.0])
dualforward(:(sin(x[1]^3) + cos(x[1]*x[2]*4)/5-2.0),[1.0,0.0])
dualforward(:(x[1]*x[2]),[3.427139283036299e-206,1.0], ignore_nan=true)


include("test_coloring.jl")
