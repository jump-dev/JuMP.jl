using ReverseDiffSparse
using Base.Test
using FactCheck


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
reverse_eval(grad,reverse_storage,partials_storage,nd,adj,[],1.0)

true_grad = [2*x[1]*cos(x[1]^2), -4*sin(x[2]*4)/5]
@test isapprox(grad,true_grad)

# subexpressions

nd_outer = [NodeData(SUBEXPRESSION,1,-1,-1)]
li = list_subexpressions(nd_outer)
@test li == [1]
li,li_individual = order_subexpressions(Vector{NodeData}[nd_outer], Vector{NodeData}[nd])
@test li == [1]
@test li_individual[1] == [1]

nd_outer2 = [NodeData(CALL,operator_to_id[:+],-1,-1),NodeData(SUBEXPRESSION,2,1,1),NodeData(SUBEXPRESSION,1,1,2)]
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
reverse_eval(grad,outer_reverse_storage,outer_storage_partials,nd_outer,adj_outer,subexpr_output,1.0)
@assert subexpr_output[1] == 1.0
reverse_eval(grad,reverse_storage,partials_storage,nd,adj,[],subexpr_output[1])
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
reverse_eval(grad,reverse_storage,partials_storage,nd,adj,[],1.0)

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
reverse_eval(grad,reverse_storage,partials_storage,nd,adj,[],1.0)
@test grad[1] == 1

fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[-0.1],[])
@test fval == -0.1
fill!(grad,0)
reverse_eval(grad,reverse_storage,partials_storage,nd,adj,[],1.0)
@test grad[1] == 1

fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[0.2],[])
@test fval == 5
fill!(grad,0)
reverse_eval(grad,reverse_storage,partials_storage,nd,adj,[],1.0)
@test grad[1] == 0

# parameters
nd = [NodeData(PARAMETER,1,-1,-1)]
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
reverse_eval(grad,reverse_storage,partials_storage,nd,adj,[],1.0)
@test grad[1] == 1.0


fval = forward_eval(storage,partials_storage,nd,adj,const_values,[],[-2.0],[])
@test fval == 2.0
grad = zeros(1)
reverse_eval(grad,reverse_storage,partials_storage,nd,adj,[],1.0)
@test grad[1] == -1.0


function test_linearity(ex,testval,IJ = [],indices=[])
    nd,const_values = expr_to_nodedata(ex)
    adj = adjmat(nd)
    linearity = classify_linearity(nd,adj,[])
    @test linearity[1] == testval
    edgelist = compute_hessian_sparsity(nd,adj,linearity,Array(Set{Tuple{Int,Int}},0), Array(Set{Int},0))
    if linearity[1] != NONLINEAR
        @test length(edgelist) == 0
    elseif length(IJ) > 0
        @test IJ == edgelist
    end
    if length(indices) > 0
        ix = sort(collect(compute_gradient_sparsity(nd)))
        @test ix == indices
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


using DualNumbers

ex = :(x[1]^2/2 + 2x[1]*x[2])
nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)
forward_storage = Array(Dual{Float64},length(nd))
reverse_storage = Array(Dual{Float64},length(nd))
partials_storage = Array(Dual{Float64},length(nd))
reverse_output_vector = Array(Dual{Float64},2)
forward_input_vector = Array(Dual{Float64},2)
x_values = [10.0,2.0]

local_to_global_idx = [1,2]
R = [1.0 0.0; 0.0 1.0]
hessmat_eval!(R, reverse_storage, forward_storage, partials_storage, nd, adj, const_values, x_values, reverse_output_vector, forward_input_vector, local_to_global_idx)
@test R == [1.0 2.0; 2.0 0.0]
# now with a permutation
local_to_global_idx = [2,1]
R = [0.0 1.0; 1.0 0.0]
hessmat_eval!(R, reverse_storage, forward_storage, partials_storage, nd, adj, const_values, x_values, reverse_output_vector, forward_input_vector, local_to_global_idx)
@test R == [2.0 0.0; 1.0 2.0]

using ForwardDiff

# dual forward test
function dualforward(ex, x)
    nd,const_values = expr_to_nodedata(ex)
    adj = adjmat(nd)
    forward_storage = zeros(length(nd))
    partials_storage = zeros(length(nd))
    reverse_storage = zeros(length(nd))

    fval = forward_eval(forward_storage,partials_storage,nd,adj,const_values,[],x,[])
    grad = zeros(length(x))
    reverse_eval(grad,reverse_storage,partials_storage,nd,adj,[],1.0)

    zero_ϵ = ForwardDiff.zero_partials(NTuple{1,Float64},1)
    forward_storage_ϵ = fill(zero_ϵ,length(nd))
    partials_storage_ϵ = fill(zero_ϵ,length(nd))
    x_values_ϵ = fill(ForwardDiff.Partials((1.0,)),length(x))
    fval_ϵ = forward_eval_ϵ(forward_storage,forward_storage_ϵ,partials_storage,partials_storage_ϵ,nd,adj,x_values_ϵ,[])
    @test_approx_eq fval_ϵ.data[1] dot(grad,ones(length(x)))

    # compare with running dual numbers
    forward_dual_storage = zeros(Dual{Float64},length(nd))
    partials_dual_storage = zeros(Dual{Float64},length(nd))
    x_dual = [Dual(x[i],1.0) for i in 1:length(x)]
    fval = forward_eval(forward_dual_storage,partials_dual_storage,nd,adj,const_values,[],x_dual,[])
    for k in 1:length(nd)
        @test_approx_eq epsilon(forward_dual_storage[k]) forward_storage_ϵ[k].data[1]
    end
    for k in 1:length(nd)
        @test_approx_eq epsilon(partials_dual_storage[k]) partials_storage_ϵ[k].data[1]
    end
end

dualforward(:(sin(x[1]^2) + cos(x[2]*4)/5-2.0),[1.0,2.0])
dualforward(:(sin(x[1]^x[2]) + cos(x[2]*4)/5-2.0),[1.0,2.0])
dualforward(:(sin(x[1]^3) + cos(x[1]*x[2]*4)/5-2.0),[1.0,0.0])


include("test_coloring.jl")
#include("test_jump.jl")
# FactCheck.exitstatus() # ignore these errors for now
