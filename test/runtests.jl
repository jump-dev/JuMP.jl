using ReverseDiffSparse2
using Base.Test
using FactCheck


ex = :(sin(x[1]^2) + cos(x[2]*4)/5-2.0)

nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)

#@show nd

storage = zeros(length(nd))
reverse_storage = zeros(length(nd))

x = [2.0,3.0]
#@show x
fval = forward_eval(storage,nd,adj,const_values,x,[])
true_val = sin(x[1]^2) + cos(x[2]*4)/5 -2.0
@test isapprox(fval,true_val)

grad = zeros(2)
reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,[])

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
fval = forward_eval(outer_storage,nd_outer,adj_outer,[],x,[fval])
@test isapprox(fval,true_val)

outer_reverse_storage = zeros(1)
fill!(grad,0.0)
subexpr_output = zeros(1)
reverse_eval(grad,outer_reverse_storage,outer_storage,nd_outer,adj_outer,[],subexpr_output)
@assert subexpr_output[1] == 1.0
reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,[],subexpr_output[1])
@test isapprox(grad,true_grad)


ex = :((1/x[1])^x[2]-x[3])

nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)

storage = zeros(length(nd))
reverse_storage = zeros(length(nd))

x = [2.5,3.5,1.0]
#@show x
fval = forward_eval(storage,nd,adj,const_values,x,[])
true_val = (1/x[1])^x[2]-x[3]
@test isapprox(fval,true_val)

grad = zeros(3)
reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,[])

true_grad = [-x[2]*x[1]^(-x[2]-1), -((1/x[1])^x[2])*log(x[1]),-1]
@test isapprox(grad,true_grad)


import ReverseDiffSparse2: CONSTANT, LINEAR, NONLINEAR

function test_linearity(ex,testval,IJ = [],indices=[])
    nd,const_values = expr_to_nodedata(ex)
    adj = adjmat(nd)
    linearity = classify_linearity(nd,adj)
    @test linearity[1] == testval
    edgelist = compute_hessian_sparsity(nd,adj,linearity)
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

using DualNumbers

ex = :(x[1]^2/2 + 2x[1]*x[2])
nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)
forward_storage = Array(Dual{Float64},length(nd))
reverse_storage = Array(Dual{Float64},length(nd))
reverse_output_vector = Array(Dual{Float64},2)
forward_input_vector = Array(Dual{Float64},2)
x_values = [10.0,2.0]

local_to_global_idx = [1,2]
R = [1.0 0.0; 0.0 1.0]
hessmat_eval!(R, reverse_storage, forward_storage, nd, adj, const_values, x_values, reverse_output_vector, forward_input_vector, local_to_global_idx)
@test R == [1.0 2.0; 2.0 0.0]
# now with a permutation
local_to_global_idx = [2,1]
R = [0.0 1.0; 1.0 0.0]
hessmat_eval!(R, reverse_storage, forward_storage, nd, adj, const_values, x_values, reverse_output_vector, forward_input_vector, local_to_global_idx)
@test R == [2.0 0.0; 1.0 2.0]


include("test_coloring.jl")
#include("test_jump.jl")
# FactCheck.exitstatus() # ignore these errors for now
