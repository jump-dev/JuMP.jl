using ReverseDiffSparse2
using Base.Test


ex = :(sin(x[1]^2) + cos(x[2]*4)/5-2.0)

nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)

#@show nd

storage = zeros(length(nd))
reverse_storage = zeros(length(nd))

for k in 1:3
    x = rand(2)
    #@show x
    fval = forward_eval(storage,nd,adj,const_values,x)
    true_val = sin(x[1]^2) + cos(x[2]*4)/5 -2.0
    @test isapprox(fval,true_val)

    grad = zeros(2)
    reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,x)

    true_grad = [2*x[1]*cos(x[1]^2), -4*sin(x[2]*4)/5]
    @test isapprox(grad,true_grad)

end

ex = :((1/x[1])^x[2]-x[3])

nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)

storage = zeros(length(nd))
reverse_storage = zeros(length(nd))

x = [2.5,3.5,1.0]
#@show x
fval = forward_eval(storage,nd,adj,const_values,x)
true_val = (1/x[1])^x[2]-x[3]
@test isapprox(fval,true_val)

grad = zeros(3)
reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,x)

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
        ix = compute_gradient_sparsity(nd, adj)
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
