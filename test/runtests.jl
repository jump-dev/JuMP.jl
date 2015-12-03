using ReverseDiffSparse2
using Base.Test


ex = :(sin(x[1]^2) + cos(x[2]*4)-2.0)

nd,const_values = expr_to_nodedata(ex)
adj = adjmat(nd)

#@show nd

storage = zeros(length(nd))
reverse_storage = zeros(length(nd))

for k in 1:3
    x = rand(2)
    #@show x
    fval = forward_eval(storage,nd,adj,const_values,x)
    true_val = sin(x[1]^2) + cos(x[2]*4) -2.0
    @test isapprox(fval,true_val)

    grad = zeros(2)
    reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,x)

    true_grad = [2*x[1]*cos(x[1]^2), -4*sin(x[2]*4)]
    @test isapprox(grad,true_grad)

end
