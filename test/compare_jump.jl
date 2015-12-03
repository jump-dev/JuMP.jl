using JuMP
using Base.Test
using ReverseDiffSparse2

function compare_jump(m::Model)

    nl = JuMPNLPEvaluator(m)
    nvar = MathProgBase.numvar(m)
    x = rand(nvar)
    println("JuMP init time:")
    @time MathProgBase.initialize(nl, [:Grad,:ExprGraph])



    objexpr = MathProgBase.obj_expr(nl)

    println("RDS2 init time:")

    @time nd,const_values = expr_to_nodedata(objexpr)
    @time adj = adjmat(nd)

    storage = zeros(length(nd))
println("RDS2 Feval time")
    forward_eval(storage,nd,adj,const_values,x)
    @time fval = forward_eval(storage,nd,adj,const_values,x)

    @show fval
    MathProgBase.eval_f(nl,x)
println("JuMP Feval time")
    @time fval = MathProgBase.eval_f(nl,x)
    @show fval
    
    reverse_storage = zeros(length(nd))
    grad = zeros(nvar)
    ∇f = zeros(nvar)
    MathProgBase.eval_grad_f(nl, ∇f, x)
    reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,x)
    @test isapprox(grad,∇f)
    println("RDS2 grad time")
    @time reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,x)
    @time reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,x)
    println("JuMP grad time")
    @time MathProgBase.eval_grad_f(nl, ∇f, x)
    @time MathProgBase.eval_grad_f(nl, ∇f, x)
    for i in 1:10
        @profile reverse_eval(grad,reverse_storage,storage,nd,adj,const_values,x)
    end
end

function mod1(N)

    m = Model()
    @defVar(m, x[1:N])
    c = 10

    @setNLObjective(m, Min, c*sum{x[1]*sum{sin(x[i]);iseven(i)},i=1:N})

    return m
end

#m = mod1(100000)
#compare_jump(m)

