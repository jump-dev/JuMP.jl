using JuMP
import Ipopt
import Random

function nlp_model()
    k = 800
    n = 12
    Random.seed!(1)
    p = rand(400:700, k, 1)
    c1 = rand(100:200, k, n)
    c2 = 0.9 .* c1
    b = rand(150:250, k, 1)
    model = Model(Ipopt.Optimizer)
    @variable(model, 0 <= x[i = 1:n] <= 1)
    @variable(model, 0 <= var1 <= 1)
    @variable(model, 0 <= var2 <= 1)
    @variable(model, 0 <= var3 <= 1)
    @objective(model, Max, var1 - var2 + var3)
    @NLexpression(model, expr, sum(x[i] * p[i] for i in 1:n))
    @NLexpression(model, expr_c1[j = 1:k], sum(x[i] * c1[j, i] for i in 1:n))
    @NLexpression(model, expr_c2[j = 1:k], sum(x[i] * c2[j, i] for i in 1:n))
    @NLconstraint(model, expr - sum(b[j] / (1 + var1)^j for j in 1:k) == 0)
    @NLconstraint(
        model,
        expr - sum(expr_c1[j] / (1 + var2)^j for j in 1:k) == 0,
    )
    @NLconstraint(
        model,
        expr - sum(expr_c2[j] / (1 + var3)^j for j in 1:k) == 0,
    )
    @NLconstraint(model, [j = 1:k], expr_c1[j] >= b[j])
    @time optimize!(model)
    return model
end

@time nlp_model();

using ProfileView
@profview nlp_model();
