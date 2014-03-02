using Base.Test
using ReverseDiffSparse

function test_sparsity(sp, H)
    I,J = sp
    Hsp = sparse(I,J,ones(length(I)))
    H = sparse(tril(H))
    @test all( Hsp.colptr .== H.colptr)
    @test all( Hsp.rowval .== H.rowval)
end

x,y,z,q = placeholders(4)

ex = @processNLExpr sin(x*y) + exp(z+2q)
sp = compute_hessian_sparsity_IJ(ex)
hfunc = gen_hessian_dense(ex)
val = [3.0, 4.0, 5.0, 6.0]
exact(x,y,z,q) = [ -y^2*sin(x*y) cos(x*y)-x*y*sin(x*y) 0 0
                   cos(x*y)-x*y*sin(x*y) -x^2*sin(x*y) 0 0
                   0 0 exp(z+2q) 2*exp(z+2q)
                   0 0 2*exp(z+2q) 4*exp(z+2q) ]
@test_approx_eq hfunc(val) exact(val...)
test_sparsity(sp, exact(val...))

sparsemat, sparsefunc = gen_hessian_sparse_mat(ex)
sparsefunc(val, sparsemat)
@test_approx_eq sparsemat tril(exact(val...))

I,J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex)
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq to_H(ex, I, J, V, 4) tril(exact(val...))

x = placeholders(5)

ex = @processNLExpr  sum{ x[i]^2, i =1:5 } + sin(x[1]*x[2])
sp = compute_hessian_sparsity_IJ(ex)
hfunc = gen_hessian_dense(ex)
exact(x) = [ -x[2]^2*sin(x[1]*x[2]) cos(x[1]*x[2])-x[1]*x[2]*sin(x[1]*x[2]) 0 0 0
            cos(x[1]*x[2])-x[1]*x[2]*sin(x[1]*x[2]) -x[1]^2*sin(x[1]*x[2]) 0 0 0
            0 0 0 0 0
            0 0 0 0 0
            0 0 0 0 0] + diagm(fill(2.0, length(x)))
val = rand(5)

@test_approx_eq hfunc(val) exact(val)
test_sparsity(sp, exact(val))

sparsemat, sparsefunc = gen_hessian_sparse_mat(ex)
sparsefunc(val, sparsemat)
@test_approx_eq sparsemat tril(exact(val))

I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex)
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq to_H(ex, I, J, V, 5) tril(exact(val))

constr = {}
for i in 1:5
    push!(constr,@processNLExpr x[i]^2)
    #println(constr[end])
end

I,J, sparsefunc_color = gen_hessian_sparse_color_parametric(constr[1])
V = zeros(length(I))
exact(i) = (x = zeros(5); x[i] = 2; diagm(x))
for i in 1:5
    sparsefunc_color(val, V, constr[i])
    @test_approx_eq to_H(constr[i], I, J, V, 5) exact(i)
end


# test linear expressions
x,y = placeholders(2)
ex = @processNLExpr 2x + y
I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex)
@assert length(I) == length(J) == 0

println("Passed tests")
