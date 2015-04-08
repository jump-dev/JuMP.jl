using Base.Test
using ReverseDiffSparse

function test_sparsity(sp, H)
    I,J = sp
    Hsp = sparse(I,J,ones(length(I)),size(H,1),size(H,2))
    H = sparse(tril(H))
    @test all( Hsp.colptr .== H.colptr)
    @test all( Hsp.rowval .== H.rowval)
end

x,y,z,q = placeholders(4)

ex = @processNLExpr sin(x*y) + exp(z+2q)
sp = compute_hessian_sparsity_IJ(ex,4)
hfunc = gen_hessian_dense(ex)
val = [3.0, 4.0, 5.0, 6.0]
exact(x,y,z,q) = [ -y^2*sin(x*y) cos(x*y)-x*y*sin(x*y) 0 0
                   cos(x*y)-x*y*sin(x*y) -x^2*sin(x*y) 0 0
                   0 0 exp(z+2q) 2*exp(z+2q)
                   0 0 2*exp(z+2q) 4*exp(z+2q) ]
@test_approx_eq hfunc(val) exact(val...)
test_sparsity(sp, exact(val...))

sparsemat, sparsefunc = gen_hessian_sparse_mat(ex,4)
sparsefunc(val, sparsemat)
@test_approx_eq sparsemat tril(exact(val...))

I,J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex, 4)
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq to_H(ex, I, J, V, 4) tril(exact(val...))

x = placeholders(5)

ex = @processNLExpr  sum{ x[i]^2, i =1:5 } + sin(x[1]*x[2])
sp = compute_hessian_sparsity_IJ(ex,5)
hfunc = gen_hessian_dense(ex)
exact(x) = [ -x[2]^2*sin(x[1]*x[2]) cos(x[1]*x[2])-x[1]*x[2]*sin(x[1]*x[2]) 0 0 0
            cos(x[1]*x[2])-x[1]*x[2]*sin(x[1]*x[2]) -x[1]^2*sin(x[1]*x[2]) 0 0 0
            0 0 0 0 0
            0 0 0 0 0
            0 0 0 0 0] + diagm(fill(2.0, length(x)))
val = rand(5)

@test_approx_eq hfunc(val) exact(val)
test_sparsity(sp, exact(val))

sparsemat, sparsefunc = gen_hessian_sparse_mat(ex,5)
sparsefunc(val, sparsemat)
@test_approx_eq sparsemat tril(exact(val))

I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex,5)
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq to_H(ex, I, J, V, 5) tril(exact(val))

# sparsity detection with constants
a = 10
b = 20
ex = @processNLExpr  sum{ a*b*x[i], i =1:5 } + sin(x[1]*x[2])
exact(x) = [ -x[2]^2*sin(x[1]*x[2]) cos(x[1]*x[2])-x[1]*x[2]*sin(x[1]*x[2]) 0 0 0
            cos(x[1]*x[2])-x[1]*x[2]*sin(x[1]*x[2]) -x[1]^2*sin(x[1]*x[2]) 0 0 0
            0 0 0 0 0
            0 0 0 0 0
            0 0 0 0 0]
sp = compute_hessian_sparsity_IJ(ex,5)
test_sparsity(sp, exact(val))


# Expr list
exlist = ExprList()
for i in 1:5
    push!(exlist,@processNLExpr x[i]^3/6)
end

I,J = prep_sparse_hessians(exlist,5)
prep_expression_output(exlist)
V = zeros(length(I))
lambda = rand(5)
eval_hess!(V, exlist, val, lambda)
@test_approx_eq sparse(I,J,V) diagm(lambda.*val)
@test to_flat_expr(exlist, 4) == :(x[4]^3/6)

# irregular indices
exlist = ExprList()
idx = Any[[1,1],[1,2]]
for i in 1:2
    push!(exlist,@processNLExpr sum{x[k]^2/2, k in idx[i]})
end
I,J = prep_sparse_hessians(exlist,2)
V = zeros(length(I))
lambda = [1.0,2.0]
eval_hess!(V, exlist, val, lambda)
@test_approx_eq sparse(I,J,V) diagm([4.,2.])

# test linear expressions
x,y = placeholders(2)
ex = @processNLExpr 2x + y
prepare_indexmap(ex)
I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex,2)
@assert length(I) == length(J) == 0

# constant expressions
a = 10
ex = @processNLExpr (1/a+a)*x^2*y
prepare_indexmap(ex)
I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex,2)
exact(x,y) = [2y*(1/a+a) 0; 2x*(1/a+a) 0]
val = [4.5,2.3]
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq to_H(ex, I, J, V, 2) tril(exact(val...))

# subexpressions
subex = @parametricExpr (1/a+a)*x^2*y
ex = @processNLExpr subex
prepare_indexmap(ex)
I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex,2)
exact(x,y) = [2y*(1/a+a) 0; 2x*(1/a+a) 0]
val = [4.5,2.3]
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq to_H(ex, I, J, V, 2) tril(exact(val...))

# prod{}
x = placeholders(2)
ex = @processNLExpr prod{x[i], i = 1:2}
prepare_indexmap(ex)
I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex,2)
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq to_H(ex, I, J, V, 2) tril([ 0.0 1.0; 1.0 0.0 ])

x = placeholders(3)
val = [4.5,2.3,6.5]
ex = @processNLExpr prod{x[i], i = 1:3}
prepare_indexmap(ex)
I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex,3)
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq full(to_H(ex, I, J, V, 3)) [0 0 0; val[3] 0 0; val[2] val[1] 0]

x = placeholders(4)
val = [4.5,2.3,6.5,3.2]
ex = @processNLExpr prod{x[i], i = 1:4}
prepare_indexmap(ex)
I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex,4)
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq full(to_H(ex, I, J, V, 4)) [0 0 0 0; val[3]*val[4] 0 0 0; val[2]*val[4] val[1]*val[4] 0 0; val[2]*val[3] val[1]*val[3] val[1]*val[2] 0]

# hs071
x = placeholders(4)
ex = @processNLExpr x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3]
prepare_indexmap(ex)
I, J, sparsefunc_color = gen_hessian_sparse_color_parametric(ex,4)
exact(x) = [
2x[4] 0 0 0;
x[4] 0 0 0;
x[4] 0 0 0;
2x[1]+x[2]+x[3] x[1] x[1] 0]
val = [1.0,5.0,5.0,1.0]
V = zeros(length(I))
sparsefunc_color(val, V, ex)
@test_approx_eq full(to_H(ex, I, J, V, 4)) exact(val)

exlist = ExprList()
push!(exlist, @processNLExpr x[1]*x[2]*x[3]*x[4])
#push!(exlist, @processNLExpr sum{x[i]^2,i=1:4})

I,J = prep_sparse_hessians(exlist,4)
V = zeros(length(I))
lambda = 2.0
eval_hess!(V, exlist, val, lambda)
@test_approx_eq full(sparse(I,J,V)) lambda*[0 0 0 0; val[3]*val[4] 0 0 0; val[2]*val[4] val[1]*val[4] 0 0; val[2]*val[3] val[1]*val[3] val[1]*val[2] 0]


# ifelse
ex = @processNLExpr ifelse(x[1] >= x[2], x[1],x[2])
I,J = compute_hessian_sparsity_IJ(ex,4)
@test length(I) == 0 && length(J) == 0

ex = @processNLExpr x[1]*ifelse(x[1] >= x[2], x[1],x[2])
sp = compute_hessian_sparsity_IJ(ex,4)
test_sparsity(sp, [1.0 1.0; 1.0 1.0])

ex = @processNLExpr x[1]*ifelse(x[1] >= x[2] && true, x[1],x[2])
sp = compute_hessian_sparsity_IJ(ex,4)
test_sparsity(sp, [1.0 1.0; 1.0 1.0])

println("Passed tests")
