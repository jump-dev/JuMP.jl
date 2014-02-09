using Base.Test
using ReverseDiffSparse

x,y,z,q = placeholders(4)

ex = @process sin(x*y) + exp(z+2q)
sp = compute_hessian_sparsity(ex)
hfunc = gen_hessian_dense(ex)
val = [3.0, 4.0, 5.0, 6.0]
exact(x,y,z,q) = [ -y^2*sin(x*y) cos(x*y)-x*y*sin(x*y) 0 0
                   cos(x*y)-x*y*sin(x*y) -x^2*sin(x*y) 0 0
                   0 0 exp(z+2q) 2*exp(z+2q)
                   0 0 2*exp(z+2q) 4*exp(z+2q) ]
@test_approx_eq normfro(hfunc(val) - exact(val...)) 0.0

sparsemat, sparsefunc = gen_hessian_sparse(ex)
sparsefunc(val, sparsemat)
@test_approx_eq normfro(sparsemat - tril(exact(val...))) 0.0

