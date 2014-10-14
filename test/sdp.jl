#############
# Expressions
#############

# Test: * SDP bounds
#       * Changing objective sense
#       * getObjectiveValue
#       * getValue(d::MatrixFuncVar)
m = Model()
@defSDPVar(m, 0 <= X[3] <= 1/2*eye(3,3))
@defSDPVar(m, -ones(5,5) <= Y[5] <= 2*ones(5,5))
@defSDPVar(m, Z[4] <= ones(4,4))

addConstraint(m, trace(X) == 1)
addConstraint(m, trace(Y) == 3)
addConstraint(m, trace(Z) == -1)
setObjective(m, :Max, X[1,2] + Y[1,2] + Z[1,2])
solve(m)

@test_approx_eq_eps getValue(X[1,2]) 0.25 1e-3
@test_approx_eq_eps getValue(Y[1,2]) 0.6 1e-3
@test_approx_eq_eps getValue(Z[1,2]) 3.5 1e-3
@test_approx_eq_eps getObjectiveValue(m) 4.35 1e-3

setObjective(m, :Min, X[1,2] + Y[1,2] + Z[1,2])
solve(m)

@test_approx_eq_eps getValue(X[1,2]) -0.25 1e-3
@test_approx_eq_eps getValue(Y[1,2]) 0.6 1e-3
@test_approx_eq_eps getValue(Z[1,2]) -1.5 1e-3
@test_approx_eq_eps getObjectiveValue(m) -1.15 1e-3

# Test: * getValue(d::MatrixFuncExpr)
#       * getValue(d::MatrixExpr)
#       * getValue(d::SDPVar)
A = [0.25 -0.25 0.0; -0.25 0.25 0.0; 0.0 0.0 0.5]
@test_approx_eq_eps norm(getValue(X)-A) 0.0 1e-3
@test_approx_eq_eps norm(getValue(Y)-0.6*ones(5,5)) 0.0 1e-3
B = ones(4,4)
B[1:2,1:2] = -1.5*ones(2,2)
@test_approx_eq_eps norm(getValue(Z) - B) 0.0 1e-3
@test_approx_eq_eps norm(getValue([2X eye(3,4); eye(4,3) -3Z+ones(4,4)] + eye(7,7)) -
        ([2A eye(3,4); eye(4,3) -3B+ones(4,4)] + eye(7,7))) 0.0 1e-3

@test_approx_eq_eps getValue(2X[1,1]+3.5) (2A[1,1]+3.5) 1e-3
@test_approx_eq_eps getValue(-X[1,3]+2Z[2,4]-1.0) (-A[1,3]+2B[2,4]-1.0) 1e-3

# Test: * getLower(x::SDPVar), getUpper
@test getLower(X) == 0.0
@test getUpper(X) == 1/2*eye(3,3)
@test getLower(Y) == -ones(5,5)
@test getUpper(Y) == 2*ones(5,5)
@test getLower(Z) == -Inf
@test getUpper(Z) == ones(4,4)

# Test * getName(x::SDPVar), setName
@test getName(X) == "X"
setName(X, "my new name")
@test getName(X) == "my new name"

# Test: * @defSDPVar
@test_throws Exception @defSDPVar(m, psd[2] <= rand(2,2))
@test_throws Exception @defSDPVar(m, -Inf <= unbounded[3] <= Inf)
@test_throws Exception @defSDPVar(m, ones(4,4) <= constant[4] <= ones(4,4))
@test_throws Exception @defSDPVar(m, rand(5,5) <= nonsymmetric[5] <= rand(5,5))
@test_throws Exception @defSDPVar(m, -1.0 <= nonzero[6] <= 1.0)

###########
# 
###########
# transpose
# size
# issym
# getindex
# hcat, vcat, hvcat 

###########
# Operators
###########
const 𝕀 = UniformScaling(1)

# Number--Number
# Number--Variable
# Number--AffExpr
# Number--DualExpr
let 
    m = Model()
    @defVar(m,x)
    dual = ones(3,3)*x - 2eye(3,3)
    @test_throws MethodError 2 + dual
    @test_throws MethodError 2 - dual
    expr3 = 2 * dual
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {x})
    @test isequal(expr3.coeffs, {2ones(3,3)})
    @test expr3.constant == -4eye(3,3)
    @test_throws Exception 2 / dual
    @test_throws Exception (2 >= dual)
    con1 = (eye(3,3) >= dual)
    @test isequal(con1.terms, eye(3,3)-dual)
    @test con1.sense == :>=
    con2 = (eye(3,3) <= dual)
    @test isequal(con2.terms, eye(3,3)-dual)
    @test con2.sense == :<=
    con3 = (eye(3,3) == dual)
    @test isequal(con3.terms, eye(3,3)-dual)
    @test con3.sense == :(==)
end
# Number--AbstractArray
# Number--MatrixVar
let 
    m = Model()
    @defMatrixVar(m,X[2,2])
    @test_throws MethodError 2 + X
    @test_throws MethodError 2 - X
    expr3 = -2X
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, {X})
    @test isequal(expr3.pre, {-2𝕀})
    @test isequal(expr3.post, {𝕀})
    @test isequal(expr3.constant, spzeros(2,2))
    @test_throws MethodError 2/X
    @test_throws Exception (2 >= X)
    con1 = (0 >= X)
    @test isequal(con1.terms, -X)
    @test con1.sense == :>=
    con2 = (0 <= X)
    @test isequal(con2.terms, -X)
    @test con2.sense == :<=
    con3 = (0 == X)
    @test isequal(con3.terms, -X)
    @test con3.sense == :(==)
end
# Number--SDPVar
let 
    m = Model()
    @defSDPVar(m,X[2,2])
    @test_throws MethodError 2 + X
    @test_throws MethodError 2 - X
    expr3 = -2X
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, {X})
    @test isequal(expr3.pre, {-2𝕀})
    @test isequal(expr3.post, {𝕀})
    @test isequal(expr3.constant, spzeros(2,2))
    @test_throws MethodError 2/X
    @test_throws Exception (2 >= X)
    con1 = (0 >= X)
    @test isequal(con1.terms, -X)
    @test con1.sense == :>=
    con2 = (0 <= X)
    @test isequal(con2.terms, -X)
    @test con2.sense == :<=
    con3 = (0 == X)
    @test isequal(con3.terms, -X)
    @test con3.sense == :(==)
end
# Number--MatrixExpr
let 
    m = Model()
    @defSDPVar(m,X[2,2])
    expr = [-X ones(2); ones(1,2) -3]
    @test_throws MethodError 2 + expr
    @test_throws MethodError 2 - expr
    expr3 = -2expr
    @test isa(expr3, MatrixExpr)
    @test length(expr3.elem) == 1
    @test isequal(expr3.elem[1][1,1], -X[1,1]+0) # lift to ScalarExpr
    @test isequal(expr3.elem[1][1,2], -X[1,2]+0)
    @test isequal(expr3.elem[1][2,1], -X[2,1]+0)
    @test isequal(expr3.elem[1][2,2], -X[2,2]+0)
    @test isequal(expr3.elem[1][1,3], 1)
    @test isequal(expr3.elem[1][2,3], 1)
    @test isequal(expr3.elem[1][3,1], 1)
    @test isequal(expr3.elem[1][3,2], 1)
    @test isequal(expr3.elem[1][3,3], -3)
    @test isequal(expr3.pre,  {-2𝕀})
    @test isequal(expr3.post, {𝕀})
    @test isequal(expr3.constant, spzeros(3,3))
    @test_throws MethodError 2/expr
    @test_throws MethodError 2/X
    @test_throws Exception (2 >= expr)
    con1 = (0 >= expr)
    @test isequal(con1.terms, -expr)
    @test con1.sense == :>=
    con2 = (0 <= expr)
    @test isequal(con2.terms, expr)
    @test con2.sense == :>=
    con3 = (0 == expr)
    @test isequal(con3.terms, expr)
    @test con3.sense == :(==)
end
# Number--MatrixFuncVar
let 
    m = Model()
    @defSDPVar(m,X[2,2])
    expr1 = 2 + trace(X)
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, AffExpr(2.0))
    @test isequal(expr1.matvars, {trace(X)})
    @test isequal(expr1.normexpr, NormExpr())
    expr2 = 2 - trace(X)
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, AffExpr(2.0))
    @test isequal(expr2.matvars, {-trace(X)})
    @test isequal(expr2.normexpr, NormExpr())
    expr3 = 2 * trace(X)
    @test isa(expr3, MatrixFuncVar)
    @test isequal(expr3.expr, 2X)
    @test expr3.func == :trace
    @test_throws MethodError 2 / trace(X)
    con1 = (1 >= trace(X))
    @test isequal(con1.terms, 1-trace(X))
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (1 <= trace(X))
    @test isequal(con2.terms, 1-trace(X))
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (1 == trace(X))
    @test isequal(con3.terms, 1-trace(X))
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# Number--NormExpr
let 
    m = Model()
    @defMatrixVar(m,X[2])
    expr = norm(2X - ones(2))
    expr1 = 2 + expr
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, AffExpr(2.0))
    @test isempty(expr1.matvars)
    @test isequal(expr1.normexpr, expr)
    expr2 = 2 - expr
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, AffExpr(2.0))
    @test isempty(expr2.matvars)
    @test isequal(expr2.normexpr, -expr)
    expr3 = 2 * expr
    @test isa(expr3, NormExpr)
    @test isequal(expr3.vars, {2X-ones(2)})
    @test isequal(expr3.coeffs, [2.0])
    @test isequal(expr3.form, [:norm2])
    @test_throws MethodError 2 / expr
    con1 = (1 >=  expr)
    @test isequal(con1.terms, 1-expr)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (1 <= expr)
    @test isequal(con2.terms, 1-expr)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (1 == expr)
    @test isequal(con3.terms, 1-expr)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# Number--ScalarExpr
let 
    m = Model()
    @defVar(m,x)
    @defSDPVar(m,X[2])
    expr = 2x + 3 - 4trace(X) + norm(2X - ones(2,2))
    expr1 = 2 + expr
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, 2x + 5)
    @test isequal(expr1.matvars, {-4trace(X)})
    @test isequal(expr1.normexpr, norm(2X - ones(2,2)))
    expr2 = 2 - expr
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, -2x - 1)
    @test isequal(expr2.matvars, {4trace(X)})
    @test isequal(expr2.normexpr, -norm(2X - ones(2,2)))
    expr3 = 2 * expr
    @test isa(expr3, ScalarExpr)
    @test isequal(expr3.aff, 4x + 6)
    @test isequal(expr3.matvars, {-8trace(X)})
    @test isequal(expr3.normexpr, 2norm(2X - ones(2,2)))
    @test_throws MethodError 2 / expr
    con1 = (1 >=  expr)
    @test isequal(con1.terms, 1-expr)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (1 <= expr)
    @test isequal(con2.terms, 1-expr)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (1 == expr)
    @test isequal(con3.terms, 1-expr)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end

# Variable--Number
# Variable--Variable
# Variable--AffExpr
# Variable--DualExpr
# Variable--AbstractArray
let
    m = Model()
    @defVar(m,x)
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    @test_throws MethodError x + arr
    @test_throws MethodError x - arr
    expr3 = x * arr
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {x})
    @test isequal(expr3.coeffs, {arr})
    @test isequal(expr3.constant, spzeros(4,4))
    @test_throws MethodError x / arr
end
# Variable--MatrixVar
# Variable--SDPVar
# Variable--MatrixExpr
# Variable--MatrixFuncVar
let
    m = Model()
    @defVar(m,x)
    @defSDPVar(m,X[2])
    mfv = trace(X)
    expr1 = x + mfv
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, AffExpr(x))
    @test isequal(expr1.matvars, {trace(X)})
    @test isequal(expr1.normexpr, NormExpr())
    expr2 = x - mfv
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, AffExpr(x))
    @test isequal(expr2.matvars, {-trace(X)})
    @test isequal(expr2.normexpr, NormExpr())
    @test_throws MethodError x * mfv
    @test_throws MethodError x / mfv
    con1 = (x >=  trace(X))
    @test isequal(con1.terms, x-trace(X))
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (x <=  trace(X))
    @test isequal(con2.terms, x-trace(X))
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (x ==  trace(X))
    @test isequal(con3.terms, x-trace(X))
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# Variable--NormExpr
let
    m = Model()
    @defVar(m,x)
    @defMatrixVar(m,X[2])
    expr = norm(2X - ones(2))
    expr1 = x + expr
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, AffExpr(x))
    @test isempty(expr1.matvars)
    @test isequal(expr1.normexpr, expr)
    expr2 = x - expr
    @test isa(expr1, ScalarExpr)
    @test isequal(expr2.aff, AffExpr(x))
    @test isempty(expr2.matvars)
    @test isequal(expr2.normexpr, -expr)
    @test_throws MethodError x * expr
    @test_throws MethodError x / expr
    con1 = (x >=  expr)
    @test isequal(con1.terms, x-expr)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (x <=  expr)
    @test isequal(con2.terms, x-expr)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (x ==  expr)
    @test isequal(con3.terms, x-expr)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# Variable--ScalarExpr
let
    m = Model()
    @defVar(m,x)
    @defVar(m,y)
    @defSDPVar(m,X[2])
    @defMatrixVar(m,Y[2])
    expr = 3 - x + 2trace(X) - 4norm(2Y+ones(2))
    expr1 = y + expr
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, 3 - x + y)
    @test isequal(expr1.matvars, {2trace(X)})
    @test isequal(expr1.normexpr, -4norm(2Y+ones(2)))
    expr2 = y - expr
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, -3 + x + y)
    @test isequal(expr2.matvars, {-2trace(X)})
    @test isequal(expr2.normexpr, 4norm(2Y+ones(2)))
    @test_throws MethodError x * expr
    @test_throws MethodError x / expr
    con1 = (x >=  expr)
    @test isequal(con1.terms, x-expr)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (x <=  expr)
    @test isequal(con2.terms, x-expr)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (x ==  expr)
    @test isequal(con3.terms, x-expr)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end

# AffExpr--Number
# AffExpr--Variable
# AffExpr--AffExpr
# AffExpr--DualExpr
# AffExpr--AbstractArray
let
    m = Model()
    @defVar(m,x)
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    @test_throws MethodError (x+1) + arr
    @test_throws MethodError (x+1) - arr
    expr3 = (x+1) * arr
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {x})
    @test isequal(expr3.coeffs, {arr})
    @test isequal(expr3.constant, arr)
    @test_throws MethodError x / arr
end
# AffExpr--MatrixVar
# AffExpr--SDPVar
# AffExpr--MatrixExpr
# AffExpr--MatrixFuncVar
let
    m = Model()
    @defVar(m,x)
    @defSDPVar(m,X[2])
    expr1 = (x+1) + (trace(X))
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, x+1)
    @test isequal(expr1.matvars, {trace(X)})
    @test isequal(expr1.normexpr, NormExpr())
    expr2 = (x+1) - (trace(X))
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, x+1)
    @test isequal(expr2.matvars, {-trace(X)})
    @test isequal(expr2.normexpr, NormExpr())
    @test_throws MethodError (x+1) * (trace(X))
    @test_throws MethodError (x+1) / (trace(X))
    con1 = ((x+1) >=  trace(X))
    @test isequal(con1.terms, (x+1)-trace(X))
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = ((x+1) <=  trace(X))
    @test isequal(con2.terms, (x+1)-trace(X))
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = ((x+1) ==  trace(X))
    @test isequal(con3.terms, (x+1)-trace(X))
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# AffExpr--NormExpr
let
    m = Model()
    @defVar(m,x)
    @defMatrixVar(m,Y[2])
    expr1 = (x+1) + norm(2Y-ones(2))
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, x+1)
    @test isempty(expr1.matvars)
    @test isequal(expr1.normexpr, norm(2Y-ones(2)))
    expr2 = (x+1) - norm(2Y-ones(2))
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, x+1)
    @test isempty(expr2.matvars)
    @test isequal(expr2.normexpr, -norm(2Y-ones(2)))
    @test_throws MethodError (x+1) * norm(2Y-ones(2))
    @test_throws MethodError (x+1) / norm(2Y-ones(2))
    con1 = ((x+1) >=  norm(2Y-ones(2)))
    @test isequal(con1.terms, (x+1)-norm(2Y-ones(2)))
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = ((x+1) <=  norm(2Y-ones(2)))
    @test isequal(con2.terms, (x+1)-norm(2Y-ones(2)))
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = ((x+1) ==  norm(2Y-ones(2)))
    @test isequal(con3.terms, (x+1)-norm(2Y-ones(2)))
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# AffExpr--ScalarExpr
let
    m = Model()
    @defVar(m,x)
    @defVar(m,y)
    @defSDPVar(m,X[2])
    @defMatrixVar(m,Y[2])
    expr = (-2 - y - 3trace(X) + 2norm(2Y-ones(2)))
    expr1 = (x+1) + expr
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, x-y-1)
    @test isequal(expr1.matvars, {-3trace(X)})
    @test isequal(expr1.normexpr, 2norm(2Y-ones(2)))
    expr2 = (x+1) - expr
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, x+y+3)
    @test isequal(expr2.matvars, {3trace(X)})
    @test isequal(expr2.normexpr, -2norm(2Y-ones(2)))
    @test_throws MethodError (x+1) * expr
    @test_throws MethodError (x+1) / expr
    con1 = ((x+1) >=  expr)
    @test isequal(con1.terms, (x+1)-expr)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = ((x+1) <=  expr)
    @test isequal(con2.terms, (x+1)-expr)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = ((x+1) ==  expr)
    @test isequal(con3.terms, (x+1)-expr)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end

# DualExpr--Number
# DualExpr--Variable
# DualExpr--AffExpr
# DualExpr--DualExpr
let
    m = Model()
    @defVar(m,x)
    @defVar(m,y)
    dual1 = x*ones(3,3) - 2eye(3,3)
    dual2 = -3eye(3,3)*y + 2ones(3,3)
    expr1 = dual1 + dual2
    @test isa(expr1, DualExpr)
    @test isequal(expr1.vars, {x,y})
    @test isequal(expr1.coeffs, {ones(3,3), -3eye(3,3)})
    @test expr1.constant == -2eye(3,3)+2ones(3,3)
    expr2 = dual1 - dual2
    @test isa(expr2, DualExpr)
    @test isequal(expr2.vars, {x,y})
    @test isequal(expr2.coeffs, {ones(3,3), 3eye(3,3)})
    @test expr2.constant == -2eye(3,3)-2ones(3,3)
    @test_throws MethodError dual1 * dual2
    @test_throws MethodError dual1 / dual2
    con1 = (dual1 >= dual2)
    @test isequal(con1.terms, dual1-dual2)
    @test con1.sense == :>=
    con2 = (dual1 <= dual2)
    @test isequal(con2.terms, dual1-dual2)
    @test con2.sense == :<=
    con3 = (dual1 == dual2)
    @test isequal(con3.terms, dual1-dual2)
    @test con3.sense == :(==)
end
# DualExpr--AbstractArray
let
    m = Model()
    @defVar(m,x)
    @defVar(m,y)
    dual = x*ones(4,4) - 2eye(4,4)
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr1 = dual + arr
    @test isa(expr1, DualExpr)
    @test isequal(expr1.vars, {x})
    @test isequal(expr1.coeffs, {ones(4,4)})
    @test expr1.constant == -2eye(4,4)+arr
    expr2 = dual - arr
    @test isa(expr2, DualExpr)
    @test isequal(expr2.vars, {x})
    @test isequal(expr2.coeffs, {ones(4,4)})
    @test expr2.constant == -2eye(4,4)-arr
    expr3 = dual * arr
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {x})
    @test isequal(expr3.coeffs, {ones(4,4)*arr})
    @test expr3.constant == -2eye(4,4)*arr
    @test_throws MethodError dual / arr
    con1 = (dual >= arr)
    @test isequal(con1.terms, dual-arr)
    @test con1.sense == :>=
    con2 = (dual <= arr)
    @test isequal(con2.terms, dual-arr)
    @test con2.sense == :<=
    con3 = (dual == arr)
    @test isequal(con3.terms, dual-arr)
    @test con3.sense == :(==)
end
# DualExpr--MatrixVar
# DualExpr--SDPVar
# DualExpr--MatrixExpr
# DualExpr--MatrixFuncVar
# DualExpr--NormExpr
# DualExpr--ScalarExpr

# AbstractArray--Number
# AbstractArray--Variable
let
    m = Model()
    @defVar(m,x)
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    @test_throws MethodError arr + x
    @test_throws MethodError arr - x
    expr3 = arr * x
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {x})
    @test isequal(expr3.coeffs, {arr})
    @test isequal(expr3.constant, spzeros(4,4))
    @test_throws MethodError arr / x
end
# AbstractArray--AffExpr
let
    m = Model()
    @defVar(m,x)
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    @test_throws MethodError arr + (x+1)
    @test_throws MethodError arr - (x+1)
    expr3 = arr * (x+1)
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {x})
    @test isequal(expr3.coeffs, {arr})
    @test isequal(expr3.constant, arr)
    @test_throws MethodError arr / x
end
# AbstractArray--DualExpr
let
    m = Model()
    @defVar(m,x)
    @defVar(m,y)
    dual = x*ones(4,4) - 2eye(4,4)
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr1 = arr + dual
    @test isa(expr1, DualExpr)
    @test isequal(expr1.vars, {x})
    @test isequal(expr1.coeffs, {ones(4,4)})
    @test expr1.constant == arr - 2eye(4,4)
    expr2 = arr - dual
    @test isa(expr2, DualExpr)
    @test isequal(expr2.vars, {x})
    @test isequal(expr2.coeffs, {-ones(4,4)})
    @test expr2.constant == arr + 2eye(4,4)
    expr3 = arr * dual
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {x})
    @test isequal(expr3.coeffs, {arr*ones(4,4)})
    @test expr3.constant == arr*(-2eye(4,4))
    @test_throws MethodError arr / dual
    con1 = (arr >= dual)
    @test isequal(con1.terms, arr-dual)
    @test con1.sense == :>=
    con2 = (arr <= dual)
    @test isequal(con2.terms, arr-dual)
    @test con2.sense == :<=
    con3 = (arr == dual)
    @test isequal(con3.terms, arr-dual)
    @test con3.sense == :(==)
end
# AbstractArray--AbstractArray
# AbstractArray--MatrixVar
let
    m = Model()
    @defMatrixVar(m, X[4,4])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr1 = arr + X
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X})
    @test isequal(expr1.pre, {𝕀})
    @test isequal(expr1.post, {𝕀})
    @test expr1.constant == arr
    expr2 = arr - X
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X})
    @test isequal(expr2.pre, {-𝕀})
    @test isequal(expr2.post, {𝕀})
    @test expr2.constant == arr
    expr3 = arr * X
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, {X})
    @test isequal(expr3.pre, {arr})
    @test isequal(expr3.post, {𝕀})
    @test expr3.constant == spzeros(4,4)
    @test_throws MethodError arr / X
    con1 = (arr >= X)
    @test isequal(con1.terms, arr-X)
    @test con1.sense == :>=
    con2 = (arr <= X)
    @test isequal(con2.terms, arr-X)
    @test con2.sense == :<=
    con3 = (arr == X)
    @test isequal(con3.terms, arr-X)
    @test con3.sense == :(==)
end
# AbstractArray--SDPVar
let
    m = Model()
    @defSDPVar(m, X[4,4])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr1 = arr + X
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X})
    @test isequal(expr1.pre, {𝕀})
    @test isequal(expr1.post, {𝕀})
    @test expr1.constant == arr
    expr2 = arr - X
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X})
    @test isequal(expr2.pre, {-𝕀})
    @test isequal(expr2.post, {𝕀})
    @test expr2.constant == arr
    expr3 = arr * X
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, {X})
    @test isequal(expr3.pre, {arr})
    @test isequal(expr3.post, {𝕀})
    @test expr3.constant == spzeros(4,4)
    @test_throws MethodError arr / X
    con1 = (arr >= X)
    @test isequal(con1.terms, arr-X)
    @test con1.sense == :>=
    con2 = (arr <= X)
    @test isequal(con2.terms, arr-X)
    @test con2.sense == :<=
    con3 = (arr == X)
    @test isequal(con3.terms, arr-X)
    @test con3.sense == :(==)
end
# AbstractArray--MatrixExpr
let
    m = Model()
    @defMatrixVar(m, X[3,3])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr = [-2 ones(3)'; ones(3) -3X]
    expr1 = arr + expr
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, expr.elem)
    @test isequal(expr1.pre, {𝕀})
    @test isequal(expr1.post, {𝕀})
    @test expr1.constant == arr + expr.constant
    expr2 = arr - expr
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, expr.elem)
    @test isequal(expr2.pre, {-𝕀})
    @test isequal(expr2.post, {𝕀})
    @test expr2.constant == arr - expr.constant
    expr3 = arr * expr
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, expr.elem)
    @test isequal(expr3.pre, {arr})
    @test isequal(expr3.post, {𝕀})
    @test expr3.constant == arr*expr.constant
    @test_throws MethodError arr / expr
    con1 = (arr >= expr)
    @test isequal(con1.terms, arr-expr)
    @test con1.sense == :>=
    con2 = (arr <= expr)
    @test isequal(con2.terms, arr-expr)
    @test con2.sense == :<=
    con3 = (arr == expr)
    @test isequal(con3.terms, arr-expr)
    @test con3.sense == :(==)
end
# AbstractArray--MatrixFuncVar
let
    m = Model()
    @defSDPVar(m, X[4])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    @test_throws MethodError arr + trace(X)
    @test_throws MethodError arr - trace(X)
    expr3 = arr * trace(X)
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {trace(X)})
    @test isequal(expr3.coeffs, {arr})
    @test expr3.constant == spzeros(4,4)
    @test_throws MethodError arr / trace(X)
end
# AbstractArray--NormExpr
# AbstractArray--ScalarExpr
let
    m = Model()
    @defVar(m,x)
    @defSDPVar(m, X[2])
    @defMatrixVar(m, Y[2])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr = 3 - x + 2trace(X)
    @test_throws MethodError arr + expr
    @test_throws MethodError arr - expr
    expr3 = arr * expr
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {2trace(X),x})
    @test length(expr3.coeffs) == 2
    @test expr3.coeffs[1] ==  arr
    @test expr3.coeffs[2] == -arr
    @test expr3.constant == 3arr
    nexpr = expr + vecnorm(Y)
    @test_throws Exception arr * nexpr
    @test_throws MethodError arr / trace(X)
end

# MatrixVar--Number
let 
    m = Model()
    @defMatrixVar(m,X[2,2])
    @test_throws MethodError X + 2
    @test_throws MethodError X - 2
    expr3 = X*(-2)
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, {X})
    @test isequal(expr3.pre, {-2𝕀})
    @test isequal(expr3.post, {𝕀})
    @test isequal(expr3.constant, spzeros(2,2))
    expr4 = X / (-2)
    @test isa(expr4, MatrixExpr)
    @test isequal(expr4.elem, {X})
    @test isequal(expr4.pre, {-0.5𝕀})
    @test isequal(expr4.post, {𝕀})
    @test isequal(expr4.constant, spzeros(2,2))
    @test_throws Exception (X >= 2)
    con1 = (X >= 0)
    @test isequal(con1.terms, 1X)
    @test con1.sense == :>=
    con2 = (X <= 0)
    @test isequal(con2.terms, 1X)
    @test con2.sense == :<=
    con3 = (X == 0)
    @test isequal(con3.terms, 1X)
    @test con3.sense == :(==)
end
# MatrixVar--Variable
# MatrixVar--AffExpr
# MatrixVar--DualExpr
# MatrixVar--AbstractArray
let
    m = Model()
    @defMatrixVar(m, X[4,4])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr1 = X + arr
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X})
    @test isequal(expr1.pre, {𝕀})
    @test isequal(expr1.post, {𝕀})
    @test expr1.constant == arr
    expr2 = X - arr
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X})
    @test isequal(expr2.pre, {𝕀})
    @test isequal(expr2.post, {𝕀})
    @test expr2.constant == -arr
    expr3 = X * arr
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, {X})
    @test isequal(expr3.pre, {𝕀})
    @test isequal(expr3.post, {arr})
    @test expr3.constant == spzeros(4,4)
    @test_throws MethodError X / arr
    con1 = (X >= arr)
    @test isequal(con1.terms, X-arr)
    @test con1.sense == :>=
    con2 = (X <= arr)
    @test isequal(con2.terms, X-arr)
    @test con2.sense == :<=
    con3 = (X == arr)
    @test isequal(con3.terms, X-arr)
    @test con3.sense == :(==)
end
# MatrixVar--MatrixVar
let
    m = Model()
    @defMatrixVar(m, X[2])
    @defMatrixVar(m, Y[2])
    expr1 = X + Y
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X,Y})
    @test isequal(expr1.pre, {𝕀,𝕀})
    @test isequal(expr1.post, {𝕀,𝕀})
    @test expr1.constant == spzeros(2,1)
    expr2 = X - Y
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X,Y})
    @test isequal(expr2.pre, {𝕀,-𝕀})
    @test isequal(expr2.post, {𝕀,𝕀})
    @test expr2.constant == spzeros(2,1)
    @defMatrixVar(m, Z[3,2])
    @defMatrixVar(m, W[2,3])
    # expr3 = Z * W
    # @test isa(expr3, MatrixExpr)
    # for i in 1:3, j in 1:3
    #     @test isequal(expr3[i,j], Z[i,j]*W[j,i])
    # end
    @test_throws MethodError X / Y
    con1 = (X >= Y)
    @test isequal(con1.terms, X-Y)
    @test con1.sense == :>=
    con2 = (X <= Y)
    @test isequal(con2.terms, X-Y)
    @test con2.sense == :<=
    con3 = (X == Y)
    @test isequal(con3.terms, X-Y)
    @test con3.sense == :(==)
end
# MatrixVar--SDPVar
let
    m = Model()
    @defMatrixVar(m, X[2,2])
    @defSDPVar(m, Y[2])
    expr1 = X + Y
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X,Y})
    @test isequal(expr1.pre, {𝕀,𝕀})
    @test isequal(expr1.post, {𝕀,𝕀})
    @test expr1.constant == spzeros(2,2)
    expr2 = X - Y
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X,Y})
    @test isequal(expr2.pre, {𝕀,-𝕀})
    @test isequal(expr2.post, {𝕀,𝕀})
    @test expr2.constant == spzeros(2,2)
    @defMatrixVar(m, Z[3,2])
    @defMatrixVar(m, W[2,3])
    # expr3 = Z * W
    # @test isa(expr3, MatrixExpr)
    # for i in 1:3, j in 1:3
    #     @test isequal(expr3[i,j], Z[i,j]*W[j,i])
    # end
    @test_throws MethodError X / Y
    con1 = (X >= Y)
    @test isequal(con1.terms, X-Y)
    @test con1.sense == :>=
    con2 = (X <= Y)
    @test isequal(con2.terms, X-Y)
    @test con2.sense == :<=
    con3 = (X == Y)
    @test isequal(con3.terms, X-Y)
    @test con3.sense == :(==)
end
# MatrixVar--MatrixExpr
let
    m = Model()
    @defMatrixVar(m, X[2,2])
    @defMatrixVar(m, Y[2,2])
    expr = 2Y - eye(2,2)
    expr1 = X + expr
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X,Y})
    @test isequal(expr1.pre, {𝕀,2𝕀})
    @test isequal(expr1.post, {𝕀,𝕀})
    @test expr1.constant == -eye(2,2)
    expr2 = X - expr
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X,Y})
    @test isequal(expr2.pre, {𝕀,-2𝕀})
    @test isequal(expr2.post, {𝕀,𝕀})
    @test expr2.constant == eye(2,2)
    @defMatrixVar(m, Z[3,2])
    @defMatrixVar(m, W[2,3])
    # expr3 = Z * W
    # @test isa(expr3, MatrixExpr)
    # for i in 1:3, j in 1:3
    #     @test isequal(expr3[i,j], Z[i,j]*W[j,i])
    # end
    @test_throws MethodError X / expr
end
# MatrixVar--MatrixFuncVar
# MatrixVar--NormExpr
# MatrixVar--ScalarExpr

# SDPVar--Number
let 
    m = Model()
    @defSDPVar(m,X[2,2])
    @test_throws MethodError X + 2
    @test_throws MethodError X - 2
    expr3 = X*(-2)
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, {X})
    @test isequal(expr3.pre, {-2𝕀})
    @test isequal(expr3.post, {𝕀})
    @test isequal(expr3.constant, spzeros(2,2))
    expr4 = X / (-2)
    @test isa(expr4, MatrixExpr)
    @test isequal(expr4.elem, {X})
    @test isequal(expr4.pre, {-0.5𝕀})
    @test isequal(expr4.post, {𝕀})
    @test isequal(expr4.constant, spzeros(2,2))
    @test_throws Exception (X >= 2)
    con1 = (X >= 0)
    @test isequal(con1.terms, 1X)
    @test con1.sense == :>=
    con2 = (X <= 0)
    @test isequal(con2.terms, 1X)
    @test con2.sense == :<=
    con3 = (X == 0)
    @test isequal(con3.terms, 1X)
    @test con3.sense == :(==)
end
# SDPVar--Variable
# SDPVar--AffExpr
# SDPVar--DualExpr
# SDPVar--AbstractArray
let
    m = Model()
    @defSDPVar(m, X[4,4])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr1 = X + arr
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X})
    @test isequal(expr1.pre, {𝕀})
    @test isequal(expr1.post, {𝕀})
    @test expr1.constant == arr
    expr2 = X - arr
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X})
    @test isequal(expr2.pre, {𝕀})
    @test isequal(expr2.post, {𝕀})
    @test expr2.constant == -arr
    expr3 = X * arr
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, {X})
    @test isequal(expr3.pre, {𝕀})
    @test isequal(expr3.post, {arr})
    @test expr3.constant == spzeros(4,4)
    @test_throws MethodError X / arr
    con1 = (X >= arr)
    @test isequal(con1.terms, X-arr)
    @test con1.sense == :>=
    con2 = (X <= arr)
    @test isequal(con2.terms, X-arr)
    @test con2.sense == :<=
    con3 = (X == arr)
    @test isequal(con3.terms, X-arr)
    @test con3.sense == :(==)
end
# SDPVar--MatrixVar
let
    m = Model()
    @defSDPVar(m, X[2,2])
    @defMatrixVar(m, Y[2,2])
    expr1 = X + Y
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X,Y})
    @test isequal(expr1.pre, {𝕀,𝕀})
    @test isequal(expr1.post, {𝕀,𝕀})
    @test expr1.constant == spzeros(2,2)
    expr2 = X - Y
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X,Y})
    @test isequal(expr2.pre, {𝕀,-𝕀})
    @test isequal(expr2.post, {𝕀,𝕀})
    @test expr2.constant == spzeros(2,2)
    @defMatrixVar(m, Z[3,2])
    @defMatrixVar(m, W[2,3])
    # expr3 = Z * W
    # @test isa(expr3, MatrixExpr)
    # for i in 1:3, j in 1:3
    #     @test isequal(expr3[i,j], Z[i,j]*W[j,i])
    # end
    @test_throws MethodError X / Y
end
# SDPVar--SDPVar
let
    m = Model()
    @defSDPVar(m, X[2,2])
    @defSDPVar(m, Y[2,2])
    expr1 = X + Y
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X,Y})
    @test isequal(expr1.pre, {𝕀,𝕀})
    @test isequal(expr1.post, {𝕀,𝕀})
    @test expr1.constant == spzeros(2,2)
    expr2 = X - Y
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X,Y})
    @test isequal(expr2.pre, {𝕀,-𝕀})
    @test isequal(expr2.post, {𝕀,𝕀})
    @test expr2.constant == spzeros(2,2)
    @defMatrixVar(m, Z[3,2])
    @defMatrixVar(m, W[2,3])
    # expr3 = Z * W
    # @test isa(expr3, MatrixExpr)
    # for i in 1:3, j in 1:3
    #     @test isequal(expr3[i,j], Z[i,j]*W[j,i])
    # end
    @test_throws MethodError X / Y
    con1 = (X >= Y)
    @test isequal(con1.terms, X-Y)
    @test con1.sense == :>=
    con2 = (X <= Y)
    @test isequal(con2.terms, X-Y)
    @test con2.sense == :<=
    con3 = (X == Y)
    @test isequal(con3.terms, X-Y)
    @test con3.sense == :(==)
end
# SDPVar--MatrixExpr
let
    m = Model()
    @defSDPVar(m, X[2,2])
    @defMatrixVar(m, Y[2,2])
    expr = 2Y - eye(2,2)
    expr1 = X + expr
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X,Y})
    @test isequal(expr1.pre, {𝕀,2𝕀})
    @test isequal(expr1.post, {𝕀,𝕀})
    @test expr1.constant == -eye(2,2)
    expr2 = X - expr
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X,Y})
    @test isequal(expr2.pre, {𝕀,-2𝕀})
    @test isequal(expr2.post, {𝕀,𝕀})
    @test expr2.constant == eye(2,2)
    @defMatrixVar(m, Z[3,2])
    @defMatrixVar(m, W[2,3])
    # expr3 = Z * W
    # @test isa(expr3, MatrixExpr)
    # for i in 1:3, j in 1:3
    #     @test isequal(expr3[i,j], Z[i,j]*W[j,i])
    # end
    @test_throws MethodError X / expr
    con1 = (X >= expr)
    @test isequal(con1.terms, X-expr)
    @test con1.sense == :>=
    con2 = (X <= expr)
    @test isequal(con2.terms, X-expr)
    @test con2.sense == :<=
    con3 = (X == expr)
    @test isequal(con3.terms, X-expr)
    @test con3.sense == :(==)
end
# SDPVar--MatrixFuncVar
# SDPVar--NormExpr
# SDPVar--ScalarExpr

# MatrixExpr--Number
let 
    m = Model()
    @defSDPVar(m,X[2,2])
    expr = [-X ones(2); ones(1,2) -3]
    @test_throws MethodError expr + 2
    @test_throws MethodError expr - 2
    expr3 = expr*(-2)
    @test isa(expr3, MatrixExpr)
    @test length(expr3.elem) == 1
    @test isequal(expr3.elem[1][1,1], -X[1,1]+0) # lift to ScalarExpr
    @test isequal(expr3.elem[1][1,2], -X[1,2]+0)
    @test isequal(expr3.elem[1][2,1], -X[2,1]+0)
    @test isequal(expr3.elem[1][2,2], -X[2,2]+0)
    @test isequal(expr3.elem[1][1,3], 1)
    @test isequal(expr3.elem[1][2,3], 1)
    @test isequal(expr3.elem[1][3,1], 1)
    @test isequal(expr3.elem[1][3,2], 1)
    @test isequal(expr3.elem[1][3,3], -3)
    @test isequal(expr3.pre,  {-2𝕀})
    @test isequal(expr3.post, {𝕀})
    @test isequal(expr3.constant, spzeros(3,3))
    expr4 = expr / -2
    @test isa(expr4, MatrixExpr)
    @test length(expr4.elem) == 1
    @test isequal(expr4.elem[1][1,1], -X[1,1]+0) # lift to ScalarExpr
    @test isequal(expr4.elem[1][1,2], -X[1,2]+0)
    @test isequal(expr4.elem[1][2,1], -X[2,1]+0)
    @test isequal(expr4.elem[1][2,2], -X[2,2]+0)
    @test isequal(expr4.elem[1][1,3], 1)
    @test isequal(expr4.elem[1][2,3], 1)
    @test isequal(expr4.elem[1][3,1], 1)
    @test isequal(expr4.elem[1][3,2], 1)
    @test isequal(expr4.elem[1][3,3], -3)
    @test isequal(expr4.pre,  {-0.5𝕀})
    @test isequal(expr4.post, {𝕀})
    @test isequal(expr4.constant, spzeros(3,3))
end
# MatrixExpr--Variable
# MatrixExpr--AffExpr
# MatrixExpr--DualExpr
# MatrixExpr--AbstractArray
let
    m = Model()
    @defMatrixVar(m, X[3,3])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr = [-2 ones(3)'; ones(3) -3X]
    expr1 = expr + arr
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, expr.elem)
    @test isequal(expr1.pre, {𝕀})
    @test isequal(expr1.post, {𝕀})
    @test expr1.constant == expr.constant + arr
    expr2 = expr - arr
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, expr.elem)
    @test isequal(expr2.pre, {𝕀})
    @test isequal(expr2.post, {𝕀})
    @test expr2.constant == expr.constant - arr
    expr3 = expr * arr
    @test isa(expr3, MatrixExpr)
    @test isequal(expr3.elem, expr.elem)
    @test isequal(expr3.pre, {𝕀})
    @test isequal(expr3.post, {arr})
    @test expr3.constant == expr.constant*arr
    @test_throws MethodError expr / arr
end
# MatrixExpr--MatrixVar
let
    m = Model()
    @defMatrixVar(m, X[2,2])
    @defMatrixVar(m, Y[2,2])
    expr = 2Y - eye(2,2)
    expr1 = expr + X
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {Y,X})
    @test isequal(expr1.pre, {2𝕀,𝕀})
    @test isequal(expr1.post, {𝕀,𝕀})
    @test expr1.constant == -eye(2,2)
    expr2 = expr - X
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {Y,X})
    @test isequal(expr2.pre, {2𝕀,-𝕀})
    @test isequal(expr2.post, {𝕀,𝕀})
    @test expr2.constant == -eye(2,2)
    @defMatrixVar(m, Z[3,2])
    @defMatrixVar(m, W[2,3])
    # expr3 = Z * W
    # @test isa(expr3, MatrixExpr)
    # for i in 1:3, j in 1:3
    #     @test isequal(expr3[i,j], Z[i,j]*W[j,i])
    # end
    @test_throws MethodError expr / X
end
# MatrixExpr--SDPVar
let
    m = Model()
    @defSDPVar(m, X[2,2])
    @defMatrixVar(m, Y[2,2])
    expr = 2Y - eye(2,2)
    expr1 = expr + X
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {Y,X})
    @test isequal(expr1.pre, {2𝕀,𝕀})
    @test isequal(expr1.post, {𝕀,𝕀})
    @test expr1.constant == -eye(2,2)
    expr2 = expr - X
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {Y,X})
    @test isequal(expr2.pre, {2𝕀,-𝕀})
    @test isequal(expr2.post, {𝕀,𝕀})
    @test expr2.constant == -eye(2,2)
    @defMatrixVar(m, Z[3,2])
    @defMatrixVar(m, W[2,3])
    # expr3 = Z * W
    # @test isa(expr3, MatrixExpr)
    # for i in 1:3, j in 1:3
    #     @test isequal(expr3[i,j], Z[i,j]*W[j,i])
    # end
    @test_throws MethodError expr / X
end
# MatrixExpr--MatrixExpr
let
    m = Model()
    @defMatrixVar(m, X[3,3])
    @defSDPVar(m, Y[3,3])
    mex1 = 2X + eye(3,3)
    mex2 = -Y + 2ones(3,3)
    expr1 = mex1 + mex2
    @test isa(expr1, MatrixExpr)
    @test isequal(expr1.elem, {X,Y})
    @test isequal(expr1.pre, {2𝕀,-𝕀})
    @test isequal(expr1.post, {𝕀,𝕀})
    @test expr1.constant == eye(3,3) + 2ones(3,3)
    expr2 = mex1 - mex2
    @test isa(expr2, MatrixExpr)
    @test isequal(expr2.elem, {X,Y})
    @test isequal(expr2.pre, {2𝕀,𝕀})
    @test isequal(expr2.post, {𝕀,𝕀})
    @test expr2.constant == eye(3,3) - 2ones(3,3)
    @test_throws Exception mex1 * mex2
    @test_throws MethodError mex1 / mex2
end
# MatrixExpr--MatrixFuncVar
# MatrixExpr--NormExpr
# MatrixExpr--ScalarExpr

# MatrixFuncVar--Number
let 
    m = Model()
    @defSDPVar(m,X[2,2])
    expr1 = trace(X) + 2
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, AffExpr(2.0))
    @test isequal(expr1.matvars, {trace(X)})
    @test isequal(expr1.normexpr, NormExpr())
    expr2 = trace(X) - 2
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, AffExpr(-2.0))
    @test isequal(expr2.matvars, {trace(X)})
    @test isequal(expr2.normexpr, NormExpr())
    expr3 = trace(X) * 2
    @test isa(expr3, MatrixFuncVar)
    @test isequal(expr3.expr, 2X)
    @test expr3.func == :trace
    expr4 = trace(X) / 2
    @test isa(expr4, MatrixFuncVar)
    @test isequal(expr4.expr, 0.5X)
    @test expr4.func == :trace
    con1 = (trace(X) >= 1)
    @test isequal(con1.terms, trace(X)-1)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (trace(X) <= 1)
    @test isequal(con2.terms, trace(X)-1)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (trace(X) == 1)
    @test isequal(con3.terms, trace(X)-1)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# MatrixFuncVar--Variable
let
    m = Model()
    @defVar(m,x)
    @defSDPVar(m,X[2])
    mfv = trace(X)
    expr1 = mfv + x
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, AffExpr(x))
    @test isequal(expr1.matvars, {trace(X)})
    @test isequal(expr1.normexpr, NormExpr())
    expr2 = mfv - x
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, -x)
    @test isequal(expr2.matvars, {trace(X)})
    @test isequal(expr2.normexpr, NormExpr())
    @test_throws MethodError mfv * x
    @test_throws MethodError mfv / x
    con1 = (trace(X) >= x)
    @test isequal(con1.terms, trace(X)-x)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (trace(X) <= x)
    @test isequal(con2.terms, trace(X)-x)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (trace(X) == x)
    @test isequal(con3.terms, trace(X)-x)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# MatrixFuncVar--AffExpr
let
    m = Model()
    @defVar(m,x)
    @defSDPVar(m,X[2])
    expr1 = (trace(X)) + (x+1)
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, x+1)
    @test isequal(expr1.matvars, {trace(X)})
    @test isequal(expr1.normexpr, NormExpr())
    expr2 = (trace(X)) - (x+1)
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, -x-1)
    @test isequal(expr2.matvars, {trace(X)})
    @test isequal(expr2.normexpr, NormExpr())
    @test_throws MethodError (trace(X)) * (x+1)
    @test_throws MethodError (trace(X)) / (x+1)
    con1 = (trace(X) >= (x+1))
    @test isequal(con1.terms, trace(X)-(x+1))
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (trace(X) <= (x+1))
    @test isequal(con2.terms, trace(X)-(x+1))
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (trace(X) == (x+1))
    @test isequal(con3.terms, trace(X)-(x+1))
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# MatrixFuncVar--DualExpr
# MatrixFuncVar--AbstractArray
let
    m = Model()
    @defSDPVar(m, X[4])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    @test_throws MethodError trace(X) + arr
    @test_throws MethodError trace(X) - arr
    expr3 = trace(X) * arr
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {trace(X)})
    @test isequal(expr3.coeffs, {arr})
    @test expr3.constant == spzeros(4,4)
    @test_throws MethodError trace(X) / arr
end
# MatrixFuncVar--MatrixVar
# MatrixFuncVar--SDPVar
# MatrixFuncVar--MatrixExpr
# MatrixFuncVar--MatrixFuncVar
# MatrixFuncVar--NormExpr
let 
    m = Model()
    @defSDPVar(m, X[2,2])
    @defMatrixVar(m, Y[3])
    mfv = trace(X)
    nexpr = norm(2Y-ones(3))
    expr1 = mfv + nexpr
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, AffExpr())
    @test isequal(expr1.matvars, {mfv})
    @test isequal(expr1.normexpr, nexpr)
    expr2 = mfv - nexpr
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, AffExpr())
    @test isequal(expr2.matvars, {mfv})
    @test isequal(expr2.normexpr, -nexpr)   
    @test_throws MethodError mfv * nexpr
    @test_throws MethodError mfv / nexpr
    con1 = (mfv >= nexpr)
    @test isequal(con1.terms, mfv-nexpr)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (mfv <= nexpr)
    @test isequal(con2.terms, mfv-nexpr)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (mfv == nexpr)
    @test isequal(con3.terms, mfv-nexpr)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# MatrixFuncVar--ScalarExpr
let 
    m = Model()
    @defVar(m,x)
    @defSDPVar(m, X[2,2])
    @defMatrixVar(m, Y[3])
    @defSDPVar(m, Z[2])
    mfv = trace(X)
    sexpr = norm(2Y-ones(3)) - 2trace(Z) + 3x - 1
    expr1 = mfv + sexpr
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, 3x-1)
    @test isequal(expr1.matvars, {mfv,-2trace(Z)})
    @test isequal(expr1.normexpr, norm(2Y-ones(3)))
    expr2 = mfv - sexpr
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, -3x+1)
    @test isequal(expr2.matvars, {mfv,2trace(Z)})
    @test isequal(expr2.normexpr, -norm(2Y-ones(3)))   
    @test_throws MethodError mfv * sexpr
    @test_throws MethodError mfv / sexpr
    con1 = (mfv >= sexpr)
    @test isequal(con1.terms, mfv-sexpr)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (mfv <= sexpr)
    @test isequal(con2.terms, mfv-sexpr)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (mfv == sexpr)
    @test isequal(con3.terms, mfv-sexpr)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end

# NormExpr--Number
let 
    m = Model()
    @defMatrixVar(m,X[2])
    expr = norm(2X - ones(2))
    expr1 = expr + 2
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, AffExpr(2.0))
    @test isempty(expr1.matvars)
    @test isequal(expr1.normexpr, expr)
    expr2 = expr - 2
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, AffExpr(-2.0))
    @test isempty(expr2.matvars)
    @test isequal(expr2.normexpr, expr)
    expr3 = expr * 2
    @test isa(expr3, NormExpr)
    @test isequal(expr3.vars, {2X-ones(2)})
    @test isequal(expr3.coeffs, [2.0])
    @test isequal(expr3.form, [:norm2])
    expr4 = expr / 2
    @test isa(expr4, NormExpr)
    @test isequal(expr4.vars, {2X-ones(2)})
    @test isequal(expr4.coeffs, [0.5])
    @test isequal(expr4.form, [:norm2])
    con1 = (expr >= 1)
    @test isequal(con1.terms, expr-1)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (expr <= 1)
    @test isequal(con2.terms, expr-1)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (expr == 1)
    @test isequal(con3.terms, expr-1)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# NormExpr--Variable
let
    m = Model()
    @defVar(m,x)
    @defMatrixVar(m,X[2])
    expr = norm(2X - ones(2))
    expr1 = expr + x
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, AffExpr(x))
    @test isempty(expr1.matvars)
    @test isequal(expr1.normexpr, expr)
    expr2 = expr - x
    @test isa(expr1, ScalarExpr)
    @test isequal(expr2.aff, -x)
    @test isempty(expr2.matvars)
    @test isequal(expr2.normexpr, expr)
    @test_throws MethodError expr * x
    @test_throws MethodError expr / x
    con1 = (expr >= x)
    @test isequal(con1.terms, expr-x)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (expr <= x)
    @test isequal(con2.terms, expr-x)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (expr == x)
    @test isequal(con3.terms, expr-x)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# NormExpr--AffExpr
let
    m = Model()
    @defVar(m,x)
    @defMatrixVar(m,Y[2])
    expr1 = norm(2Y-ones(2)) + (x+1)
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, x+1)
    @test isempty(expr1.matvars)
    @test isequal(expr1.normexpr, norm(2Y-ones(2)))
    expr2 = norm(2Y-ones(2)) - (x+1)
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, -x-1)
    @test isempty(expr2.matvars)
    @test isequal(expr2.normexpr, norm(2Y-ones(2)))
    @test_throws MethodError norm(2Y-ones(2)) * (x+1)
    @test_throws MethodError norm(2Y-ones(2)) / (x+1)
    con1 = (expr1 >= (x+1))
    @test isequal(con1.terms, expr1-(x+1))
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (expr1 <= (x+1))
    @test isequal(con2.terms, expr1-(x+1))
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (expr1 == (x+1))
    @test isequal(con3.terms, expr1-(x+1))
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# NormExpr--DualExpr
# NormExpr--AbstractArray
# NormExpr--MatrixVar
# NormExpr--SDPVar
# NormExpr--MatrixExpr
# NormExpr--MatrixFuncVar
# NormExpr--NormExpr
let
    m = Model()
    @defMatrixVar(m, X[2])
    @defMatrixVar(m, Y[3,3])
    n1 = 2norm(2X-ones(2))
    n2 = -vecnorm(-Y+2ones(3,3))
    expr1 = n1 + n2
    @test isa(expr1, NormExpr)
    @test isequal(expr1.vars, {2X-ones(2),-Y+2ones(3,3)})
    @test isequal(expr1.coeffs, {2,-1})
    @test expr1.form == [:norm2,:normfrob]
    expr2 = n1 - n2
    @test isa(expr2, NormExpr)
    @test isequal(expr2.vars, {2X-ones(2),-Y+2ones(3,3)})
    @test isequal(expr2.coeffs, {2,1})
    @test expr2.form == [:norm2,:normfrob]
    @test_throws MethodError n1 * n2
    @test_throws MethodError n1 / n2
    con1 = (n1 >= n2)
    @test isequal(con1.terms, n1-n2+0)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (n1 <= n2)
    @test isequal(con2.terms, n1-n2+0)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (n1 == n2)
    @test isequal(con3.terms, n1-n2+0)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# NormExpr--ScalarExpr
let
    m = Model()
    @defVar(m,x)
    @defMatrixVar(m, X[2])
    @defSDPVar(m, Y[3,3])
    n1 = 2norm(2X-ones(2))
    n2 = 2x + 1 - vecnorm(-Y+2ones(3,3)) + 2trace(Y)
    expr1 = n1 + n2
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, 2x+1)
    @test isequal(expr1.matvars, {2trace(Y)})
    @test isequal(expr1.normexpr, 2norm(2X-ones(2)) - vecnorm(-Y+2ones(3,3)))
    expr2 = n1 - n2
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, -2x-1)
    @test isequal(expr2.matvars, {-2trace(Y)})
    @test isequal(expr2.normexpr, 2norm(2X-ones(2)) + vecnorm(-Y+2ones(3,3)))
    @test_throws MethodError n1 * n2
    @test_throws MethodError n1 / n2
    con1 = (n1 >= n2)
    @test isequal(con1.terms, n1-n2)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (n1 <= n2)
    @test isequal(con2.terms, n1-n2)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (n1 == n2)
    @test isequal(con3.terms, n1-n2)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end

# ScalarExpr--Number
let 
    m = Model()
    @defVar(m,x)
    @defSDPVar(m,X[2])
    expr = 2x + 3 - 4trace(X) + norm(2X - ones(2,2))
    expr1 = expr + 2
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, 2x + 5)
    @test isequal(expr1.matvars, {-4trace(X)})
    @test isequal(expr1.normexpr, norm(2X - ones(2,2)))
    expr2 = expr - 2
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, 2x + 1)
    @test isequal(expr2.matvars, {-4trace(X)})
    @test isequal(expr2.normexpr, norm(2X - ones(2,2)))
    expr3 = expr * 2
    @test isa(expr3, ScalarExpr)
    @test isequal(expr3.aff, 4x + 6)
    @test isequal(expr3.matvars, {-8trace(X)})
    @test isequal(expr3.normexpr, 2norm(2X - ones(2,2)))
    expr4 = expr / 2
    @test isa(expr4, ScalarExpr)
    @test isequal(expr4.aff, x + 1.5)
    @test isequal(expr4.matvars, {-2trace(X)})
    @test isequal(expr4.normexpr, 0.5norm(2X - ones(2,2)))
    con1 = (expr >=  1)
    @test isequal(con1.terms, expr-1)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (expr <= 1)
    @test isequal(con2.terms, expr-1)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (expr == 1)
    @test isequal(con3.terms, expr-1)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# ScalarExpr--Variable
let
    m = Model()
    @defVar(m,x)
    @defVar(m,y)
    @defSDPVar(m,X[2])
    @defMatrixVar(m,Y[2])
    expr = 3 - x + 2trace(X) - 4norm(2Y+ones(2))
    expr1 = expr + y
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, 3 - x + y)
    @test isequal(expr1.matvars, {2trace(X)})
    @test isequal(expr1.normexpr, -4norm(2Y+ones(2)))
    expr2 = expr - y
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, 3 - x - y)
    @test isequal(expr2.matvars, {2trace(X)})
    @test isequal(expr2.normexpr, -4norm(2Y+ones(2)))
    @test_throws MethodError x * expr
    @test_throws MethodError x / expr
    con1 = (expr >=  x)
    @test isequal(con1.terms, expr-x)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (expr <= x)
    @test isequal(con2.terms, expr-x)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (expr == x)
    @test isequal(con3.terms, expr-x)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end

# ScalarExpr--AffExpr
let
    m = Model()
    @defVar(m,x)
    @defVar(m,y)
    @defSDPVar(m,X[2])
    @defMatrixVar(m,Y[2])
    expr1 = (-2 - y - 3trace(X) + 2norm(2Y-ones(2))) + (x+1)
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, -y+x-1)
    @test isequal(expr1.matvars, {-3trace(X)})
    @test isequal(expr1.normexpr, 2norm(2Y-ones(2)))
    expr2 = (-2 - y - 3trace(X) + 2norm(2Y-ones(2))) - (x+1)
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, -y-x-3)
    @test isequal(expr2.matvars, {-3trace(X)})
    @test isequal(expr2.normexpr, 2norm(2Y-ones(2)))
    @test_throws MethodError (-2 - y - 3trace(X) + 2norm(2Y-ones(2))) * (x+1)
    @test_throws MethodError (-2 - y - 3trace(X) + 2norm(2Y-ones(2))) / (x+1)
    con1 = (expr1 >=  (x+1))
    @test isequal(con1.terms, expr1-(x+1))
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (expr1 <= (x+1))
    @test isequal(con2.terms, expr1-(x+1))
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (expr1 == (x+1))
    @test isequal(con3.terms, expr1-(x+1))
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# ScalarExpr--DualExpr
# ScalarExpr--AbstractArray
let
    m = Model()
    @defVar(m,x)
    @defSDPVar(m, X[2])
    @defMatrixVar(m, Y[2])
    arr = 2eye(4,4) - diagm(ones(3), 1) - diagm(ones(3), -1)
    expr = 3 - x + 2trace(X)
    @test_throws MethodError expr + arr
    @test_throws MethodError expr - arr
    expr3 = expr * arr
    @test isa(expr3, DualExpr)
    @test isequal(expr3.vars, {2trace(X),x})
    @test length(expr3.coeffs) == 2
    @test expr3.coeffs[1] ==  arr
    @test expr3.coeffs[2] == -arr
    @test expr3.constant == 3arr
    nexpr = expr + vecnorm(Y)
    @test_throws Exception nexpr * arr
    @test_throws MethodError arr / trace(X)
end

# ScalarExpr--MatrixVar
# ScalarExpr--SDPVar
# ScalarExpr--MatrixExpr
# ScalarExpr--MatrixFuncVar
let 
    m = Model()
    @defVar(m,x)
    @defSDPVar(m, X[2,2])
    @defMatrixVar(m, Y[3])
    @defSDPVar(m, Z[2])
    mfv = trace(X)
    sexpr = norm(2Y-ones(3)) - 2trace(Z) + 3x - 1
    expr1 = sexpr + mfv
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, 3x-1)
    @test isequal(expr1.matvars, {-2trace(Z),mfv})
    @test isequal(expr1.normexpr, norm(2Y-ones(3)))
    expr2 = sexpr - mfv
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, 3x-1)
    @test isequal(expr2.matvars, {-2trace(Z),-mfv})
    @test isequal(expr2.normexpr, norm(2Y-ones(3)))   
    @test_throws MethodError sexpr * mfv
    @test_throws MethodError sexpr / mfv
    con1 = (mfv >=  sexpr)
    @test isequal(con1.terms, mfv-sexpr)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (mfv <= sexpr)
    @test isequal(con2.terms, mfv-sexpr)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (mfv == sexpr)
    @test isequal(con3.terms, mfv-sexpr)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end

# ScalarExpr--NormExpr
let
    m = Model()
    @defVar(m,x)
    @defMatrixVar(m, X[2])
    @defSDPVar(m, Y[3,3])
    n1 = 2x + 1 - vecnorm(-Y+2ones(3,3)) + 2trace(Y)
    n2 = 2norm(2X-ones(2))
    expr1 = n1 + n2
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, 2x+1)
    @test isequal(expr1.matvars, {2trace(Y)})
    @test isequal(expr1.normexpr, -vecnorm(-Y+2ones(3,3)) + 2norm(2X-ones(2)))
    expr2 = n1 - n2
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, 2x+1)
    @test isequal(expr2.matvars, {2trace(Y)})
    @test isequal(expr2.normexpr, - vecnorm(-Y+2ones(3,3)) - 2norm(2X-ones(2)))
    @test_throws MethodError n1 * n2
    @test_throws MethodError n1 / n2
    con1 = (n1 >=  n2)
    @test isequal(con1.terms, n1-n2)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (n1 <= n2)
    @test isequal(con2.terms, n1-n2)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (n1 == n2)
    @test isequal(con3.terms, n1-n2)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
# ScalarExpr--ScalarExpr
let
    m = Model()
    @defVar(m,x)
    @defMatrixVar(m,X[2,2])
    @defSDPVar(m,Y[2,2])
    s1 = 2x - 1 + vecnorm(2X) - 2trace(Y)
    s2 = trace(Y) - 3 - x - vecnorm(X)
    expr1 = s1 + s2
    @test isa(expr1, ScalarExpr)
    @test isequal(expr1.aff, 2x - x - 4)
    @test isequal(expr1.matvars, {-2trace(Y),trace(Y)})
    @test isequal(expr1.normexpr, vecnorm(2X) - vecnorm(X))
    expr2 = s1 - s2
    @test isa(expr2, ScalarExpr)
    @test isequal(expr2.aff, 2x + x + 2)
    @test isequal(expr2.matvars, {-2trace(Y),-trace(Y)})
    @test isequal(expr2.normexpr, vecnorm(2X) + vecnorm(X))
    @test_throws MethodError s1 * s2
    @test_throws MethodError s1 / s2
    con1 = (s1 >=  s2)
    @test isequal(con1.terms, s1-s2)
    @test con1.lb == 0.0
    @test con1.ub == Inf
    con2 = (s1 <= s2)
    @test isequal(con2.terms, s1-s2)
    @test con2.lb == -Inf
    @test con2.ub ==  0.0
    con3 = (s1 == s2)
    @test isequal(con3.terms, s1-s2)
    @test con3.lb == 0.0
    @test con3.ub == 0.0
end
