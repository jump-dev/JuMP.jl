using JuMP

m = Model()
@defSDPVar(m, X[3] >= zeros(3,3))
@defSDPVar(m, Y[2] >= zeros(2,2))

C = eye(3,3)
A1 = zeros(3,3)
A1[1,1] = 1.0
A2 = zeros(3,3)
A2[2,2] = 1.0
A3 = zeros(3,3)
A3[3,3] = 1.0
D = eye(2,2)
B1 = ones(2,2)
B2 = zeros(2,2)
B2[1,1] = 1
B3 = zeros(2,2)
B3[1,2] = B3[2,1] = 2

setObjective(m, :Min, trace(C*X)+1+trace(D*Y))
addConstraint(m, trace(A1*X-eye(3,3)/3) == 0)
addConstraint(m, 2*trace(A2*X) == 1)
addConstraint(m, trace(A3*X) >= 2)
addConstraint(m, trace(B1*Y) == 1)
addConstraint(m, trace(B2*Y) == 0)
addConstraint(m, trace(B3*Y) <= 0)
addConstraint(m, trace(A1*X)+trace(B1*Y) >= 1)
addConstraint(m, Y[2,2] == 1)

solve(m)

println("objective value = $(getObjectiveValue(m))")
println("X = $(getValue(X))")
println("Y = $(getValue(Y))")
