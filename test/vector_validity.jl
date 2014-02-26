using JuMP
using Base.Test


m = Model()
@defVar(m,x[1:3])
c = [1:3]
@test affToStr(dot(c,x)) == "1.0 x[1] + 2.0 x[2] + 3.0 x[3]"
@test affToStr(dot(x,c)) == "1.0 x[1] + 2.0 x[2] + 3.0 x[3]"

A = [1 2 ; 3 4]
@defVar(m,y[1:2,1:2])
@test affToStr(bigdot(A,y)) == "1.0 y[1,1] + 2.0 y[1,2] + 3.0 y[2,1] + 4.0 y[2,2]"
@test affToStr(bigdot(y,A)) == "1.0 y[1,1] + 2.0 y[1,2] + 3.0 y[2,1] + 4.0 y[2,2]"

B = ones(2,2,2)
@defVar(m,z[1:2,1:2,1:2])
@test affToStr(bigdot(B,z)) == "1.0 z[1,1,1] + 1.0 z[1,1,2] + 1.0 z[1,2,1] + 1.0 z[1,2,2] + 1.0 z[2,1,1] + 1.0 z[2,1,2] + 1.0 z[2,2,1] + 1.0 z[2,2,2]"
@test affToStr(bigdot(z,B)) == "1.0 z[1,1,1] + 1.0 z[1,1,2] + 1.0 z[1,2,1] + 1.0 z[1,2,2] + 1.0 z[2,1,1] + 1.0 z[2,1,2] + 1.0 z[2,2,1] + 1.0 z[2,2,2]"
