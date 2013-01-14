printf "Time at start: %f\n",time();
param N default 10000;
param M default 1000;

param A {j in 1..N, i in 1..M} := Uniform(1,10);
param B {i in 1..M} := Uniform(N,N*10);
param C {j in 1..N} := Uniform(1,10);

var x {j in 1..N} >= 0, <= 1;

maximize OBJ:
  sum {j in 1..N} C[j] * x[j];

subject to CON {i in 1..M}:
  sum {j in 1..N} A[j,i] * x[j] <= B[i];

write mtest;
printf "Time at end: %f\n",time();
