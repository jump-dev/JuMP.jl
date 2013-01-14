param numFacility default 100;
param numCustomer default 100;
param numLocation default 1000;

param custLoc {a in 1..numCustomer} := 1 + floor(1000 * Uniform01());

var s {i in 1..numLocation}, >= 0, <= 1;
var x {i in 1..numLocation, a in 1..numCustomer}, >= 0, <= 1;

maximize OBJ:
  sum {i in 1..numLocation, a in 1..numCustomer} abs(custLoc[a] - i)*x[i,a];

subject to XSUM {a in 1..numCustomer}:
  sum {i in 1..numLocation} x[i,a] == 1;

subject to XSLINK {a in 1..numCustomer, i in 1..numLocation}:
  x[i,a] <= s[i];

subject to SSUM:
  sum {i in 1..numLocation} s[i] == numFacility;

write mpmedian;

