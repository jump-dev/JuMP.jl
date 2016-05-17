#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/print.jl
# Testing $fa pretty-printing-related functionality
#############################################################################
using JuMP, FactCheck
import JuMP.REPLMode, JuMP.IJuliaMode
import JuMP.repl, JuMP.ijulia

# Helper function to test IO methods work correctly
function io_test(mode, obj, exp_str; repl=:both)
    if mode == REPLMode
        repl != :show  && @fact sprint(print, obj) --> exp_str
        repl != :print && @fact sprint(show,  obj) --> exp_str
    else
        @fact sprint(writemime, "text/latex", obj) --> "\$\$ $exp_str \$\$"
    end
end


facts("[print] JuMPContainer{Variable}") do
    le, ge, fa = repl[:leq], repl[:geq], repl[:for_all]
    inset, dots = repl[:in], repl[:dots]
    infty, union = repl[:infty], repl[:union]

    m = Model()

    #------------------------------------------------------------------
    # Test bound printing
    context("bound printing") do
    @variable(m,      bnd_free[2:5])
    @variable(m,      bnd_lowb[2:5] >= 2)
    @variable(m,      bnd_high[2:5] <= 5)
    @variable(m, 2 <= bnd_both[2:5] <= 5)
    @variable(m,      bnd_difflo[i=2:5] >= i)
    @variable(m,      bnd_diffup[i=2:5] <= i)
    @variable(m, i <= bnd_diffbo[i=2:5] <= 2i)
    @variable(m, i <= bnd_difflo_with_up[i=2:5] <= 5)
    @variable(m, 2 <= bnd_diffup_with_lo[i=2:5] <= i)

    io_test(REPLMode, bnd_free, "bnd_free[i] free $fa i $inset {2,3,4,5}")
    io_test(REPLMode, bnd_lowb, "bnd_lowb[i] $ge 2 $fa i $inset {2,3,4,5}")
    io_test(REPLMode, bnd_high, "bnd_high[i] $le 5 $fa i $inset {2,3,4,5}")
    io_test(REPLMode, bnd_both, "2 $le bnd_both[i] $le 5 $fa i $inset {2,3,4,5}")
    io_test(REPLMode, bnd_difflo, "bnd_difflo[i] $ge $dots $fa i $inset {2,3,4,5}")
    io_test(REPLMode, bnd_diffup, "bnd_diffup[i] $le $dots $fa i $inset {2,3,4,5}")
    io_test(REPLMode, bnd_diffbo, "$dots $le bnd_diffbo[i] $le $dots $fa i $inset {2,3,4,5}")
    io_test(REPLMode, bnd_difflo_with_up, "$dots $le bnd_difflo_with_up[i] $le 5 $fa i $inset {2,3,4,5}")
    io_test(REPLMode, bnd_diffup_with_lo, "2 $le bnd_diffup_with_lo[i] $le $dots $fa i $inset {2,3,4,5}")

    io_test(IJuliaMode, bnd_free, "bnd_free_{i} free \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_lowb, "bnd_lowb_{i} \\geq 2 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_high, "bnd_high_{i} \\leq 5 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_both, "2 \\leq bnd_both_{i} \\leq 5 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_difflo, "bnd_difflo_{i} \\geq \\dots \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_diffup, "bnd_diffup_{i} \\leq \\dots \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_diffbo, "\\dots \\leq bnd_diffbo_{i} \\leq \\dots \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_difflo_with_up, "\\dots \\leq bnd_difflo_with_up_{i} \\leq 5 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_diffup_with_lo, "2 \\leq bnd_diffup_with_lo_{i} \\leq \\dots \\quad\\forall i \\in \\{2,3,4,5\\}")
    end

    #------------------------------------------------------------------
    # Test index set printing
    context("index set printing") do
    @variable(m, rng_unit1[1:10])  # Array{Variable}
    @variable(m, rng_unit2[-2:3])  # JuMPArray
    @variable(m, rng_unit3[[1:10;]])  # JuMPDict
    @variable(m, rng_step1[1:2:10])
    @variable(m, rng_step2[-2:5:10])
    @variable(m, rng_step3[1:5:3])
    @variable(m, rng_step4[0:2:2])
    @variable(m, arr_1[[:a,:b,:c]])
    @variable(m, arr_2[[:a,1,"test"]])
    @variable(m, arr_3[[:apple,:banana,:carrot,:diamonds]])
    @variable(m, rng2_1[1:10,[:a,:b,:c]])
    @variable(m, tri_1[i=1:3,j=i:3])
    @variable(m, tri_2[i=1:3,j=-i])
    @variable(m, tri_3[(i,j)=[(i,i+2) for i in 1:5],k=i:j])

    io_test(REPLMode, rng_unit1, "rng_unit1[i] free $fa i $inset {1,2,$dots,9,10}")
    io_test(REPLMode, rng_unit2, "rng_unit2[i] free $fa i $inset {-2,-1,$dots,2,3}")
    io_test(REPLMode, rng_unit3, "rng_unit3[i] free $fa i $inset {1,2,$dots,9,10}")
    io_test(REPLMode, rng_step1, "rng_step1[i] free $fa i $inset {1,3,$dots,7,9}")
    io_test(REPLMode, rng_step2, "rng_step2[i] free $fa i $inset {-2,3,8}")
    io_test(REPLMode, rng_step3, "rng_step3[i] free $fa i $inset {1}")
    io_test(REPLMode, rng_step4, "rng_step4[i] free $fa i $inset {0,2}")
    io_test(REPLMode, arr_1, "arr_1[i] free $fa i $inset {a,b,c}")
    io_test(REPLMode, arr_2, "arr_2[i] free $fa i $inset {a,1,test}")
    io_test(REPLMode, arr_3, "arr_3[i] free $fa i $inset {apple,banana,carrot,diamonds}")
    io_test(REPLMode, rng2_1, "rng2_1[i,j] free $fa i $inset {1,2,$dots,9,10}, j $inset {a,b,c}")
    io_test(REPLMode, tri_1, "tri_1[i,j] free $fa i $inset {1,2,3}, j $inset {$dots}")
    io_test(REPLMode, tri_2, "tri_2[i,j] free $fa i $inset {1,2,3}, j $inset {$dots}")
    io_test(REPLMode, tri_3, "tri_3[(i,j),k] free $fa (i,j) $inset {(1,3),(2,4),$dots,(4,6),(5,7)}, k $inset {$dots}")

    io_test(IJuliaMode, rng_unit1, "rng_unit1_{i} free \\quad\\forall i \\in \\{1,2,\\dots,9,10\\}")
    io_test(IJuliaMode, rng_unit2, "rng_unit2_{i} free \\quad\\forall i \\in \\{-2,-1,\\dots,2,3\\}")
    io_test(IJuliaMode, rng_unit3, "rng_unit3_{i} free \\quad\\forall i \\in \\{1,2,\\dots,9,10\\}")
    io_test(IJuliaMode, rng_step1, "rng_step1_{i} free \\quad\\forall i \\in \\{1,3,\\dots,7,9\\}")
    io_test(IJuliaMode, rng_step2, "rng_step2_{i} free \\quad\\forall i \\in \\{-2,3,8\\}")
    io_test(IJuliaMode, rng_step3, "rng_step3_{i} free \\quad\\forall i \\in \\{1\\}")
    io_test(IJuliaMode, rng_step4, "rng_step4_{i} free \\quad\\forall i \\in \\{0,2\\}")
    io_test(IJuliaMode, arr_1, "arr_1_{i} free \\quad\\forall i \\in \\{a,b,c\\}")
    io_test(IJuliaMode, arr_2, "arr_2_{i} free \\quad\\forall i \\in \\{a,1,test\\}")
    io_test(IJuliaMode, arr_3, "arr_3_{i} free \\quad\\forall i \\in \\{apple,banana,carrot,diamonds\\}")
    io_test(IJuliaMode, rng2_1, "rng2_1_{i,j} free \\quad\\forall i \\in \\{1,2,\\dots,9,10\\}, j \\in \\{a,b,c\\}")
    io_test(IJuliaMode, tri_1, "tri_1_{i,j} free \\quad\\forall i \\in \\{1,2,3\\}, j \\in \\{\\dots\\}")
    io_test(IJuliaMode, tri_2, "tri_2_{i,j} free \\quad\\forall i \\in \\{1,2,3\\}, j \\in \\{\\dots\\}")
    io_test(IJuliaMode, tri_3, "tri_3_{(i,j),k} free \\quad\\forall (i,j) \\in \\{(1,3),(2,4),\\dots,(4,6),(5,7)\\}, k \\in \\{\\dots\\}")
    end

    #------------------------------------------------------------------
    # Test category printing
    context("category printing") do
    @variable(m, cat_bin[1:3], Bin)
    @variable(m, 2 <= cat_int[1:3] <= 5, Int)
    @variable(m, cat_semiint_both[2:3] >= 2, SemiInt)
    @variable(m, i <= cat_semiint_difflow[i=2:3] <= 4, SemiInt)
    @variable(m, 2 <= cat_semiint_diffup[i=2:3] <= i, SemiInt)
    @variable(m, i <= cat_semiint_none[i=2:3] <= 2i, SemiInt)
    @variable(m, cat_semicont_both[2:3] >= 2, SemiCont)
    @variable(m, i <= cat_semicont_difflow[i=2:3] <= 4, SemiCont)
    @variable(m, 2 <= cat_semicont_diffup[i=2:3] <= i, SemiCont)
    @variable(m, i <= cat_semicont_none[i=2:3] <= 2i, SemiCont)
    @variable(m, fixed_var[i=2:3] == i)

    io_test(REPLMode, cat_bin, "cat_bin[i] $inset {0,1} $fa i $inset {1,2,3}")
    io_test(REPLMode, cat_int, "2 $le cat_int[i] $le 5, integer, $fa i $inset {1,2,3}")
    io_test(REPLMode, cat_semiint_both, "cat_semiint_both[i] $inset {2,$dots,$infty} $union {0} $fa i $inset {2,3}")
    io_test(REPLMode, cat_semiint_difflow, "cat_semiint_difflow[i] $inset {$dots,$dots,4} $union {0} $fa i $inset {2,3}")
    io_test(REPLMode, cat_semiint_diffup, "cat_semiint_diffup[i] $inset {2,$dots,$dots} $union {0} $fa i $inset {2,3}")
    io_test(REPLMode, cat_semiint_none, "cat_semiint_none[i] $inset {$dots,$dots,$dots} $union {0} $fa i $inset {2,3}")
    io_test(REPLMode, cat_semicont_both, "cat_semicont_both[i] $inset [2,$infty] $union {0} $fa i $inset {2,3}")
    io_test(REPLMode, cat_semicont_difflow, "cat_semicont_difflow[i] $inset [$dots,4] $union {0} $fa i $inset {2,3}")
    io_test(REPLMode, cat_semicont_diffup, "cat_semicont_diffup[i] $inset [2,$dots] $union {0} $fa i $inset {2,3}")
    io_test(REPLMode, cat_semicont_none, "cat_semicont_none[i] $inset [$dots,$dots] $union {0} $fa i $inset {2,3}")
    io_test(REPLMode, fixed_var, "fixed_var[i] = $dots $fa i $inset {2,3}")

    io_test(IJuliaMode, cat_bin, "cat_bin_{i} \\in \\{0,1\\} \\quad\\forall i \\in \\{1,2,3\\}")
    io_test(IJuliaMode, cat_int, "2 \\leq cat_int_{i} \\leq 5, \\in \\mathbb{Z}, \\quad\\forall i \\in \\{1,2,3\\}")
    io_test(IJuliaMode, cat_semiint_both, "cat_semiint_both_{i} \\in \\{2,\\dots,\\infty\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semiint_difflow, "cat_semiint_difflow_{i} \\in \\{\\dots,\\dots,4\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semiint_diffup, "cat_semiint_diffup_{i} \\in \\{2,\\dots,\\dots\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semiint_none, "cat_semiint_none_{i} \\in \\{\\dots,\\dots,\\dots\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_both, "cat_semicont_both_{i} \\in \\[2,\\infty\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_difflow, "cat_semicont_difflow_{i} \\in \\[\\dots,4\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_diffup, "cat_semicont_diffup_{i} \\in \\[2,\\dots\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_none, "cat_semicont_none_{i} \\in \\[\\dots,\\dots\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, fixed_var, "fixed_var_{i} = \\dots \\quad\\forall i \\in \\{2,3\\}")
    end

    #------------------------------------------------------------------
    # Tests for particular issues
    context("Empty JuMPContainer printing (#124)") do
    @variable(m, empty_free[1:0])
    io_test(REPLMode, empty_free, "Empty Array{Variable} (no indices)")
    io_test(IJuliaMode, empty_free, "Empty Array{Variable} (no indices)")
    @variable(m, empty_set[[]])
    io_test(REPLMode, empty_set, "empty_set (no indices)")
    io_test(IJuliaMode, empty_set, "empty_set (no indices)")
    end
end



facts("[print] JuMPContainer{Number}") do
    # The same output for REPL and IJulia, so only testing one
    mod = Model()
    @variable(mod, i*j <= w[i=9:10, [:Apple,5,:Banana], j=-1:+1] <= i*j)
    @variable(mod, i*j*k <= x[i=9:11,j=99:101,k=3:4] <= i*j*k)
    @variable(mod, i*j <= y[i=9:11,j=i:11] <= i*j)
    @variable(mod, j <= z[i=[:a,'b'],j=1:3] <= j)
    solve(mod)

    # Deal with hashing variations
    if hash(5) < hash(:Apple)
        io_test(REPLMode, getvalue(w), """
w: 3 dimensions:
[ 9,:,:]
  [ 9, Apple,:]
    [ 9, Apple,-1] = -9.0
    [ 9, Apple, 0] = 0.0
    [ 9, Apple, 1] = 9.0
  [ 9,     5,:]
    [ 9,     5,-1] = -9.0
    [ 9,     5, 0] = 0.0
    [ 9,     5, 1] = 9.0
  [ 9,Banana,:]
    [ 9,Banana,-1] = -9.0
    [ 9,Banana, 0] = 0.0
    [ 9,Banana, 1] = 9.0
[10,:,:]
  [10, Apple,:]
    [10, Apple,-1] = -10.0
    [10, Apple, 0] = 0.0
    [10, Apple, 1] = 10.0
  [10,     5,:]
    [10,     5,-1] = -10.0
    [10,     5, 0] = 0.0
    [10,     5, 1] = 10.0
  [10,Banana,:]
    [10,Banana,-1] = -10.0
    [10,Banana, 0] = 0.0
    [10,Banana, 1] = 10.0
""", repl=:print)
    else
        io_test(REPLMode, getvalue(w), """
w: 3 dimensions:
[ 9,:,:]
  [ 9, Apple,:]
    [ 9, Apple,-1] = -9.0
    [ 9, Apple, 0] = 0.0
    [ 9, Apple, 1] = 9.0
  [ 9,     5,:]
    [ 9,     5,-1] = -9.0
    [ 9,     5, 0] = 0.0
    [ 9,     5, 1] = 9.0
  [ 9,Banana,:]
    [ 9,Banana,-1] = -9.0
    [ 9,Banana, 0] = 0.0
    [ 9,Banana, 1] = 9.0
[10,:,:]
  [10, Apple,:]
    [10, Apple,-1] = -10.0
    [10, Apple, 0] = 0.0
    [10, Apple, 1] = 10.0
  [10,     5,:]
    [10,     5,-1] = -10.0
    [10,     5, 0] = 0.0
    [10,     5, 1] = 10.0
  [10,Banana,:]
    [10,Banana,-1] = -10.0
    [10,Banana, 0] = 0.0
    [10,Banana, 1] = 10.0
""", repl=:print)
    end

    io_test(REPLMode, getvalue(x), """
x: 3 dimensions:
[ 9,:,:]
  [ 9, 99,:]
    [ 9, 99,3] = 2673.0
    [ 9, 99,4] = 3564.0
  [ 9,100,:]
    [ 9,100,3] = 2700.0
    [ 9,100,4] = 3600.0
  [ 9,101,:]
    [ 9,101,3] = 2727.0
    [ 9,101,4] = 3636.0
[10,:,:]
  [10, 99,:]
    [10, 99,3] = 2970.0
    [10, 99,4] = 3960.0
  [10,100,:]
    [10,100,3] = 3000.0
    [10,100,4] = 4000.0
  [10,101,:]
    [10,101,3] = 3030.0
    [10,101,4] = 4040.0
[11,:,:]
  [11, 99,:]
    [11, 99,3] = 3267.0
    [11, 99,4] = 4356.0
  [11,100,:]
    [11,100,3] = 3300.0
    [11,100,4] = 4400.0
  [11,101,:]
    [11,101,3] = 3333.0
    [11,101,4] = 4444.0
""", repl=:print)

    io_test(REPLMode, getvalue(y), """
y: 2 dimensions, 6 entries:
 [ 9, 9] = 81.0
 [ 9,10] = 90.0
 [ 9,11] = 99.0
 [10,10] = 100.0
 [10,11] = 110.0
 [11,11] = 121.0""")

    io_test(REPLMode, getvalue(z), """
z: 2 dimensions:
[a,:]
  [a,1] = 1.0
  [a,2] = 2.0
  [a,3] = 3.0
[b,:]
  [b,1] = 1.0
  [b,2] = 2.0
  [b,3] = 3.0
""")

end



facts("[print] SOS constraints") do
    modS = Model()
    a = [1,2,3]
    @variable(modS, x[1:3], Bin)
    addSOS1(modS, [a[i]x[i] for i in 1:3])
    s1 = JuMP.SOSConstraint([x[i] for i in 1:3],
                            [a[i] for i in 1:3], :SOS1)
    io_test(REPLMode, s1, "SOS1: {1 x[1], 2 x[2], 3 x[3]}")
    io_test(IJuliaMode, s1, "SOS1: \\{1 x[1], 2 x[2], 3 x[3]\\}")

    b = [5,4,7,2,1]
    @variable(modS, y[1:5], Bin)
    s2 = JuMP.SOSConstraint([y[i] for i in 1:5],
                            [b[i] for i in 1:5], :SOS2)
    io_test(REPLMode, s2, "SOS2: {5 y[1], 4 y[2], 7 y[3], 2 y[4], 1 y[5]}")
    io_test(IJuliaMode, s2, "SOS2: \\{5 y[1], 4 y[2], 7 y[3], 2 y[4], 1 y[5]\\}")
end



facts("[print] Model") do
    le, ge, fa = repl[:leq], repl[:geq], repl[:for_all]
    inset, dots = repl[:in], repl[:dots]
    infty, union = repl[:infty], repl[:union]
    Vert, sub2 = repl[:Vert], repl[:sub2]

    #------------------------------------------------------------------

    mod_1 = Model()
    @variable(mod_1, a>=1)
    @variable(mod_1, b<=1)
    @variable(mod_1, -1<=c<=1)
    @variable(mod_1, a1>=1,Int)
    @variable(mod_1, b1<=1,Int)
    @variable(mod_1, -1<=c1<=1,Int)
    @variable(mod_1, x, Bin)
    @variable(mod_1, y)
    @variable(mod_1, z, Int)
    @variable(mod_1, sos[1:3], Bin)
    @variable(mod_1, 2 <= si <= 3, SemiInt)
    @variable(mod_1, 2 <= sc <= 3, SemiCont)
    @variable(mod_1, fi == 9)
    @objective(mod_1, Max, a - b + 2a1 - 10x)
    @constraint(mod_1, a + b - 10c - 2x + c1 <= 1)
    @constraint(mod_1, a*b <= 2)
    addSOS1(mod_1, [i*sos[i] for i in 1:3])
    @constraint(mod_1, norm(sos) + a <= 1)

    io_test(REPLMode, mod_1, """
Max a - b + 2 a1 - 10 x
Subject to
 a + b - 10 c - 2 x + c1 $le 1
 a*b - 2 $le 0
 SOS1: {1 sos[1], 2 sos[2], 3 sos[3]}
 $(Vert)[sos[1],sos[2],sos[3]]$(Vert)$(sub2) $le -a + 1
 sos[i] $inset {0,1} $fa i $inset {1,2,3}
 a $ge 1
 b $le 1
 -1 $le c $le 1
 a1 $ge 1, integer
 b1 $le 1, integer
 -1 $le c1 $le 1, integer
 x $inset {0,1}
 y free
 z free, integer
 si $inset {2,$dots,3} $union {0}
 sc $inset [2,3] $union {0}
 fi = 9
""", repl=:print)

    io_test(REPLMode, mod_1, """
Maximization problem with:
 * 1 linear constraint
 * 1 quadratic constraint
 * 1 SOS constraint
 * 1 SOC constraint
 * 15 variables: 4 binary, 4 integer, 1 semicontinuous, 1 semi-integer
Solver is default solver""", repl=:show)

    io_test(IJuliaMode, mod_1, """
\\begin{alignat*}{1}\\max\\quad & a - b + 2 a1 - 10 x\\\\
\\text{Subject to} \\quad & a + b - 10 c - 2 x + c1 \\leq 1\\\\
 & a\\times b - 2 \\leq 0\\\\
 & SOS1: \\{1 sos[1], 2 sos[2], 3 sos[3]\\}\\\\
 & \\Vert[sos_{1},sos_{2},sos_{3}]\\Vert_2 $le -a + 1\\\\
 & sos_{i} \\in \\{0,1\\} \\quad\\forall i \\in \\{1,2,3\\}\\\\
 & a \\geq 1\\\\
 & b \\leq 1\\\\
 & -1 \\leq c \\leq 1\\\\
 & a1 \\geq 1, \\in \\mathbb{Z}\\\\
 & b1 \\leq 1, \\in \\mathbb{Z}\\\\
 & -1 \\leq c1 \\leq 1, \\in \\mathbb{Z}\\\\
 & x \\in \\{0,1\\}\\\\
 & y free\\\\
 & z free, \\in \\mathbb{Z}\\\\
 & si \\in \\{2,\\dots,3\\} \\cup \\{0\\}\\\\
 & sc \\in \\[2,3\\] \\cup \\{0\\}\\\\
 & fi = 9\\\\
\\end{alignat*}
""")

    #------------------------------------------------------------------

    mod_2 = Model()
    @variable(mod_2, x, Bin)
    @variable(mod_2, y, Int)
    @constraint(mod_2, x*y <= 1)

    io_test(REPLMode, mod_2, """
Feasibility problem with:
 * 0 linear constraints
 * 1 quadratic constraint
 * 2 variables: 1 binary, 1 integer
Solver is default solver""", repl=:show)

    mod_2 = Model()
    @variable(mod_2, x)
    @constraint(mod_2, x <= 3)

    io_test(REPLMode, mod_2, """
Feasibility problem with:
 * 1 linear constraint
 * 1 variable
Solver is default solver""", repl=:show)

    #------------------------------------------------------------------

    mod_3 = Model()

    @variable(mod_3, x[1:5])
    @NLconstraint(mod_3, x[1]*x[2] == 1)
    @NLconstraint(mod_3, x[3]*x[4] == 1)
    @NLconstraint(mod_3, x[5]*x[1] == 1)
    @NLobjective(mod_3, Min, x[1]*x[3])

    io_test(REPLMode, mod_3, """
Min (nonlinear expression)
Subject to
 3 nonlinear constraints
 x[i] free $fa i $inset {1,2,$dots,4,5}
""", repl=:print)
    io_test(REPLMode, mod_3, """
Minimization problem with:
 * 0 linear constraints
 * 3 nonlinear constraints
 * 5 variables
Solver is default solver""", repl=:show)
    io_test(IJuliaMode, mod_3, """
\\begin{alignat*}{1}\\min\\quad & (nonlinear expression)\\\\
\\text{Subject to} \\quad & 3 nonlinear constraints\\\\
 & x_{i} free \\quad\\forall i \\in \\{1,2,\\dots,4,5\\}\\\\
\\end{alignat*}
""", repl=:print)
end

facts("[print] changing variable categories") do
    le, ge, fa = repl[:leq], repl[:geq], repl[:for_all]
    inset, dots = repl[:in], repl[:dots]
    infty, union = repl[:infty], repl[:union]

    mod = Model()
    @variable(mod, x[1:3])
    @variable(mod, y[i=1:3,i:3])
    setcategory(x[3], :SemiCont)
    setcategory(y[1,3], :Int)

    io_test(REPLMode, mod, """
Min 0
Subject to
 x[1] free
 x[2] free
 x[3] $inset [-$infty,$infty] $union {0}
 y[1,1] free
 y[1,2] free
 y[1,3] free, integer
 y[2,2] free
 y[2,3] free
 y[3,3] free
""", repl=:print)

    io_test(IJuliaMode, mod, """
\\begin{alignat*}{1}\\min\\quad & 0\\\\
\\text{Subject to} \\quad & x_{1} free\\\\
 & x_{2} free\\\\
 & x_{3} \\in \\[-\\infty,\\infty\\] \\cup \\{0\\}\\\\
 & y_{1,1} free\\\\
 & y_{1,2} free\\\\
 & y_{1,3} free, \\in \\mathbb{Z}\\\\
 & y_{2,2} free\\\\
 & y_{2,3} free\\\\
 & y_{3,3} free\\\\
\\end{alignat*}
""")
end

facts("[print] expressions") do
    # Most of the expression logic is well covered by test/operator.jl
    # This is really just to check IJulia printing for expressions
    le, ge = repl[:leq], repl[:geq]

    #------------------------------------------------------------------
    mod = Model()
    @variable(mod, x[1:5])
    @variable(mod, y[i=2:4,j=i:5])
    @variable(mod, z)

    @constraint(mod, x[1] + 2*y[2,3] <= 3)
    io_test(REPLMode, mod.linconstr[end], "x[1] + 2 y[2,3] $le 3")
    io_test(IJuliaMode, mod.linconstr[end], "x_{1} + 2 y_{2,3} \\leq 3")

    @constraint(mod, (x[1]+x[2])*(y[2,2]+3.0) <= 1)
    io_test(REPLMode, mod.quadconstr[end], "x[1]*y[2,2] + x[2]*y[2,2] + 3 x[1] + 3 x[2] - 1 $le 0")
    io_test(IJuliaMode, mod.quadconstr[end], "x_{1}\\times y_{2,2} + x_{2}\\times y_{2,2} + 3 x_{1} + 3 x_{2} - 1 \\leq 0")

    @constraint(mod, (y[2,2]+3.0)*(x[1]+x[2]) <= 1)
    io_test(REPLMode, mod.quadconstr[end], "x[1]*y[2,2] + x[2]*y[2,2] + 3 x[1] + 3 x[2] - 1 $le 0")
    io_test(IJuliaMode, mod.quadconstr[end], "x_{1}\\times y_{2,2} + x_{2}\\times y_{2,2} + 3 x_{1} + 3 x_{2} - 1 \\leq 0")
end



facts("[print] Variable") do
    m = Model()
    @variable(m, 0 <= x <= 2, inconstraints=ConstraintRef{Model,LinearConstraint}[], objective=0.0, coefficients=Float64[] )

    @fact    getname(x) --> "x"
    io_test(REPLMode,   x, "x")
    io_test(IJuliaMode, x, "x")

    setname(x, "x2")
    @fact    getname(x) --> "x2"
    io_test(REPLMode,   x, "x2")
    io_test(IJuliaMode, x, "x2")

    setname(x, "")
    @fact    getname(x) --> "col_1"
    io_test(REPLMode,   x, "col_1")
    io_test(IJuliaMode, x, "col_1")

    @variable(m, z[1:2,3:5])
    @fact       getname(z[1,3]) --> "z[1,3]"
    io_test(REPLMode,   z[1,3],    "z[1,3]")
    io_test(IJuliaMode, z[1,3],    "z_{1,3}")
    @fact       getname(z[2,4]) --> "z[2,4]"
    io_test(REPLMode,   z[2,4],    "z[2,4]")
    io_test(IJuliaMode, z[2,4],    "z_{2,4}")
    @fact       getname(z[2,5]) --> "z[2,5]"
    io_test(REPLMode,   z[2,5],    "z[2,5]")
    io_test(IJuliaMode, z[2,5],    "z_{2,5}")

    @variable(m, w[3:9,["red","blue","green"]])
    @fact    getname(w[7,"green"]) --> "w[7,green]"
    io_test(REPLMode,   w[7,"green"], "w[7,green]")
    io_test(IJuliaMode, w[7,"green"], "w_{7,green}")

    rng = 2:5
    @variable(m, v[rng,rng,rng,rng,rng,rng,rng])
    a_v = v[4,5,2,3,2,2,4]
    @fact    getname(a_v) --> "v[4,5,2,3,2,2,4]"
    io_test(REPLMode,   a_v, "v[4,5,2,3,2,2,4]")
    io_test(IJuliaMode, a_v, "v_{4,5,2,3,2,2,4}")
end

facts("[print] User-created Array{Variable}") do
    m = Model()
    @variable(m, x)
    @variable(m, y)

    v = [x,y,x]
    A = [x y; y x]
    io_test(REPLMode,   v, "[x,y,x]")
    io_test(IJuliaMode, v, "[x,y,x]")

    if VERSION >= v"0.5.0-dev+3642"
        io_test(REPLMode,   A, """
2×2 Array{JuMP.Variable,2}:
 x  y
 y  x""")
        io_test(IJuliaMode, A, """
2×2 Array{JuMP.Variable,2}:
 x  y
 y  x""")
    else
        io_test(REPLMode,   A, """
2x2 Array{JuMP.Variable,2}:
 x  y
 y  x""")
        io_test(IJuliaMode, A, """
2x2 Array{JuMP.Variable,2}:
 x  y
 y  x""")
    end
end

facts("[print] basename keyword argument") do
    m = Model()
    @variable(m, x, basename="foo")
    @variable(m, y[1:3], basename=:bar)
    num = 123
    @variable(m, z[[:red,:blue]], basename="color_$num")
    @variable(m, v[1:2,1:2], SDP, basename=string("i","$num",num))
    @variable(m, w[1:3,1:3], Symmetric, basename="symm")

    io_test(REPLMode,   x, "foo")
    io_test(IJuliaMode, x, "foo")
    io_test(REPLMode,   y[2], "bar[2]")
    io_test(IJuliaMode, y[2], "bar_{2}")
    io_test(REPLMode,   z[:red], "color_123[red]")
    io_test(IJuliaMode, z[:red], "color_123_{red}")
    io_test(REPLMode,   v[2,1], "i123123[1,2]")
    io_test(IJuliaMode, v[2,1], "i123123_{1,2}")
    io_test(REPLMode,   w[1,3], "symm[1,3]")
    io_test(IJuliaMode, w[1,3], "symm_{1,3}")
end
