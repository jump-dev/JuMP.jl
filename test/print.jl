#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/print.jl
# Testing for all pretty-printing-related functionality
# Defines:
#  - test_print_JuMPContainer()
#  - 
#############################################################################
using JuMP
using Base.Test

import JuMP.REPLMode, JuMP.IJuliaMode

const THROW_ERROR = true

# Helper function to test IO methods work correctly, and to provide
# useful outputs if they don't
function io_test(mode, obj, exp_str; repl=:both)
    if mode == REPLMode
        buf_print = IOBuffer()
        print(buf_print, obj)
        seek(buf_print, 0);
        print_str = readall(buf_print)
        expct_str = exp_str
        if !(print_str == exp_str) && repl != :show
            @show print_str
            @show expct_str
            THROW_ERROR && error()
        end 

        buf_show = IOBuffer()
        show(buf_show, obj)
        seek(buf_show, 0)
        show_str = readall(buf_show)
        expt_str = "\""*exp_str*"\""
        if !(show_str == expt_str) && repl != :print
            @show show_str
            @show expt_str
            THROW_ERROR && error()
        end
    else
        buf_display = IOBuffer()
        writemime(buf_display, "text/latex", obj)
        seek(buf_display,0)
        display_str = readall(buf_display)
        expectd_str = "\$\$ "*exp_str*" \$\$"
        if !(display_str == expectd_str)
            @show display_str
            @show expectd_str
            THROW_ERROR && error()
        end
    end
end


function test_print_JuMPContainer()
    println("  test_print_JuMPContainer")
    le, ge = JuMP.repl_leq, JuMP.repl_geq
    m = Model()
    
    #------------------------------------------------------------------
    # Test bound printing
    println("    bound printing")
    @defVar(m,      bnd_free[2:5])
    @defVar(m,      bnd_lowb[2:5] >= 2)
    @defVar(m,      bnd_high[2:5] <= 5)
    @defVar(m, 2 <= bnd_both[2:5] <= 5)
    @defVar(m,      bnd_difflo[i=2:5] >= i)
    @defVar(m,      bnd_diffup[i=2:5] <= i)
    @defVar(m, i <= bnd_diffbo[i=2:5] <= 2i)
    @defVar(m, i <= bnd_difflo_with_up[i=2:5] <= 5)
    @defVar(m, 2 <= bnd_diffup_with_lo[i=2:5] <= i)

    io_test(REPLMode, bnd_free, "bnd_free[i] free for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_lowb, "bnd_lowb[i] $ge 2 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_high, "bnd_high[i] $le 5 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_both, "2 $le bnd_both[i] $le 5 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_difflo, "bnd_difflo[i] $ge .. for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_diffup, "bnd_diffup[i] $le .. for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_diffbo, ".. $le bnd_diffbo[i] $le .. for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_difflo_with_up, ".. $le bnd_difflo_with_up[i] $le 5 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_diffup_with_lo, "2 $le bnd_diffup_with_lo[i] $le .. for all i in {2,3,4,5}")

    io_test(IJuliaMode, bnd_free, "bnd_free_{i} free \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_lowb, "bnd_lowb_{i} \\geq 2 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_high, "bnd_high_{i} \\leq 5 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_both, "2 \\leq bnd_both_{i} \\leq 5 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_difflo, "bnd_difflo_{i} \\geq .. \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_diffup, "bnd_diffup_{i} \\leq .. \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_diffbo, ".. \\leq bnd_diffbo_{i} \\leq .. \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_difflo_with_up, ".. \\leq bnd_difflo_with_up_{i} \\leq 5 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_diffup_with_lo, "2 \\leq bnd_diffup_with_lo_{i} \\leq .. \\quad\\forall i \\in \\{2,3,4,5\\}")

    #------------------------------------------------------------------
    # Test index set printing
    println("    index set printing")
    @defVar(m, rng_unit1[1:10])  # JuMPArray
    @defVar(m, rng_unit2[-2:3])  # JuMPArray
    @defVar(m, rng_unit3[[1:10]])  # JuMPDict
    @defVar(m, rng_step1[1:2:10])
    @defVar(m, rng_step2[-2:5:10])
    @defVar(m, rng_step3[1:5:3])
    @defVar(m, rng_step4[0:2:2])
    @defVar(m, arr_1[[:a,:b,:c]])
    @defVar(m, arr_2[[:a,1,"test"]])
    @defVar(m, arr_3[[:apple,:banana,:carrot,:diamonds]])
    @defVar(m, rng2_1[1:10,[:a,:b,:c]])
    @defVar(m, tri_1[i=1:3,j=i:3])
    @defVar(m, tri_2[i=1:3,j=-i])
    @defVar(m, tri_3[(i,j)={(i,i+2) for i in 1:5},k=i:j])
    
    io_test(REPLMode, rng_unit1, "rng_unit1[i] free for all i in {1,2..9,10}")
    io_test(REPLMode, rng_unit2, "rng_unit2[i] free for all i in {-2,-1..2,3}")
    io_test(REPLMode, rng_unit3, "rng_unit3[i] free for all i in {1,2..9,10}")
    io_test(REPLMode, rng_step1, "rng_step1[i] free for all i in {1,3..7,9}")
    io_test(REPLMode, rng_step2, "rng_step2[i] free for all i in {-2,3,8}")
    io_test(REPLMode, rng_step3, "rng_step3[i] free for all i in {1}")
    io_test(REPLMode, rng_step4, "rng_step4[i] free for all i in {0,2}")
    io_test(REPLMode, arr_1, "arr_1[i] free for all i in {a,b,c}")
    io_test(REPLMode, arr_2, "arr_2[i] free for all i in {a,1,test}")
    io_test(REPLMode, arr_3, "arr_3[i] free for all i in {apple,banana,carrot,diamonds}")
    io_test(REPLMode, rng2_1, "rng2_1[i,j] free for all i in {1,2..9,10}, j in {a,b,c}")
    io_test(REPLMode, tri_1, "tri_1[i,j] free for all i in {1,2,3}, j in {..}")
    io_test(REPLMode, tri_2, "tri_2[i,j] free for all i in {1,2,3}, j in {..}")
    io_test(REPLMode, tri_3, "tri_3[i,j] free for all i in {(1,3),(2,4)..(4,6),(5,7)}, j in {..}")

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
    io_test(IJuliaMode, tri_1, "tri_1_{i,j} free \\quad\\forall i \\in \\{1,2,3\\}, j \\in \\{..\\}")
    io_test(IJuliaMode, tri_2, "tri_2_{i,j} free \\quad\\forall i \\in \\{1,2,3\\}, j \\in \\{..\\}")
    io_test(IJuliaMode, tri_3, "tri_3_{i,j} free \\quad\\forall i \\in \\{(1,3),(2,4),\\dots,(4,6),(5,7)\\}, j \\in \\{..\\}")

    #------------------------------------------------------------------
    # Test category printing
    println("    category printing")
    @defVar(m, cat_bin[1:3], Bin)
    @defVar(m, 2 <= cat_int[1:3] <= 5, Int)
    @defVar(m, cat_semiint_both[2:3] >= 2, SemiInt)
    @defVar(m, i <= cat_semiint_difflow[i=2:3] <= 4, SemiInt)
    @defVar(m, 2 <= cat_semiint_diffup[i=2:3] <= i, SemiInt)
    @defVar(m, i <= cat_semiint_none[i=2:3] <= 2i, SemiInt)
    @defVar(m, cat_semicont_both[2:3] >= 2, SemiCont)
    @defVar(m, i <= cat_semicont_difflow[i=2:3] <= 4, SemiCont)
    @defVar(m, 2 <= cat_semicont_diffup[i=2:3] <= i, SemiCont)
    @defVar(m, i <= cat_semicont_none[i=2:3] <= 2i, SemiCont)

    io_test(REPLMode, cat_bin, "cat_bin[i] in {0,1} for all i in {1,2,3}")
    io_test(REPLMode, cat_int, "2 $le cat_int[i] $le 5, integer, for all i in {1,2,3}")
    io_test(REPLMode, cat_semiint_both, "cat_semiint_both[i] in {2..Inf} or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semiint_difflow, "cat_semiint_difflow[i] in {....4} or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semiint_diffup, "cat_semiint_diffup[i] in {2....} or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semiint_none, "cat_semiint_none[i] in {......} or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semicont_both, "cat_semicont_both[i] in [2,Inf] or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semicont_difflow, "cat_semicont_difflow[i] in [..,4] or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semicont_diffup, "cat_semicont_diffup[i] in [2,..] or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semicont_none, "cat_semicont_none[i] in [..,..] or {0} for all i in {2,3}")

    io_test(IJuliaMode, cat_bin, "cat_bin_{i} \\in \\{0,1\\} \\quad\\forall i \\in \\{1,2,3\\}")
    io_test(IJuliaMode, cat_int, "2 \\leq cat_int_{i} \\leq 5, \\in \\mathbb{Z}, \\quad\\forall i \\in \\{1,2,3\\}")
    io_test(IJuliaMode, cat_semiint_both, "cat_semiint_both_{i} \\in \\{2,\\dots,\\intfy\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semiint_difflow, "cat_semiint_difflow_{i} \\in \\{..,\\dots,4\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semiint_diffup, "cat_semiint_diffup_{i} \\in \\{2,\\dots,..\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semiint_none, "cat_semiint_none_{i} \\in \\{..,\\dots,..\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_both, "cat_semicont_both_{i} \\in \\[2,\\intfy\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_difflow, "cat_semicont_difflow_{i} \\in \\[..,4\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_diffup, "cat_semicont_diffup_{i} \\in \\[2,..\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_none, "cat_semicont_none_{i} \\in \\[..,..\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")

    #------------------------------------------------------------------
    # Tests for particular issues
    println("    issue testing")
    # Empty JuMPContainer printing (#124)
    @defVar(m, empty_free[1:0])
    io_test(REPLMode, empty_free, "empty_free (no indices)")
    io_test(IJuliaMode, empty_free, "empty_free (no indices)")
end

function test_print_SOS()
    println("  test_print_SOS")
    modS = Model()
    a = [1,2,3]
    @defVar(modS, x[1:3], Bin)
    addSOS1(modS, [a[i]x[i] for i in 1:3])    
    s1 = JuMP.SOSConstraint([x[i] for i in 1:3],
                            [a[i] for i in 1:3], :SOS1)
    io_test(REPLMode, s1, "SOS1: {1 x[1], 2 x[2], 3 x[3]}")
    io_test(IJuliaMode, s1, "SOS1: \\{1 x[1], 2 x[2], 3 x[3]\\}")

    b = [5,4,7,2,1]
    @defVar(modS, y[1:5], Bin)
    s2 = JuMP.SOSConstraint([y[i] for i in 1:5],
                            [b[i] for i in 1:5], :SOS2)
    io_test(REPLMode, s2, "SOS2: {5 y[1], 4 y[2], 7 y[3], 2 y[4], 1 y[5]}")
    io_test(IJuliaMode, s2, "SOS2: \\{5 y[1], 4 y[2], 7 y[3], 2 y[4], 1 y[5]\\}")
end

function test_print_Model()
    println("  test_print_Model")
    le, ge = JuMP.repl_leq, JuMP.repl_geq

    #------------------------------------------------------------------

    mod_1 = Model()
    @defVar(mod_1, a>=1)
    @defVar(mod_1, b<=1)
    @defVar(mod_1, -1<=c<=1)
    @defVar(mod_1, a1>=1,Int)
    @defVar(mod_1, b1<=1,Int)
    @defVar(mod_1, -1<=c1<=1,Int)
    @defVar(mod_1, x, Bin)
    @defVar(mod_1, y)
    @defVar(mod_1, z, Int)
    @defVar(mod_1, sos[1:3], Bin)
    @defVar(mod_1, 2 <= si <= 3, SemiInt)
    @defVar(mod_1, 2 <= sc <= 3, SemiCont)
    @setObjective(mod_1, Max, a - b + 2a1 - 10x)
    @addConstraint(mod_1, a + b - 10c - 2x + c1 <= 1)
    @addConstraint(mod_1, a*b <= 2)
    addSOS1(mod_1, [i*sos[i] for i in 1:3])

    io_test(REPLMode, mod_1, """
Max a - b + 2 a1 - 10 x
Subject to
 a + b - 10 c - 2 x + c1 $le 1
 a*b - 2 $le 0
 SOS1: {1 sos[1], 2 sos[2], 3 sos[3]}
 sos[i] in {0,1} for all i in {1,2,3}
 a $ge 1
 b $le 1
 -1 $le c $le 1
 a1 $ge 1, integer
 b1 $le 1, integer
 -1 $le c1 $le 1, integer
 x in {0,1}
 y free
 z free, integer
 si in {2..3} or {0}
 sc in [2,3] or {0}
""", repl=:print)

    io_test(IJuliaMode, mod_1, """
\\begin{alignat*}{1}\\max\\quad & a - b + 2 a1 - 10 x\\\\
\\text{Subject to} \\quad & a + b - 10 c - 2 x + c1 \\leq 1\\\\
 & a\\timesb - 2 \\leq 0\\\\
 & SOS1: \\{1 sos[1], 2 sos[2], 3 sos[3]\\}\\\\
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
\\end{alignat*}
""")

    #------------------------------------------------------------------

    mod_2 = Model()
    @defVar(mod_2, x[1:5])
    @addNLConstraint(mod_2, x[1]*x[2] == 1)
    @addNLConstraint(mod_2, x[3]*x[4] == 1)
    @addNLConstraint(mod_2, x[5]*x[1] == 1)
    @setNLObjective(mod_2, Min, x[1]*x[3])
    
    io_test(REPLMode, mod_2, """
Min (nonlinear expression)
Subject to
 3 nonlinear constraints
 x[i] free for all i in {1,2..4,5}
""", repl=:print)
    io_test(IJuliaMode, mod_2, """
\\begin{alignat*}{1}\\min\\quad & (nonlinear expression)\\\\
\\text{Subject to} \\quad & 3 nonlinear constraints\\\\
 & x_{i} free \\quad\\forall i \\in \\{1,2,\\dots,4,5\\}\\\\
\\end{alignat*}
""", repl=:print)
end

function test_print_expr()
    # Most of the expression logic is well covered by test/operator.jl
    # This is really just to check IJulia printing for expressions
    println("  test_print_expr")
    le, ge = JuMP.repl_leq, JuMP.repl_geq

    #------------------------------------------------------------------
    mod = Model()
    @defVar(mod, x[1:5])
    @defVar(mod, y[i=2:4,j=i:5])
    @defVar(mod, z)

    @addConstraint(mod, x[1] + 2*y[2,3] <= 3)
    io_test(REPLMode, mod.linconstr[end], "x[1] + 2 y[2,3] $le 3")
    io_test(IJuliaMode, mod.linconstr[end], "x_{1} + 2 y_{2,3} \\leq 3")
end

function test_print_Variable()
    println("  test_print_Variable")
    
    m = Model()

    @defVar(m, 0 <= x <= 2)
    
    @test    getName(x) == "x"
    io_test(REPLMode,   x, "x")
    io_test(IJuliaMode, x, "x")

    setName(x, "x2")
    @test    getName(x) == "x2"
    io_test(REPLMode,   x, "x2")
    io_test(IJuliaMode, x, "x2")

    setName(x, "")
    @test    getName(x) == "col_1"
    io_test(REPLMode,   x, "col_1")
    io_test(IJuliaMode, x, "col_1")

    @defVar(m, z[1:2,3:5])
    @test       getName(z[1,3]) == "z[1,3]"
    io_test(REPLMode,   z[1,3],    "z[1,3]")
    io_test(IJuliaMode, z[1,3],    "z_{1,3}")
    @test       getName(z[2,4]) == "z[2,4]"
    io_test(REPLMode,   z[2,4],    "z[2,4]")
    io_test(IJuliaMode, z[2,4],    "z_{2,4}")
    @test       getName(z[2,5]) == "z[2,5]"
    io_test(REPLMode,   z[2,5],    "z[2,5]")
    io_test(IJuliaMode, z[2,5],    "z_{2,5}")

    @defVar(m, w[3:9,["red","blue","green"]])
    @test    getName(w[7,"green"]) == "w[7,green]"
    io_test(REPLMode,   w[7,"green"], "w[7,green]")
    io_test(IJuliaMode, w[7,"green"], "w_{7,green}")

    rng = 2:5
    @defVar(m, v[rng,rng,rng,rng,rng,rng,rng])
    a_v = v[4,5,2,3,2,2,4]
    @test    getName(a_v) == "v[4,5,2,3,2,2,4]"
    io_test(REPLMode,   a_v, "v[4,5,2,3,2,2,4]")
    io_test(IJuliaMode, a_v, "v_{4,5,2,3,2,2,4}")
end

test_print_JuMPContainer()
test_print_SOS()
test_print_Model()
test_print_expr()
test_print_Variable()