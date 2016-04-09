#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# deprecated.jl
# List of deprecated macros and functions, mostly from the "big renaming"
macro deprecate_macro(old,new)
    oldmac = symbol(string("@",old))
    newmac = symbol(string("@",new))
    s = string(oldmac," is deprecated, use ", newmac, " instead.")
    if VERSION > v"0.5-"
        # backtraces are ok on 0.5
        depwarn = :(Base.depwarn($s,$(quot(oldmac))))
    else
        # backtraces are junk on 0.4
        depwarn = :(Base.warn_once($s))
    end
    @eval macro $old(args...)
        return Expr(:block, $depwarn, Expr(:macrocall, $(quot(newmac)), [esc(x) for x in args]...))
    end
    eval(Expr(:export,oldmac))
    return
end

@deprecate_macro defVar variable
@deprecate_macro defVars variables
@deprecate_macro setObjective objective
@deprecate_macro setNLObjective NLobjective
@deprecate_macro addConstraint constraint
@deprecate_macro addConstraints constraints
@deprecate_macro addNLConstraint NLconstraint
@deprecate_macro addNLConstraints NLconstraints
@deprecate_macro addSDPConstraint SDconstraint
@deprecate_macro defExpr expression
@deprecate_macro defNLExpr NLexpression
@deprecate_macro defNLParam NLparameter
@deprecate_macro defConstrRef constraintref

# callback macros
@deprecate_macro addLazyConstraint lazyconstraint
@deprecate_macro addUserCut usercut
