#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
    oldmac = Symbol(string("@",old))
    newmac = Symbol(string("@",new))
    s = string(oldmac," is deprecated, use ", newmac, " instead.")
    depwarn = :(Base.depwarn($s,$(quot(oldmac))))
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
@deprecate_macro addusercut usercut

@deprecate JuMPNLPEvaluator JuMP.NLPEvaluator
@deprecate addConstraint JuMP.addconstraint
@deprecate addCutCallback addcutcallback
@deprecate addHeuristicCallback addheuristiccallback
@deprecate addInfoCallback addinfocallback
@deprecate addLazyCallback addlazycallback
@deprecate addLazyConstraint JuMP.addlazyconstraint
@deprecate addSolution addsolution
@deprecate addToExpression JuMP.addtoexpr
@deprecate addUserCut JuMP.addusercut
@deprecate affToStr string
@deprecate buildInternalModel JuMP.build
@deprecate chgConstrRHS JuMP.setRHS
@deprecate conToStr string
@deprecate exprToStr string
@deprecate getCategory getcategory
@deprecate getConstraintBounds JuMP.constraintbounds
@deprecate getDual getdual
@deprecate getInternalModel internalmodel
@deprecate getLinearIndex linearindex
@deprecate getLower getlowerbound
@deprecate getName getname
@deprecate getObjective getobjective
@deprecate getObjectiveSense getobjectivesense
@deprecate getObjectiveValue getobjectivevalue
@deprecate getUpper getupperbound
@deprecate getValue getvalue
@deprecate getVar getvariable
@deprecate quadToStr string
@deprecate setCategory setcategory
@deprecate setLower setlowerbound
@deprecate setName setname
@deprecate setObjective JuMP.setobjective
@deprecate setObjectiveSense setobjectivesense
@deprecate setPrintHook JuMP.setprinthook
@deprecate setSolutionValue! setsolutionvalue
@deprecate setSolveHook JuMP.setsolvehook
@deprecate setSolver setsolver
@deprecate setUpper setupperbound
@deprecate setValue setvalue
@deprecate getvariable getindex
@deprecate getconstraint getindex

function registerNLFunction(args...;autodiff=false)
    Base.depwarn("registerNLFunction is deprecated, use JuMP.register instead.",:registerNLFunction)
    JuMP.register(args...,autodiff=autodiff)
end
export registerNLFunction
