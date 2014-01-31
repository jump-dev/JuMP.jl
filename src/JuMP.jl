#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

import Base.getindex
import Base.setindex!
import Base.print
import Base.show

module JuMP

# Use the standard solver interface for LPs and MIPs
using MathProgBase
require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
using MathProgSolverInterface

importall Base

export
# Objects
    Model, Variable, AffExpr, QuadExpr, LinearConstraint, QuadConstraint, MultivarDict,
# Functions
    # Relevant to all
    print,show,
    # Model related
    getNumVars, getNumConstraints, getObjectiveValue, getObjective,
    getObjectiveSense, setObjectiveSense, writeLP, writeMPS, setObjective,
    addConstraint, addVar, addVars, solve, copy,
    # Variable
    setName, getName, setLower, setUpper, getLower, getUpper,
    getValue, setValue, getDual,
    # Expressions and constraints
    affToStr, quadToStr, conToStr, chgConstrRHS,
    
# Macros and support functions
    @addConstraint, @defVar, 
    @defConstrRef, @setObjective, addToExpression

include("JuMPDict.jl")
include("utils.jl")

###############################################################################
# Constants
const CONTINUOUS = 0
const INTEGER = 1
export CONTINUOUS, INTEGER

###############################################################################
# Model class
# Keeps track of all model and column info
type Model
    obj#::QuadExpr
    objSense::Symbol
    
    linconstr#::Vector{LinearConstraint}
    quadconstr
    
    # Column data
    numCols::Int
    colNames::Vector{String}
    colLower::Vector{Float64}
    colUpper::Vector{Float64}
    colCat::Vector{Int}

    # Solution data
    objVal
    colVal::Vector{Float64}
    redCosts::Vector{Float64}
    linconstrDuals::Vector{Float64}
    # internal solver model object
    internalModel
    # Solver+option object from MathProgBase
    solver::AbstractMathProgSolver
    # true if we haven't solved yet
    firstsolve::Bool
    # callbacks
    lazycallback
    cutcallback

    # JuMPDict list
    dictList::Vector
end

# Default constructor
function Model(sense::Symbol;lpsolver=MathProgBase.defaultLPsolver,mipsolver=MathProgBase.defaultMIPsolver,solver=nothing)
    Base.warn_once("Model(:$sense) syntax is deprecated. The sense should be passed to setObjective, e.g. @setObjective(model, :$sense, ...)")
    if lpsolver != MathProgBase.defaultLPsolver || mipsolver != MathProgBase.defaultMIPsolver
        error("lpsolver and mipsolver keywords have been merged. Use 'solver' instead, for example, Model(solver=ClpSolver())")
    end
    if solver == nothing
        # use default solvers
        Model(QuadExpr(),sense,LinearConstraint[], QuadConstraint[],
              0,String[],Float64[],Float64[],Int[],
              0,Float64[],Float64[],Float64[],nothing,MathProgBase.MissingSolver("",Symbol[]),true,
              nothing,nothing,JuMPDict[])
    else
        if !isa(solver,AbstractMathProgSolver)
            error("solver argument ($solver) must be an AbstractMathProgSolver")
        end
        # user-provided solver must support problem class
        Model(QuadExpr(),sense,LinearConstraint[], QuadConstraint[],
              0,String[],Float64[],Float64[],Int[],
              0,Float64[],Float64[],Float64[],nothing,solver,true,
              nothing,nothing,JuMPDict[])
    end
end

function Model(;solver=nothing,lpsolver=MathProgBase.defaultLPsolver,mipsolver=MathProgBase.defaultMIPsolver)
    if lpsolver != MathProgBase.defaultLPsolver || mipsolver != MathProgBase.defaultMIPsolver
        error("lpsolver and mipsolver keywords have been merged. Use 'solver' instead, for example, Model(solver=ClpSolver())")
    end
    
    if solver == nothing
        # use default solvers
        Model(QuadExpr(),:Min,LinearConstraint[], QuadConstraint[],
              0,String[],Float64[],Float64[],Int[],
              0,Float64[],Float64[],Float64[],nothing,MathProgBase.MissingSolver("",Symbol[]),true,
              nothing,nothing,JuMPDict[])
    else
        if !isa(solver,AbstractMathProgSolver)
            error("solver argument ($solver) must be an AbstractMathProgSolver")
        end
        # user-provided solver must support problem class
        Model(QuadExpr(),:Min,LinearConstraint[], QuadConstraint[],
              0,String[],Float64[],Float64[],Int[],
              0,Float64[],Float64[],Float64[],nothing,solver,true,
              nothing,nothing,JuMPDict[])
    end
end

# Getters/setters
getNumVars(m::Model) = m.numCols
getNumConstraints(m::Model) = length(m.linconstr)
getObjectiveValue(m::Model) = m.objVal
getObjectiveSense(m::Model) = m.objSense
function setObjectiveSense(m::Model, newSense::Symbol)
    if (newSense != :Max && newSense != :Min)
        error("Model sense must be :Max or :Min")
    end
    m.objSense = newSense
end

macro fill_names(m, cur, name, indexsets)
    m = esc(m)
    cur = esc(cur)
    name = esc(name)
    refcall = :( $m.colNames[$cur] )
    idxvars = {}
    idxsets = {}
    for s in indexsets.args
        if isa(s,Expr) && s.head == :(=)
            idxvar = s.args[1]
            idxset = s.args[2]
        else
            idxvar = gensym()
            idxset = s
        end
        push!(idxvars, idxvar)
        push!(idxsets, idxset)
    end
    tup = Expr(:tuple, [x for x in idxvars]...)
    code =  quote
                $(refcall) = $(name)*string($tup)
                # $(refcall) = $(name)
                $cur += 1
            end
    for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
        code =  quote
                    for $idxvar in $idxset
                        $code
                    end
                end
    end   
    println(code)
    return code
end

function fill_var_names(m)
    cnt, cur = 0, 1
    while cur <= m.numCols
        if m.colNames[cur] == ""
            cnt += 1
            dict = m.dictList[cnt]
            println("dict.indexsets=$(dict.indexsets)")
            println("cur=$cur")
            println("dict.name=$(dict.name)")
            @fill_names(m, cur, dict.name, dict.indexsets)
            println("cur=$cur")
        else
            cur += 1
        end
    end
end

# Pretty print
function print(io::IO, m::Model)
    (length(m.dictList) > 0) && fill_var_names(m)

    println(io, string(m.objSense," ",quadToStr(m.obj)))
    println(io, "Subject to: ")
    for c in m.linconstr
        println(io, conToStr(c))
    end
    for c in m.quadconstr
        println(io, conToStr(c))
    end
    for i in 1:m.numCols
        if m.colCat[i] == INTEGER && m.colLower[i] == 0 && m.colUpper[i] == 1
            print(io, (m.colNames[i] == "" ? string("_col",i) : m.colNames[i]))
            # println(" \u220a {0,1}")
            println(io, " binary")
        elseif m.colLower[i] == -Inf && m.colUpper[i] == Inf
            if m.colCat[i] == INTEGER
                print(io, (m.colNames[i] == "" ? string("_col",i) : m.colNames[i]))
                # println(io, "\u220a\u2124")
                println(io, " free integer")
            else
                print(io, (m.colNames[i] == "" ? string("_col",i) : m.colNames[i]))
                println(io, " free")
            end
        elseif m.colLower[i] == -Inf
            print(io, (m.colNames[i] == "" ? string("_col",i) : m.colNames[i]))
            print(io, " \u2264 $(m.colUpper[i])")
            # println(io, m.colCat[i] == INTEGER ? ", $(m.colNames[i])\u220a\u2124" : "")
            println(io, m.colCat[i] == INTEGER ? ", $(m.colNames[i]) integer" : "")
        elseif m.colUpper[i] == Inf
            print(io, (m.colNames[i] == "" ? string("_col",i) : m.colNames[i]))
            print(io, " \u2265 $(m.colLower[i])")
            # println(io, m.colCat[i] == INTEGER ? ", $(m.colNames[i])\u220a\u2124" : "")
            println(io, m.colCat[i] == INTEGER ? ", $(m.colNames[i]) integer" : "")
        else
            print(io, "$(m.colLower[i]) \u2264 ")
            print(io, (m.colNames[i] == "" ? string("_col",i) : m.colNames[i]))
            print(io, " \u2264 $(m.colUpper[i])")
            # println(io, m.colCat[i] == INTEGER ? ", $(m.colNames[i])\u220a\u2124" : "")
            println(io, m.colCat[i] == INTEGER ? ", $(m.colNames[i]) integer" : "")
        end
    end
end
show(io::IO, m::Model) = print(m.objSense == :Max ? "Maximization problem" :
                                                    "Minimization problem") 
                                                    # What looks good here?

# Deep copy the model
function copy(source::Model)
    
    dest = Model()
    dest.solver = source.solver  # The two models are linked by this
    dest.lazycallback = source.lazycallback
    dest.cutcallback = source.cutcallback
    
    # Objective
    dest.obj = copy(source.obj, dest)
    dest.objSense = source.objSense

    # Constraints
    dest.linconstr = [copy(c, dest) for c in source.linconstr]
    dest.quadconstr = [copy(c, dest) for c in source.quadconstr]

    # Variables
    dest.numCols = source.numCols
    dest.colNames = source.colNames[:]
    dest.colLower = source.colLower[:]
    dest.colUpper = source.colUpper[:]
    dest.colCat = source.colCat[:]

    return dest
end

###############################################################################
# Variable class
# Doesn't actually do much, just a pointer back to the model
type Variable
    m::Model
    col::Int
end

function Variable(m::Model,lower::Number,upper::Number,cat::Int,name::String)
    m.numCols += 1
    push!(m.colNames, name)
    push!(m.colLower, convert(Float64,lower))
    push!(m.colUpper, convert(Float64,upper))
    push!(m.colCat, cat)
    push!(m.colVal,NaN)
    return Variable(m, m.numCols)
end

Variable(m::Model,lower::Number,upper::Number,cat::Int) =
    Variable(m,lower,upper,cat,"")



# Name setter/getters
setName(v::Variable,n::String) = (v.m.colNames[v.col] = n)

function getName(v::Variable) 
    if length(v.m.colNames) > 0
        return (v.m.colNames[v.col] == "" ? string("_col",v.col) : v.m.colNames[v.col])
    end
    nothing
end

getName(m::Model, col) = (m.colNames[col] == "" ? string("_col",col) : m.colNames[col])
print(io::IO, v::Variable) = print(io, getName(v))
show(io::IO, v::Variable) = print(io, getName(v))

# Bound setter/getters
setLower(v::Variable,lower::Number) = (v.m.colLower[v.col] = convert(Float64,lower))
setUpper(v::Variable,upper::Number) = (v.m.colUpper[v.col] = convert(Float64,upper))
getLower(v::Variable) = v.m.colLower[v.col]
getUpper(v::Variable) = v.m.colUpper[v.col]

# Value setter/getter
function setValue(v::Variable, val::Number)
    v.m.colVal[v.col] = val
end

function getValue(v::Variable) 
    if length(v.m.colVal) < getNumVars(v.m)
        error("Variable values not available. Check that the model was properly solved.")
    end
    return v.m.colVal[v.col]
end

# Dual value (reduced cost) getter
function getDual(v::Variable) 
    if length(v.m.redCosts) < getNumVars(v.m)
        error("Variable bound duals (reduced costs) not available. Check that the model was properly solved and no integer variables are present.")
    end
    return v.m.redCosts[v.col]
end

###############################################################################
# Generic affine expression class
# Holds a vector of tuples (Var, Coeff)
type GenericAffExpr{CoefType,VarType}
    vars::Array{VarType,1}
    coeffs::Array{CoefType,1}
    constant::CoefType
end
print(io::IO, a::GenericAffExpr) = print(io, affToStr(a))
show(io::IO, a::GenericAffExpr) = print(io, affToStr(a))

typealias AffExpr GenericAffExpr{Float64,Variable}
AffExpr() = AffExpr(Variable[],Float64[],0.)

function setObjective(m::Model, a::AffExpr)
    Base.warn_once("Calling setObjective without specifying an objective sense is deprecated. Use setObjective(model, sense, expr) (or @setObjective(model, sense, expr)).")
    m.obj = QuadExpr()
    m.obj.aff = a
end

function setObjective(m::Model, sense::Symbol, a::AffExpr)
    setObjectiveSense(m, sense)
    m.obj = QuadExpr()
    m.obj.aff = a
end

function affToStr(a::AffExpr, showConstant=true)
    if length(a.vars) == 0
        if showConstant
            return string(a.constant)
        else
            return "0.0"
        end
    end

    # Get reference to model
    m = a.vars[1].m

    # Collect like terms
    indvec = IndexedVector(Float64,m.numCols)
    for ind in 1:length(a.vars)
        addelt(indvec, a.vars[ind].col, a.coeffs[ind])
    end

    elm = 0
    termStrings = Array(UTF8String, 2*length(a.vars))
    for i in 1:indvec.nnz
        idx = indvec.nzidx[i]
        if abs(indvec.elts[idx]) > 1e-20
            if elm == 0
                elm += 1
                termStrings[1] = "$(indvec.elts[idx]) $(getName(m,idx))"
            else 
                if indvec.elts[idx] < 0
                    termStrings[2*elm] = " - "
                else
                    termStrings[2*elm] = " + "
                end
                termStrings[2*elm+1] = "$(abs(indvec.elts[idx])) $(getName(m,idx))"
                elm += 1
            end
        end
    end

    if elm == 0
        ret = "0.0"
    else
        # And then connect them up with +s
        ret = join(termStrings[1:(2*elm-1)])
    end
    
    if abs(a.constant) >= 0.000001 && showConstant
        if a.constant < 0
            ret = string(ret, " - ", abs(a.constant))
        else
            ret = string(ret, " + ", a.constant)
        end
    end
    return ret
end

# Copy utility function, not exported
function copy(a::AffExpr, new_model::Model)
    return AffExpr([Variable(new_model, v.col) for v in a.vars],
                                 a.coeffs[:], a.constant)
end


###############################################################################
# QuadExpr class
# Holds a vector of tuples (Var, Var, Coeff), as well as an AffExpr
type QuadExpr
    qvars1::Vector{Variable}
    qvars2::Vector{Variable}
    qcoeffs::Vector{Float64}
    aff::AffExpr
end

QuadExpr() = QuadExpr(Variable[],Variable[],Float64[],AffExpr())

function setObjective(m::Model, q::QuadExpr)
    Base.warn_once("Calling setObjective without specifying an objective sense is deprecated. Use setObjective(model, sense, expr).")
    m.obj = q
end

function setObjective(m::Model, sense::Symbol, q::QuadExpr)
    m.obj = q
    setObjectiveSense(m, sense)
end


print(io::IO, q::QuadExpr) = print(io, quadToStr(q))
show(io::IO, q::QuadExpr) = print(io, quadToStr(q))

function quadToStr(q::QuadExpr)
    if length(q.qvars1) == 0
        return affToStr(q.aff)
    end

    m::Model = q.qvars1[1].m
    # canonicalize and merge duplicates
    for ind in 1:length(q.qvars1)
            if q.qvars2[ind].col < q.qvars1[ind].col
                    q.qvars1[ind],q.qvars2[ind] = q.qvars2[ind],q.qvars1[ind]
            end
    end
    Q = sparse([v.col for v in q.qvars1], [v.col for v in q.qvars2], q.qcoeffs)
    I,J,V = findnz(Q)

    termStrings = Array(UTF8String, 2*nnz(Q))
    if nnz(Q) > 0
        if V[1] < 0
            termStrings[1] = "-"
        else
            termStrings[1] = ""
        end
        for ind in 1:nnz(Q)
            if ind >= 2
                if V[ind] < 0
                    termStrings[2*ind-1] = " - "
                else 
                    termStrings[2*ind-1] = " + "
                end
            end
            x = Variable(m,I[ind])
            if I[ind] == J[ind]
                # Squared term
                termStrings[2*ind] = string(abs(V[ind])," ",
                                            getName(x),"²")
            else
                # Normal term
                y = Variable(m,J[ind])
                termStrings[2*ind] = string(abs(V[ind])," ",
                                            getName(x),"*",
                                            getName(y))
            end
        end
    end
    ret = join(termStrings)

    if q.aff.constant == 0 && length(q.aff.vars) == 0
        return ret
    else
        aff = affToStr(q.aff)
        if aff[1] == '-'
            return string(ret, " - ", aff[2:end])
        else
            return string(ret, " + ", aff)
        end
    end
end

# Copy utility function, not exported
function copy(q::QuadExpr, new_model::Model)
    return QuadExpr([Variable(new_model, v.col) for v in q.qvars1],
                    [Variable(new_model, v.col) for v in q.qvars2],
                    q.qcoeffs[:], copy(q.aff, new_model))
end

##########################################################################
# JuMPConstraint
# abstract base for constraint types
abstract JuMPConstraint

##########################################################################
# LinearConstraint class
# An affine expression with lower bound (possibly -Inf) and upper bound (possibly Inf).
type LinearConstraint <: JuMPConstraint
    terms::AffExpr
    lb::Float64
    ub::Float64
end

if VERSION.major == 0 && VERSION.minor < 3
    LinearConstraint(terms::AffExpr,lb::Number,ub::Number) =
        LinearConstraint(terms,float(lb),float(ub))
end


function addConstraint(m::Model, c::LinearConstraint)
    push!(m.linconstr,c)
    if !m.firstsolve
        # TODO: we don't check for duplicates here
        try
            addconstr!(m.internalModel,[v.col for v in c.terms.vars],c.terms.coeffs,c.lb,c.ub)
        catch
            Base.warn_once("Solver does not appear to support adding constraints to an existing model. Hot-start is disabled.")
            m.firstsolve = true
        end
    end
    return ConstraintRef{LinearConstraint}(m,length(m.linconstr))
end

print(io::IO, c::LinearConstraint) = print(io, conToStr(c))
show(io::IO, c::LinearConstraint) = print(io, conToStr(c))

function sense(c::LinearConstraint) 
    if c.lb != -Inf
        if c.ub != Inf
            if c.ub == c.lb
                return :(==)
            else
                return :range
            end
        else
                return :>=
        end
    else
        @assert c.ub != Inf
        return :<=
    end
end

function rhs(c::LinearConstraint)
    s = sense(c)
    @assert s != :range
    if s == :<=
        return c.ub
    else
        return c.lb
    end
end

function conToStr(c::LinearConstraint)
    s = sense(c)
    if s == :range
        return string(c.lb," <= ",affToStr(c.terms,false)," <= ",c.ub)
    else
        return string(affToStr(c.terms,false)," ",s," ",rhs(c))
    end
end

# Copy utility function, not exported
function copy(c::LinearConstraint, new_model::Model)
    return LinearConstraint(copy(c.terms, new_model), c.lb, c.ub)
end

##########################################################################
# QuadConstraint class
# An quadratic constraint. Right-hand side is implicitly taken to be zero; 
# constraint is stored in the included QuadExpr.
type QuadConstraint <: JuMPConstraint
    terms::QuadExpr
    sense::Symbol
end

function addConstraint(m::Model, c::QuadConstraint)
    push!(m.quadconstr,c)
    if !m.firstsolve
        # we don't (yet) support hot-starting QCQP solutions
        m.firstsolve = true
    end
    return ConstraintRef{QuadConstraint}(m,length(m.quadconstr))
end

print(io::IO, c::QuadConstraint) = print(io, conToStr(c))
show(io::IO, c::QuadConstraint)  = print(io, conToStr(c))

conToStr(c::QuadConstraint) = string(quadToStr(c.terms), " ", c.sense, " 0")

# Copy utility function, not exported
function copy(c::QuadConstraint, new_model::Model)
    return QuadConstraint(copy(c.terms, new_model), c.sense)
end

##########################################################################
# ConstraintRef
# Reference to a constraint for retrieving solution info
immutable ConstraintRef{T<:JuMPConstraint}
    m::Model
    idx::Int
end

function getDual(c::ConstraintRef{LinearConstraint}) 
    if length(c.m.linconstrDuals) != getNumConstraints(c.m)
        error("Dual solution not available. Check that the model was properly solved and no integer variables are present.")
    end
    return c.m.linconstrDuals[c.idx]
end

function chgConstrRHS(c::ConstraintRef{LinearConstraint}, rhs::Number)
    constr = c.m.linconstr[c.idx]
    sen = sense(constr)
    if sen == :range
        error("Modifying range constraints is currently unsupported.")
    elseif sen == :(==)
        constr.lb = float(rhs)
        constr.ub = float(rhs)
    elseif sen == :>=
        constr.lb = float(rhs)
    else
        @assert sen == :<=
        constr.ub = float(rhs)
    end
end

print(io::IO, c::ConstraintRef{LinearConstraint}) = print(io, conToStr(c.m.linconstr[c.idx]))
print(io::IO, c::ConstraintRef{QuadConstraint}) = print(io, conToStr(c.m.quadconstr[c.idx]))
show{T}(io::IO, c::ConstraintRef{T}) = print(io, c)

# add variable to existing constraints
function Variable(m::Model,lower::Number,upper::Number,cat::Int,objcoef::Number,
    constraints::Vector{ConstraintRef{LinearConstraint}},coefficients::Vector{Float64};
    name::String="")
        
    v = Variable(m, lower, upper, cat, name)
    # add to existing constraints
    @assert length(constraints) == length(coefficients)
    for i in 1:length(constraints)
        c::LinearConstraint = m.linconstr[constraints[i].idx]
        coef = coefficients[i]
        push!(c.terms.vars,v)
        push!(c.terms.coeffs,coef)
    end
    push!(m.obj.aff.vars, v)
    push!(m.obj.aff.coeffs,objcoef)

    if !m.firstsolve
        try
            addvar!(m.internalModel,Int[c.idx for c in constraints],coefficients,float(lower),float(upper),float(objcoef))
        catch
            Base.warn_once("Solver does not appear to support adding variables to an existing model. Hot-start is disabled.")
            m.firstsolve = true
        end
    end

    return v
end

##########################################################################
# Operator overloads
include("operators.jl")
# Writers - we support MPS (MILP + QuadObj), LP (MILP)
include("writers.jl")
# Solvers
include("solvers.jl")
# Macros - @defVar, sum{}, etc.
include("macros.jl")

include("callbacks.jl")

##########################################################################
end
