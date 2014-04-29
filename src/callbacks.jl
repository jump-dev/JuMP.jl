export setLazyCallback, setCutCallback, setHeuristicCallback
export setlazycallback
@Base.deprecate setlazycallback setLazyCallback
setLazyCallback(m::Model, f::Function) = (m.lazycallback = f)
setCutCallback(m::Model, f::Function) = (m.cutcallback = f)
setHeuristicCallback(m::Model, f::Function) = (m.heurcallback = f)

function registercallbacks(m::Model)
    if isa(m.lazycallback, Function)
        function lazycallback(d::MathProgCallbackData)
            state = cbgetstate(d)
            if state == :MIPSol
                cbgetmipsolution(d,m.colVal)
            else
                cbgetlpsolution(d,m.colVal)
            end
            m.lazycallback(d)
        end
        #try
            setlazycallback!(m.internalModel, lazycallback)
        #catch
        #  error("Solver does not support lazy callbacks")
        #end
    end
    if isa(m.cutcallback, Function)
        function cutcallback(d::MathProgCallbackData)
            state = cbgetstate(d)
            if state == :MIPSol  # This shouldn't happen right?
                println("Is this ever called?")
                cbgetmipsolution(d,m.colVal)
            else
                cbgetlpsolution(d,m.colVal)
            end
            m.cutcallback(d)
        end
        #
            setcutcallback!(m.internalModel, cutcallback)
        #
    end
    if isa(m.heurcallback, Function)
        function heurcallback(d::MathProgCallbackData)
            state = cbgetstate(d)
            if state == :MIPSol  # This shouldn't happen right?
                println("Is this ever called?")
                cbgetmipsolution(d,m.colVal)
            else
                cbgetlpsolution(d,m.colVal)
            end
            m.heurcallback(d)
        end
        #
            setheuristiccallback!(m.internalModel, heurcallback)
        #
    end

    # prepare storage for callbacks
    m.indexedVector = IndexedVector(Float64, m.numCols)
end


# TODO: Should this be somewhere else?
const sensemap = [:(<=) => '<', :(==) => '=', :(>=) => '>']


## Lazy constraints
export addLazyConstraint, @addLazyConstraint

macro addLazyConstraint(cbdata, x)
    cbdata = esc(cbdata)
    if (x.head != :comparison)
        error("Expected comparison operator in constraint $x")
    end
    if length(x.args) == 3 # simple comparison
        lhs = :($(x.args[1]) - $(x.args[3])) # move everything to the lhs
        quote
            aff = AffExpr()
            $(parseExpr(lhs, :aff, 1.0))
            constr = $(x.args[2])(aff,0)
            addLazyConstraint($cbdata, constr)
        end
    else
        error("Syntax error (ranged constraints not permitted in callbacks)")
    end
end

function addLazyConstraint(cbdata::MathProgCallbackData, constr::LinearConstraint)
    if length(constr.terms.vars) == 0
        cbaddlazy!(cbdata, Cint[], Float64[], sensemap[sense(constr)], rhs(constr))
        return
    end
    m::Model = constr.terms.vars[1].m
    indices, coeffs = merge_duplicates(Cint, constr.terms, m.indexedVector, m)
    cbaddlazy!(cbdata, indices, coeffs, sensemap[sense(constr)], rhs(constr))
end

## User cuts
export addUserCut, @addUserCut

macro addUserCut(cbdata, x)
    cbdata = esc(cbdata)
    if (x.head != :comparison)
        error("Expected comparison operator in constraint $x")
    end
    if length(x.args) == 3 # simple comparison
        lhs = :($(x.args[1]) - $(x.args[3])) # move everything to the lhs
        quote
            aff = AffExpr()
            $(parseExpr(lhs, :aff, 1.0))
            constr = $(x.args[2])(aff,0)
            addUserCut($cbdata, constr)
        end
    else
        error("Syntax error (ranged constraints not permitted in callbacks)")
    end
end

function addUserCut(cbdata::MathProgCallbackData, constr::LinearConstraint)
    if length(constr.terms.vars) == 0
        cbaddcut!(cbdata, Cint[], Float64[], sensemap[sense(constr)], rhs(constr))
        return
    end
    m::Model = constr.terms.vars[1].m
    indices, coeffs = merge_duplicates(Cint, constr.terms, m.indexedVector, m)
    cbaddcut!(cbdata, indices, coeffs, sensemap[sense(constr)], rhs(constr))
end

## User heuristic
export addSolution, setSolutionValue!

addSolution(cbdata::MathProgCallbackData) = cbaddsolution!(cbdata)
function setSolutionValue!(cbdata::MathProgCallbackData, v::Variable, x)
    cbsetsolutionvalue!(cbdata, convert(Cint, v.col), x)
end
