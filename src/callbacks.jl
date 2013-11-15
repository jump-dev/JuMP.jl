setlazycallback(m::Model, f::Function) = (m.lazycallback = f)

function registercallbacks(m::Model)
    if isa(m.lazycallback, Function)
        function lazycallback(d::MathProgCallbackData)
            m.colVal = cbgetmipsolution(d)
            m.lazycallback(d)
        end
        #try
            setlazycallback!(m.internalModel, lazycallback)
        #catch
        #    error("Solver does not support lazy callbacks")
        #end
    end
end

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
            # don't check for duplicates yet
            cbaddlazy!($cbdata, Cint[v.col for v in aff.vars], aff.coeffs, sense(constr), rhs(constr))
        end
    else
        error("Syntax error (ranged constraints not permitted in callbacks)")
    end
end

function addLazyConstraint(cbdata::MathProgCallbackData, constr::LinearConstraint)
    cbaddlazy!(cbdata, Cint[v.col for v in constr.terms.vars], constr.terms.coeffs, sense(constr), rhs(constr))
end

export addLazyConstraint, @addLazyConstraint, setlazycallback

