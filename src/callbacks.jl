using Gurobi

setmipsolcallback(m::Model, f::Function) = (m.mipsolcallback = f)

function registergurobicallback(m::Model, grb::GurobiSolver)
    if !isa(mipsolcallback, Function)
        return
    end
    function gurobicallback(cbdata::CallbackData, where::Cint)
        if where == CB_MIPSOL
            m.colVal = cbget_mipsol_sol(cbdata, where)
            m.mipsolcallback(cbdata)
        end
    end

    set_callback_func!(grb.inner, gurobicallback)

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
            sens = sense(constr)
            local csense::Char
            if sens == :(==)
                csense = '='
            elseif sens == :(<=)
                csense = '<'
            else
                csense = '>'
            end
            cblazy($cbdata, Cint[v.col for v in aff.vars], aff.coeffs, csense, rhs(constr))
        end
    else
        error("Syntax error (ranged constraints not permitted in callbacks)")
    end
end

function addLazyConstraint(cbdata::CallbackData, constr::Constraint)
    sens = sense(constr)
    local csense::Char
    if sens == :(==)
        csense = '='
    elseif sens == :(<=)
        csense = '<'
    else
        csense = '>'
    end
    cblazy(cbdata, Cint[v.col for v in constr.terms.vars], constr.terms.coeffs, csense, rhs(constr))


end

export addLazyConstraint, @addLazyConstraint, setmipsolcallback

