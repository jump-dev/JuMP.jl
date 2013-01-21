

# extracts left operands of *
# a*b*c -> a*b
function timescoef(x::Expr)
    if x.head != :call || x.args[1] != :*
        error("In expression $x expected multiplication")
    end
    x2 = copy(x)
    # just delete last argument
    delete!(x2.args,length(x2.args))
    return x2
end
# extracts last operand of *
# a*b*c -> c
function timesvar(x::Expr)
    if x.head != :call || x.args[1] != :*
        error("In expression $x expected multiplication")
    end
    return x.args[end]
end

# returns expression for array of the coefficient terms in the comprehension
function comprehensioncoef(x::Expr)
    if x.head != :comprehension
        error("In expression $x expected comprehension")
    end
    comp = expr(symbol("typed-comprehension"),{Float64,timescoef(x.args[1]),x.args[2:end]...})
    containers = Expr[ ex.args[2] for ex in x.args[2:end] ]
    out = quote l = 1 end
    for c in containers
        push!(out.args,:(l *= length($c)))
    end
    push!(out.args,:(reshape($comp,(l,))))
    return out
end

# returns expression for array of the variable terms in the comprehension
function comprehensionvar(x::Expr)
    if x.head != :comprehension
        error("In expression $x expected comprehension")
    end
    comp = expr(symbol("typed-comprehension"),{Variable,timesvar(x.args[1]),x.args[2:end]...})
    containers = Expr[ ex.args[2] for ex in x.args[2:end] ]
    out = quote l = 1 end
    for c in containers
        push!(out.args,:(l *= length($c)))
    end
    push!(out.args,:(reshape($comp,(l,))))
    return out
end

# parses top-level expression and returns expression for array of the coefficient terms
function topcoef(x::Expr)
    if x.head == :call && x.args[1] == :+
        return expr(:vcat,map(topcoef,x.args[2:end]))
    elseif x.head == :comprehension
        return comprehensioncoef(x)
    elseif x.head == :call && x.args[1] == :*
        println("timescoef of $x is $(timescoef(x))")
        return timescoef(x)
    else
        error("Unable to parse expression $x")
    end
end

function topvar(x::Expr)
    if x.head == :call && x.args[1] == :+
        return expr(:vcat,map(topvar,x.args[2:end]))
    elseif x.head == :comprehension
        return comprehensionvar(x)
    elseif x.head == :call && x.args[1] == :*
        return timesvar(x)
    else
        error("Unable to parse expression $x")
    end
end

macro sumExpr(x)
    esc(:(AffExpr($(topvar(x)),convert(Vector{Float64},$(topcoef(x))),0.)))
end
