#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# print.jl
# All "pretty printers" for JuMP types, including IJulia LaTeX support
# where possible.
# Note that, in general, these methods are not "fast" - in fact, they are
# quite slow. This is OK, as they are normally only used on smaller models
# and don't appear in code that needs to be fast, but they should not affect
# the speed of code materially that doesn't use it.
#############################################################################

#############################################################################
#### type MODEL


# Default REPL
function Base.show(io::IO, m::Model)
    print(io, m.objSense == :Max ? "Maximization" : ((m.objSense == :Min && !isempty(m.obj)) ? "Minimization" : "Feasibility"))
    println(io, " problem with:")
    println(io, " * $(length(m.linconstr)) linear constraints")
    nquad = length(m.quadconstr)
    if nquad > 0
        println(io, " * $(nquad) quadratic constraints")
    end
    nlp = m.nlpdata
    if nlp != nothing && length(nlp.nlconstr) > 0
        println(io, " * $(length(nlp.nlconstr)) nonlinear constraints")
    end
    print(io, " * $(m.numCols) variables")
    nbin = sum(m.colCat .== :Bin)
    nint = sum(m.colCat .== :Int)
    nsc = sum(m.colCat .== :SemiCont)
    nsi = sum(m.colCat .== :SemiInt)
    varstr = {}
    nbin == 0 || push!(varstr, "$nbin binary")
    nint == 0 || push!(varstr, "$nint integer")
    nsc  == 0 || push!(varstr, "$nsc semicontinuous")
    nsi  == 0 || push!(varstr, "$nsi semi-integer")
    if isempty(varstr)
        println(io,)
    else
        println(io, ": $(join(varstr, ","))")
    end
    print(io, "Solver set to ")
    if isa(m.solver, UnsetSolver)
        solver = "Default"
    else
        solver = string(m.solver)
    end
    print(io, split(solver, "Solver")[1])
end