using Base.Meta

# extracts left operands of *
# a*b*c -> a*b
function timescoef(x::Expr)
    if x.head != :call || x.args[1] != :*
        error("In expression $x expected multiplication")
    end
    x2 = copy(x)
    # just delete last argument
    splice!(x2.args,length(x2.args))
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

function addToExpression(aff::AffExpr,c::Number,x::Variable)
    push!(aff.vars,x)
    push!(aff.coeffs,c)
end

function addToExpression(aff::AffExpr,c::Number,x::Number)
    aff.constant += c*x
end

function addToExpression(aff::AffExpr,c::Number,x::AffExpr)
    append!(aff.vars, x.vars)
    append!(aff.coeffs, c*x.coeffs)
    aff.constant += c*x.constant
end

function parseCurly(x::Expr, aff::Symbol, constantCoef)
    if (x.args[1] != :sum)
        error("Expected sum outside curly braces")
    end
    if length(x.args) < 3
        error("Need at least two arguments for sum")
    end

    # we have a filter condition
    if isexpr(x.args[2],:parameters)
        cond = x.args[2]
        if length(cond.args) != 1
            error("No commas after semicolon allowed in sum expression, use && for multiple conditions")
        end
        # generate inner loop code first and then wrap in for loops
        code = quote
            if $(esc(cond.args[1]))
                $(parseExpr(x.args[3], aff, constantCoef))
            end
        end
        for level in length(x.args):-1:4
            code = :(
            for $(esc(x.args[level].args[1])) in $(esc(x.args[level].args[2]))
                $code
            end)
        end
    else # no condition
        code = parseExpr(x.args[2], aff, constantCoef)
        for level in length(x.args):-1:3
            code = :(
            for $(esc(x.args[level].args[1])) in $(esc(x.args[level].args[2]))
                $code
            end)
        end
        len = :len
        # precompute the number of elements to add
        # this is unncessary if we're just summing constants
        preblock = :($len += length($(esc(x.args[length(x.args)].args[2]))))
        for level in (length(x.args)-1):-1:3
            preblock = Expr(:for, esc(x.args[level]),preblock)
        end
        preblock = :($len = 0; $preblock;
            sizehint($aff.vars,length($aff.vars)+$len);
            sizehint($aff.coeffs,length($aff.coeffs)+$len))
        code = :($preblock;$code)
    end

    
    return code
end

function parseExpr(x, aff::Symbol, constantCoef)
    if !isa(x,Expr)
        # at the lowest level
        :(addToExpression($aff, $(esc(constantCoef)), $(esc(x))))
    else
        if x.head == :call && x.args[1] == :+
            Expr(:block,[parseExpr(arg,aff,constantCoef) for arg in x.args[2:end]]...)
        elseif x.head == :call && x.args[1] == :-
            if length(x.args) == 2 # unary subtraction
                parseExpr(x.args[2], aff, :(-1.0*$constantCoef))
            else # a - b - c ...
                Expr(:block,vcat(parseExpr(x.args[2], aff, constantCoef),
                     {parseExpr(arg, aff, :(-1.0*$constantCoef)) for arg in x.args[3:end]})...)
            end
        elseif x.head == :call && x.args[1] == :*
            coef = timescoef(x)
            var = timesvar(x)
            parseExpr(var, aff, :($coef*$constantCoef))
        elseif x.head == :call && x.args[1] == :/
            @assert length(x.args) == 3
            numerator = x.args[2]
            denom = x.args[3]
            parseExpr(numerator, aff, :((1/$denom)*$constantCoef))
        elseif x.head == :curly
            parseCurly(x,aff,constantCoef)
        else # at lowest level?
            :(addToExpression($aff, $(esc(constantCoef)), $(esc(x))))
        end
    end
end

macro addConstraint(m, x)
    m = esc(m)
    if (x.head != :comparison)
        error("Expected comparison operator in constraint $x")
    end
    if length(x.args) == 3 # simple comparison
        lhs = :($(x.args[1]) - $(x.args[3])) # move everything to the lhs
        quote
            aff = AffExpr()
            $(parseExpr(lhs, :aff, 1.0))
            addConstraint($m, $(x.args[2])(aff,0) )
        end
    else
        # ranged row
        if length(x.args) != 5 || x.args[2] != :<= || x.args[4] != :<=
            error("Only ranged rows of the form lb <= expr <= ub are supported")
        end
        lb = x.args[1]
        ub = x.args[5]
        quote
            aff = AffExpr()
            if !isa($(esc(lb)),Number)
                error(string("Expected ",$lb," to be a number"))
            elseif !isa($(esc(ub)),Number)
                error(string("Expected ",$ub," to be a number"))
            end
            $(parseExpr(x.args[3],:aff,1.0))
            addConstraint($m, 
                LinearConstraint(aff,$(esc(lb))-aff.constant,
                    $(esc(ub))-aff.constant))
        end
    end
end

macro setObjective(m, args...)
    if length(args) == 1
        x = args[1]
        quote
            aff = AffExpr()
            $(parseExpr(x, :aff, 1.0))
            setObjective($(esc(m)), aff)
        end
    else
        @assert length(args) == 2
        sense, x = args
        if sense == :Min || sense == :Max
            sense = Expr(:quote,sense)
        end
        quote
            aff = AffExpr()
            $(parseExpr(x, :aff, 1.0))
            setObjective($(esc(m)), $(esc(sense)), aff)
        end
    end
end
        

macro defVar(m, x, extra...)
    m = esc(m)
    if isexpr(x,:comparison)
        # we have some bounds
        if x.args[2] == :>=
            if length(x.args) == 5
                error("Use the form lb <= var <= ub instead of ub >= var >= lb")
            end
            @assert length(x.args) == 3
            # lower bounds, no upper
            lb = esc(x.args[3])
            ub = Inf
            var = x.args[1]
        elseif x.args[2] == :<=
            if length(x.args) == 5
                # lb <= x <= u
                lb = esc(x.args[1])
                if (x.args[4] != :<=)
                    error("Expected <= operator")
                end
                ub = esc(x.args[5])
                var = x.args[3]
            else
                # x <= u
                ub = esc(x.args[3])
                lb = -Inf
                var = x.args[1]
            end
        end
    else
        var = x
        lb = -Inf
        ub = Inf
    end
    t = JuMP.CONTINUOUS
    if length(extra) > 0
        gottype = 0
        if extra[1] == :Int || extra[1] == :Bin
            gottype = 1
            if extra[1] == :Int
                t = JuMP.INTEGER
            else
                if lb != -Inf || ub != Inf
                    error("Bounds may not be specified for binary variables. These are always taken to have a lower bound of 0 and upper bound of 1.")
                end
                t = JuMP.INTEGER
                lb = 0.0
                ub = 1.0
            end
        end
        if length(extra) - gottype == 3
            # adding variable to existing constraints
            objcoef = esc(extra[1+gottype])
            cols = esc(extra[2+gottype])
            coeffs = esc(extra[3+gottype])
            if !isa(var,Symbol)
                error("Cannot create multiple variables when adding to existing constraints")
            end
            return quote
                $(esc(var)) = Variable($m,$lb,$ub,$t,$objcoef,$cols,$coeffs,name=$(string(var)))
                nothing
            end
        elseif length(extra) - gottype != 0
            error("Syntax error in defVar")
        end
    end

    #println("lb: $lb ub: $ub var: $var")      
    if isa(var,Symbol)
        # easy case
        return quote
            $(esc(var)) = Variable($m,$lb,$ub,$t,$(string(var)))
            nothing
        end
    else
        if !isexpr(var,:ref)
            error("Syntax error: Expected $var to be of form var[...]")
        end
        varname = esc(var.args[1])
        idxvars = {}
        idxsets = {}
        refcall = Expr(:ref,varname)
        for s in var.args[2:end]
            if isa(s,Expr) && s.head == :(=)
                idxvar = s.args[1]
                idxset = esc(s.args[2])
            else
                idxvar = gensym()
                idxset = esc(s)
            end
            push!(idxvars, idxvar)
            push!(idxsets, idxset)
            push!(refcall.args, esc(idxvar))
        end
        tup = Expr(:tuple,[esc(x) for x in idxvars]...)
        code = :( $(refcall) = Variable($m, $lb, $ub, $t, $(string(var.args[1]))*string($tup) ) )
        for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
            code = quote
                for $(esc(idxvar)) in $idxset
                    $code
                end
            end
        end
        
        mac = Expr(:macrocall,symbol("@gendict"),varname,:Variable,idxsets...)
        code = quote 
            $mac
            $code
            nothing
        end
        return code
    end
end

macro defConstrRef(var)
    
    if isa(var,Symbol)
        # easy case
        return esc(:(local $var))
    else
        if !isexpr(var,:ref)
            error("Syntax error: Expected $var to be of form var[...]")
        end
        
        varname = var.args[1]
        idxsets = var.args[2:end]
                
        mac = Expr(:macrocall,symbol("@gendict"),varname,:ConstraintRef,idxsets...)
        code = quote 
            $(esc(mac))
            nothing
        end
        return code
    end
end
