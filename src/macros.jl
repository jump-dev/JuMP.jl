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

addToExpression(aff, c, x) = error("Cannot construct an affine expression with a term of type $(typeof(x))")

function parseCurly(x::Expr, aff::Symbol, constantCoef)
    if !(x.args[1] == :sum || x.args[1] == :∑ || x.args[1] == :Σ) # allow either N-ARY SUMMATION or GREEK CAPITAL LETTER SIGMA
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

macro addConstraint(m, x, extra...)
    m = esc(m)
    if length(extra) == 1
        c, x = x, extra[1]
        if (x.head != :comparison)
            error("Expected comparison operator in constraint $x")
        end
        if isa(c,Symbol)
            refcall = esc(c)
        else
            if !isexpr(c,:ref)
                error("Syntax error: Expected $c to be of form constr[...]")
            end
            cname = esc(c.args[1])
            idxvars = {}
            idxsets = {}
            refcall = Expr(:ref,cname)
            for s in c.args[2:end]
                if isa(s,Expr) && (s.head == :(=) || s.head == :in)
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
        end
        if length(x.args) == 3 # simple comparison
            lhs = :($(x.args[1]) - $(x.args[3])) # move everything to the lhs
            code = quote
                aff = AffExpr()
                $(parseExpr(lhs, :aff, 1.0))
                $(refcall) = addConstraint($m, $(x.args[2])(aff,0) )
            end
        else
            # ranged row
            if length(x.args) != 5 || (x.args[2] != :<= && x.args[2] != :≤) || (x.args[4] != :<= && x.args[4] != :≤)
                error("Only ranged rows of the form lb <= expr <= ub are supported")
            end
            lb = x.args[1]
            ub = x.args[5]
            code = quote
                aff = AffExpr()
                if !isa($(esc(lb)),Number)
                    error(string("Expected ",$lb," to be a number"))
                elseif !isa($(esc(ub)),Number)
                    error(string("Expected ",$ub," to be a number"))
                end
                $(parseExpr(x.args[3],:aff,1.0))
                $(refcall) = addConstraint($m, 
                    LinearConstraint(aff,$(esc(lb))-aff.constant,
                        $(esc(ub))-aff.constant))
            end
        end
        if isa(c,Symbol)
            # easy case
            return code
        else  
            for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
                code = quote
                    for $(esc(idxvar)) in $idxset
                        $code
                    end
                end
            end
            mac = Expr(:macrocall,symbol("@gendict"),cname,:(ConstraintRef{LinearConstraint}),idxsets...)
            return quote 
                $mac
                $code
                nothing
            end
        end
    elseif length(extra) == 0 # the old @addConstraint
        if (x.head != :comparison)
            error("Expected comparison operator in constraint $x")
        end
        if length(x.args) == 3 # simple comparison
            lhs = :($(x.args[1]) - $(x.args[3])) # move everything to the lhs
            return quote
                aff = AffExpr()
                $(parseExpr(lhs, :aff, 1.0))
                addConstraint($m, $(x.args[2])(aff,0) )
            end
        else
            # ranged row
            if length(x.args) != 5 || (x.args[2] != :<= && x.args[2] != :≤) || (x.args[4] != :<= && x.args[4] != :≤)
                error("Only ranged rows of the form lb <= expr <= ub are supported")
            end
            lb = x.args[1]
            ub = x.args[5]
            return quote
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
    else
        error("Too many arguments to addConstraint")
    end
end

macro setObjective(m, sense, x)
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    quote
        aff = AffExpr()
        $(parseExpr(x, :aff, 1.0))
        setObjective($(esc(m)), $(esc(sense)), aff)
    end
end
        

macro defVar(m, x, extra...)
    m = esc(m)
    # Identify the variable bounds. Four (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", or just plain "x"
    if isexpr(x,:comparison)
        # We have some bounds
        if x.args[2] == :>= || x.args[2] == :≥
            # x >= lb
            var = x.args[1]
            length(x.args) == 5 &&
                error("in @defVar ($var): use the form lb <= $var <= ub instead of ub >= $var >= lb")
            @assert length(x.args) == 3
            lb = esc(x.args[3])
            ub = Inf
        elseif x.args[2] == :<= || x.args[2] == :≤
            if length(x.args) == 5
                # lb <= x <= u
                var = x.args[3]
                (x.args[4] != :<= && x.args[4] != :≤) &&
                    error("in @defVar ($var): expected <= operator after variable name")
                lb = esc(x.args[1])
                ub = esc(x.args[5])
            else
                # x <= ub
                var = x.args[1]
                # NB: May also be lb <= x, which we do not support
                #     We handle this later in the macro
                @assert length(x.args) == 3
                ub = esc(x.args[3])
                lb = -Inf
            end
        else
            # Its a comparsion, but not using <= ... <=
            error("in @defVar ($(string(x))): use the form lb <= ... <= ub")
        end
    else
        # No bounds provided - free variable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        var = x
        lb = -Inf
        ub = Inf
    end

    # Determine variable type (if present), as well as variables being 
    # added as complete columns. 
    # Types: default is continuous (reals), alternatives are Int and Bin.
    # ColGen: format is @defVar(..., [type], objcoef, constrrefs, values)
    t = JuMP.CONTINUOUS
    if length(extra) > 0
        gottype = 0
        # Try to detect an acceptable variable type
        if extra[1] == :Int || extra[1] == :Bin
            gottype = 1
            if extra[1] == :Int
                t = JuMP.INTEGER
            else
                # Bin is internally just an integer variable
                # So if Bin, either no bounds at all, or must be 0 <= x <= 1
                if (lb != -Inf || ub != Inf) && !(lb == 0.0 && ub == 1.0)
                    error("in @defVar ($var): bounds other than [0, 1] may not be specified for binary variables.\nThese are always taken to have a lower bound of 0 and upper bound of 1.")
                end
                t = JuMP.INTEGER
                lb = 0.0
                ub = 1.0
            end
        end
        # Try to detect the case where someone meant to do Int or Bin,
        # but did something else instead
        length(extra) == 1 && gottype == 0 &&
            error("in @defVar ($var): provided variable type must be Int or Bin")
        # Handle the column generation functionality
        if length(extra) - gottype == 3
            !isa(var,Symbol) &&
                error("in @defVar ($var): can only create one variable at a time when adding to existing constraints.")

            objcoef = esc(extra[1+gottype])
            cols    = esc(extra[2+gottype])
            coeffs  = esc(extra[3+gottype])
            return quote
                $(esc(var)) = Variable($m,$lb,$ub,$t,$objcoef,$cols,$coeffs,name=$(string(var)))
                nothing
            end
        end
        gottype == 0 &&
            error("in @defVar ($var): syntax error")
    end

    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    if isa(var,Symbol)
        # Easy case - a single variable
        return quote
            $(esc(var)) = Variable($m,$lb,$ub,$t,$(string(var)))
        end
    else
        # An indexed set of variables
        # Check for malformed expression
        !isexpr(var,:ref) &&
            error("in @defVar ($var): expected $var to be of form var[...]")

        varname = esc(var.args[1])
        idxvars = {}
        idxsets = {}
        refcall = Expr(:ref,varname)
        # Iterate over each index set s
        for s in var.args[2:end]
            # Is the user providing an index variable, e.g. i=1:5?
            if isa(s,Expr) && (s.head == :(=) || s.head == :in)
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
        
        tup = Expr(:tuple, [esc(x) for x in idxvars]...)
        code = :( $(refcall) = Variable($m, $lb, $ub, $t) )
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
            push!($(m).dictList, $varname)
            $varname
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

macro setNLObjective(m, sense, x)
    m = esc(m)
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    quote
        initNLP($m)
        setObjectiveSense($m, $(esc(sense)))
        ex = @processNLExpr($(esc(x)))
        $m.nlpdata.nlobj = ex
        $m.obj = QuadExpr()
    end
end 


macro addNLConstraint(m, x)
    m = esc(m)
    if (x.head != :comparison)
        error("Expected comparison operator in constraint $x")
    end
    if length(x.args) == 3 # simple comparison
        lhs = :($(x.args[1]) - $(x.args[3])) # move everything to the lhs
        id = hash(x)
        op = x.args[2]
        if op == :(==)
            lb = 0.0
            ub = 0.0
        elseif op == :(<=) || op == :(≤)
            lb = -Inf
            ub = 0.0
        else
            @assert op == :(>=) || op == :(≥)
            lb = 0.0
            ub = Inf
        end
        quote
            initNLP($m)
            c = NonlinearConstraint(@processNLExpr($(esc(lhs))), $lb, $ub)
            push!($m.nlpdata.nlconstr, c)
            push!($m.nlpdata.nlconstrlist, c.terms)
            nothing
        end
    else
        # ranged row
        error("Two-sided nonlinear constraints not yet supported")
    end
end
