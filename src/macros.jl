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

addToExpression(aff::AffExpr,x::Variable,c::Number) = addToExpression(aff,c,x)
addToExpression(aff::AffExpr,c::Number,x::Variable) = push!(aff,convert(Float64,c),x)

addToExpression(aff::AffExpr,c::Number,x::Number) = (aff.constant += c*x)

addToExpression(aff::AffExpr,c::Variable,x::Variable) = QuadExpr([c],[x],[1.0],aff)

addToExpression(aff::AffExpr,x::AffExpr,c::Number) = addToExpression(aff,c,x)
function addToExpression(aff::AffExpr,c::Number,x::AffExpr)
    append!(aff.vars, x.vars)
    append!(aff.coeffs, c*x.coeffs)
    aff.constant += c*x.constant
end

addToExpression(aff::AffExpr,x::AffExpr,c::Variable) = addToExpression(aff,c,x)
addToExpression(aff::AffExpr,c::Variable,x::AffExpr) = QuadExpr(fill(c,length(x.vars)),
                                                                x.vars,
                                                                x.coeffs,
                                                                aff)

addToExpression(aff::AffExpr, c::Number, x::QuadExpr) = QuadExpr(copy(x.qvars1),
                                                                 copy(x.qvars2), 
                                                                 c*x.qcoeffs, 
                                                                 AffExpr(vcat(aff.vars,x.aff.vars),
                                                                         vcat(aff.coeffs,c*x.aff.coeffs),
                                                                         aff.constant+c*x.aff.constant))

addToExpression(quad::QuadExpr,x::Variable,c::Float64) = addToExpression(quad,c,x)
addToExpression(quad::QuadExpr,c::Number,x::Variable) = push!(quad.aff,convert(Float64,c),x)

addToExpression(quad::QuadExpr,c::Number,x::Number) = (quad.aff.constant += c*x)

function addToExpression(quad::QuadExpr,c::Variable,x::Variable)
    push!(quad.qvars1, x)
    push!(quad.qvars2, c)
end

addToExpression(quad::QuadExpr,x::AffExpr,c::Number) = addToExpression(quad::QuadExpr,c::Number,x::AffExpr)
function addToExpression(quad::QuadExpr,c::Number,x::AffExpr)
    append!(quad.aff.vars, x.vars)
    append!(quad.aff.coeffs, c*x.coeffs)
    quad.aff.constant += c*x.constant
end

addToExpression(quad::QuadExpr,x::AffExpr,c::Variable) = addToExpression(quad,c,x)
function addToExpression(quad::QuadExpr,c::Variable,x::AffExpr)
    append!(quad.qvars1, fill(c,length(x.vars)))
    append!(quad.qvars2, x.vars)
    append!(quad.qcoeffs, x.coeffs)
    addToExpression(quad.aff,x.constant,c)
end

addToExpression(quad::QuadExpr,x::QuadExpr,c::Number) = addToExpression(quad::QuadExpr,c::Number,x::QuadExpr)
function addToExpression(quad::QuadExpr,c::Number,x::QuadExpr)
    append!(quad.qvars1,x.qvars1)
    append!(quad.qvars2,x.qvars2)
    append!(quad.qcoeffs,c*x.qcoeffs)
    addToExpression(quad,c,x.aff)
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
            sizehint($aff.aff.vars,length($aff.aff.vars)+$len);
            sizehint($aff.aff.coeffs,length($aff.aff.coeffs)+$len))
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
    # Two formats:
    # - @addConstraint(m, a*x <= 5)
    # - @addConstraint(m, myref[a=1:5], a*x <= 5)
    length(extra) > 1 && error("in @addConstraint: too many arguments.")
    # Canonicalize the arguments
    c = length(extra) == 1 ? x        : nothing
    x = length(extra) == 1 ? extra[1] : x

    (x.head != :comparison) &&
        error("in @addConstraint ($(string(x))): expected comparison operator (<=, >=, or ==).")

    # Strategy: build up the code for non-macro addconstraint, and if needed
    # we will wrap in loops to assign to the ConstraintRefs

    # Build up machinery for looping over index sets, if needed
    refcall = gensym()  # Default
    if isa(c,Symbol)
        # Just creating a simple ConstraintRefs (not indexed)
        refcall = esc(c)
    elseif isexpr(c,:ref)
        # Creating an indexed set of ConstraintRefs
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
    elseif c != nothing
        # Something in there, but we don't know what
        error("in @addConstraint ($(string(x))): expected $(string(c)) to be of form constr[...]")
    end

    # Build the constraint
    if length(x.args) == 3
        # Simple comparison - move everything to the LHS
        if !((x.args[2] == :<=) || (x.args[2] == :≤) ||
             (x.args[2] == :>=) || (x.args[2] == :≥) ||
             (x.args[2] == :(==)))
            error("in @addConstraint ($(string(x))): expected comparison operator (<=, >=, or ==).")
        end
        lhs = :($(x.args[1]) - $(x.args[3])) 
        code = quote
            q = QuadExpr()
            $(parseExpr(lhs, :q, 1.0))
            islinear = isempty(q.qvars1)
            crefflag && !islinear && error("Three argument form form of @addConstraint does not currently support quadratic constraints")
            $(refcall) = (islinear ? addConstraint($m, $(x.args[2])(q.aff,0)) :
                                     addConstraint($m, $(x.args[2])(q,0)) )
        end
    elseif length(x.args) == 5
        # Ranged row
        if (x.args[2] != :<= && x.args[2] != :≤) || (x.args[4] != :<= && x.args[4] != :≤)
            error("in @addConstraint ($(string(x))): only ranged rows of the form lb <= expr <= ub are supported.")
        end
        lb = x.args[1]
        ub = x.args[5]
        code = quote
            q = QuadExpr()
            if !isa($(esc(lb)),Number)
                error(string("in @addConstraint (",$(string(x)),"): expected ",$(string(lb))," to be a number."))
            elseif !isa($(esc(ub)),Number)
                error(string("in @addConstraint (",$(string(x)),"): expected ",$(string(ub))," to be a number."))
            end
            $(parseExpr(x.args[3],:q,1.0))
            islinear = isempty(q.qvars1)
            crefflag && !islinear && error("Three argument form form of @addConstraint does not currently support quadratic constraints")
            $(refcall) = (islinear ? addConstraint($m, LinearConstraint(q.aff,$(esc(lb))-q.aff.constant,$(esc(ub))-q.aff.constant)) :
                                     addConstraint($m,   QuadConstraint(q,    $(esc(lb))-q.aff.constant,$(esc(ub))-q.aff.constant)) )
        end
    else
        # Unknown
        error("in @addConstraint ($(string(x))): constraints must be in one of the following forms:\n" *
              "       expr1 <= expr2\n" * "       expr1 >= expr2\n" *
              "       expr1 == expr2\n" * "       lb <= expr <= ub")
    end

    # Combine indexing (if needed) with constraint code
    if isa(c,Symbol) || c == nothing
        # No indexing
        return quote crefflag = false; $code; end
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
            crefflag = true
            $code
            nothing
        end
    end
end

macro addConstraints(m, x)
    x.head == :block || error("Invalid syntax for @addConstraints")
    @assert x.args[1].head == :line
    code = quote end
    for it in x.args
        if it.head == :line
            # do nothing
        elseif it.head == :comparison # regular constraint
            mac = Expr(:macrocall,symbol("@addConstraint"), esc(m), esc(it))
            code = quote
                    $code
                    $mac
                    end
        elseif it.head == :tuple # constraint ref
            mac = Expr(:macrocall,symbol("@addConstraint"), esc(m), esc(it.args[1]), esc(it.args[2]))
            code = quote
                    $code
                    $mac
                    end
        end
    end
    return quote 
            $code 
            nothing
            end
end

macro setObjective(m, args...)
    if length(args) != 2
        # Either just an objective sene, or just an expression.
        error("in @setObjective: needs two arguments: objective sense (Max or Min) and expression.")
    end
    sense, x = args
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    quote
        q = QuadExpr()
        $(parseExpr(x, :q, 1.0))
        setObjective($(esc(m)), $(esc(sense)), q)
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
                error("in @defVar ($var): use the form lb <= $var <= ub instead of ub >= $var >= lb.")
            @assert length(x.args) == 3
            lb = esc(x.args[3])
            ub = Inf
        elseif x.args[2] == :<= || x.args[2] == :≤
            if length(x.args) == 5
                # lb <= x <= u
                var = x.args[3]
                (x.args[4] != :<= && x.args[4] != :≤) &&
                    error("in @defVar ($var): expected <= operator after variable name.")
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
            error("in @defVar ($(string(x))): use the form lb <= ... <= ub.")
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
