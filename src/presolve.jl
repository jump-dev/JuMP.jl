export presolve
const EQ_TOL = 1e-10
const REL_EQ_TOL = 1e-6

function presolve(m::Model)
    #println(m)

    var_removed = zeros(Bool, m.numCols)
    con_removed = zeros(Bool, length(m.linconstr))
    var_changed = ones(Bool, m.numCols)
    con_changed = ones(Bool, length(m.linconstr))

    any_change  = true
    iterations  = 0
    while any_change
        any_change = false
        iterations += 1

        # PASS 1
        # Check if bounds are at equality - if so, eliminate them
        # from the constraints
        any_change |= presolve_eliminate_variables(m, 
                        var_changed, var_removed,
                        con_changed, con_removed)
        #var_changed = zeros(Bool, m.numCols)

        # PASS 2
        # Turn singleton constraints into variable bounds
        any_change |= presolve_singleton_cons(m, 
                        var_changed, var_removed,
                        con_changed, con_removed)

        # PASS 3
        # Infer upper bounds on variables from constraint RHS
        any_change |= presolve_rhs_tighten(m, 
                        var_changed, var_removed,
                        con_changed, con_removed)
        #con_changed = zeros(Bool, length(m.linconstr))
        
    end

    presolve_finish(m, con_removed, var_removed)
end

#############################################################################
# PASS 1
# Check if bounds are at equality - if so, eliminate them
# from the constraints
function presolve_eliminate_variables(m::Model, var_changed, var_removed,
                                                con_changed, con_removed)
    any_change = false
    
    for j = 1:m.numCols
        (!var_changed[j] || var_removed[j]) && continue

        # Round integer bounds
        if m.colCat[j] == INTEGER
            m.colLower[j] = ceil(m.colLower[j])
            m.colUpper[j] = floor(m.colUpper[j])
        end
        
        # Check equality of bounds
        if abs(m.colUpper[j] - m.colLower[j]) < EQ_TOL
            # Bounds are equal - remove this variable
            var_removed[j] = true
            updates = presolve_remove_col(m, j, m.colLower[j], 
                                                con_removed)
            for i in updates
                con_changed[i] = true
                any_change = true
            end
        end
    end

    return any_change
end

#############################################################################
# PASS 2
# Turn singleton constraints into variable bounds
function presolve_singleton_cons(m::Model, var_changed, var_removed,
                                           con_changed, con_removed)
    any_change = false
    
    for i = 1:length(m.linconstr)
        (!con_changed[i] || con_removed[i]) && continue
        
        # Check to see if this is a singleton constraint
        lhs = m.linconstr[i].terms
        num_nonzero = 0
        nonzero_pos = 0
        for pos = 1:length(lhs.vars)
            if lhs.coeffs[pos] != 0.0
                num_nonzero += 1
                nonzero_pos = pos
            end
            num_nonzero > 1 && break
        end
        # It is!
        if num_nonzero == 1
            j = lhs.vars[nonzero_pos].col
            nz_coeff = lhs.coeffs[nonzero_pos]
            if sense(m.linconstr[i]) == :(==)
                m.colLower[j] = rhs(m.linconstr[i]) / nz_coeff
                m.colUpper[j] = rhs(m.linconstr[i]) / nz_coeff
            elseif sense(m.linconstr[i]) == :(<=)
                m.colUpper[j] = rhs(m.linconstr[i]) / nz_coeff
            else
                m.colLower[j] = rhs(m.linconstr[i]) / nz_coeff
            end
            lhs.coeffs[nonzero_pos] = 0.0
            con_removed[i] = true
            var_changed[j] = true
            any_change = true
        end
    end

    return any_change
end

#############################################################################
# PASS 3
# Infer upper bounds on variables from constraint RHS
function presolve_rhs_tighten(m::Model, var_changed, var_removed,
                                        con_changed, con_removed)
    any_change = false

    for i = 1:length(m.linconstr)
        (!con_changed[i] || con_removed[i]) && continue

        # We'll be lazy right now and just deal with constraints
        # where variables are >= 0
        lhs = m.linconstr[i].terms
        for pos = 1:length(lhs.vars)
            col = lhs.vars[pos].col
            if m.colLower[col] < 0.0
                continue
            end
        end
        # All variables are >= 0, so we can tighten
        if sense(m.linconstr[i]) == :(<=) ||
           sense(m.linconstr[i]) == :(==)
            for pos = 1:length(lhs.vars)
                col = lhs.vars[pos].col
                coeff = lhs.coeffs[pos]
                if coeff != 0.0 &&
                   rhs(m.linconstr[i])/coeff < m.colUpper[col]
                    m.colUpper[col] = rhs(m.linconstr[i])/coeff

                    var_changed[col] = true
                    any_change = true
                end
            end 
        end
    end

    return any_change
end


#############################################################################

function presolve_remove_col(m::Model, j::Int, new_value::Float64,
                             con_removed::Vector{Bool})
    con_changed = Int[]

    for i = 1:length(m.linconstr)
        # Don't consider this constraint at all if its already dead
        con_removed[i] && continue
        # Walk through terms looking for our variable
        lhs = m.linconstr[i].terms
        any_nonzero = false
        for pos = 1:length(lhs.vars)
            coeff = lhs.coeffs[pos]
            # If nonzero, consider (otherwise already removed)
            if coeff != 0.0
                if lhs.vars[pos].col == j
                    # Found our variable, remove it from the constraint
                    # and update the constraint bounds
                    m.linconstr[i].lb -= coeff * new_value
                    m.linconstr[i].ub -= coeff * new_value
                    lhs.coeffs[pos] = 0.0
                else
                    # Found a variable we won't remove that still
                    # exists in the constraint
                    any_nonzero = true
                end
            end
        end
        # If every zero, this constraint is dead
        if !any_nonzero
            con_removed[i] = true
            push!(con_changed, i)
        end
    end
    # Finally check objective
    for pos = 1:length(m.obj.aff.vars)
        if m.obj.aff.vars[pos].col == j
            m.obj.aff.constant += m.obj.aff.coeffs[pos] * new_value
            m.obj.aff.coeffs[pos] = 0.0
        end
    end

    return con_changed
end

#############################################################################

function presolve_finish(m::Model, con_removed::Vector{Bool},
                                   var_removed::Vector{Bool})
    # Reassign columns
    new_numCols = 0
    new_column_index = zeros(m.numCols)
    old_column_index = zeros(m.numCols)
    for j = 1:m.numCols
        if !var_removed[j]
            new_numCols += 1
            new_column_index[j] = new_numCols
            old_column_index[new_numCols] = j
        end
    end

    # Create new model
    newm = Model(solver=m.solver)
    newm.objSense = m.objSense
    newm.numCols  = new_numCols
    newm.colNames = [m.colNames[old_column_index[j]] for j = 1:new_numCols]
    newm.colLower = [m.colLower[old_column_index[j]] for j = 1:new_numCols]
    newm.colUpper = [m.colUpper[old_column_index[j]] for j = 1:new_numCols]
    newm.colCat   = [m.colCat[old_column_index[j]]   for j = 1:new_numCols]

    # Move objective over
    newobj = QuadExpr()
    for pos = 1:length(m.obj.aff.vars)
        if m.obj.aff.coeffs[pos] != 0.0
            new_col = new_column_index[m.obj.aff.vars[pos].col]
            push!(newobj.aff, m.obj.aff.coeffs[pos], Variable(newm, new_col))
        end
    end
    newobj.aff.constant = m.obj.aff.constant
    newm.obj = newobj

    # Move constraints over
    for i = 1:length(m.linconstr)
        # Don't move this constraint at all if its already dead
        con_removed[i] && continue

        con = m.linconstr[i]
        new_terms = AffExpr()
        for pos = 1:length(con.terms.vars)
            new_coeff = con.terms.coeffs[pos]
            if new_coeff != 0.0
                new_col = new_column_index[con.terms.vars[pos].col]
                push!(new_terms, new_coeff, Variable(newm, new_col))
            end
        end

        push!(newm.linconstr,
                LinearConstraint(new_terms, con.lb, con.ub))
    end
    new_numCons = length(newm.linconstr)

    println("JUMP PRESOLVE")
    println("New variables:   $new_numCols (removed $(m.numCols-new_numCols))")
    println("New constraints: $new_numCons (removed $(length(m.linconstr)-new_numCons))")

    println(newm)
    solve(newm)
    println("Objective value ", getObjectiveValue(newm))
end