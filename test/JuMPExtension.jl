module JuMPExtension
# Simple example of JuMP extension used in the tests to check that JuMP works well with extensions
# The main difference between `JuMP.Model` and `JuMPExtension.MyModel` is the fact that in `add_variable` (resp. `add_constraint`),
# `JuMP.Model` applies the modification to its `moi_backend` field while
# `JuMPExtension.MyModel` stores the `AbstractVariable` (resp. `AbstractConstraint`) in a list.

using MathOptInterface
const MOI = MathOptInterface
using JuMP

struct ConstraintIndex
    value::Int # Index in `model.constraints`
end

mutable struct MyModel <: JuMP.AbstractModel
    nextvaridx::Int                                 # Next variable index is nextvaridx+1
    variables::Dict{Int, JuMP.ScalarVariable}       # Map varidx -> variable
    var_to_name::Dict{Int, String}                  # Map varidx -> name
    name_to_var::Dict{String, Int}                  # Map varidx -> name
    nextconidx::Int                                 # Next constraint index is nextconidx+1
    constraints::Dict{ConstraintIndex,
                      JuMP.AbstractConstraint}      # Map conidx -> variable
    con_to_name::Dict{ConstraintIndex, String}         # Map conidx -> name
    objectivesense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    obj_dict::Dict{Symbol, Any}                     # Same that JuMP.Model's field `obj_dict`
    function MyModel()
        new(0, Dict{Int, JuMP.AbstractVariable}(),
            Dict{Int, String}(), Dict{String, Int}(),            # Variables
            0, Dict{ConstraintIndex, JuMP.AbstractConstraint}(),
            Dict{ConstraintIndex, String}(),                     # Constraints
            MOI.FEASIBILITY_SENSE,
            zero(JuMP.GenericAffExpr{Float64, MyVariableRef}),
            Dict{Symbol, Any}())
    end
end
Base.broadcastable(model::MyModel) = Ref(model)

JuMP.object_dictionary(model::MyModel) = model.obj_dict

# Variables
struct MyVariableRef <: JuMP.AbstractVariableRef
    model::MyModel # `model` owning the variable
    idx::Int       # Index in `model.variables`
end
Base.copy(v::MyVariableRef) = v
Base.:(==)(v::MyVariableRef, w::MyVariableRef) = v.model === w.model && v.idx == w.idx
Base.broadcastable(v::MyVariableRef) = Ref(v)
JuMP.isequal_canonical(v::MyVariableRef, w::MyVariableRef) = v == w
JuMP.variable_type(::MyModel) = MyVariableRef
function JuMP.add_variable(m::MyModel, v::JuMP.AbstractVariable, name::String="")
    m.nextvaridx += 1
    vref = MyVariableRef(m, m.nextvaridx)
    m.variables[vref.idx] = v
    JuMP.set_name(vref, name)
    vref
end
function JuMP.delete(model::MyModel, vref::MyVariableRef)
    @assert JuMP.is_valid(model, vref)
    delete!(model.variables, vref.idx)
    delete!(model.var_to_name, vref.idx)
end
function JuMP.is_valid(model::MyModel, vref::MyVariableRef)
    return (model === vref.model &&
            vref.idx in keys(model.variables))
end
JuMP.num_variables(m::MyModel) = length(m.variables)

# Internal functions
variable_info(vref::MyVariableRef) = vref.model.variables[vref.idx].info
function update_variable_info(vref::MyVariableRef, info::JuMP.VariableInfo)
    vref.model.variables[vref.idx] = JuMP.ScalarVariable(info)
end

JuMP.has_lower_bound(vref::MyVariableRef) = variable_info(vref).has_lb
function JuMP.lower_bound(vref::MyVariableRef)
    @assert !JuMP.is_fixed(vref)
    return variable_info(vref).lower_bound
end
function JuMP.set_lower_bound(vref::MyVariableRef, lower)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(true, lower,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
function JuMP.delete_lower_bound(vref::MyVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(false, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
JuMP.has_upper_bound(vref::MyVariableRef) = variable_info(vref).has_ub
function JuMP.upper_bound(vref::MyVariableRef)
    @assert !JuMP.is_fixed(vref)
    return variable_info(vref).upper_bound
end
function JuMP.set_upper_bound(vref::MyVariableRef, upper)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           true, upper,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
function JuMP.delete_upper_bound(vref::MyVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           false, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
JuMP.is_fixed(vref::MyVariableRef) = variable_info(vref).has_fix
JuMP.fix_value(vref::MyVariableRef) = variable_info(vref).fixed_value
function JuMP.fix(vref::MyVariableRef, value; force::Bool = false)
    info = variable_info(vref)
    if !force && (info.has_lb || info.has_ub)
        error("Unable to fix $(vref) to $(value) because it has existing bounds.")
    end
    update_variable_info(vref, JuMP.VariableInfo(
        false, info.lower_bound, false, info.upper_bound, true, value,
        info.has_start, info.start, info.binary, info.integer)
    )
    return
end
function JuMP.unfix(vref::MyVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           false, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
JuMP.start_value(vref::MyVariableRef) = variable_info(vref).start
function JuMP.set_start_value(vref::MyVariableRef, start)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           true, start,
                                           info.binary, info.integer))
end
JuMP.is_binary(vref::MyVariableRef) = variable_info(vref).binary
function JuMP.set_binary(vref::MyVariableRef)
    @assert !JuMP.is_integer(vref)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           true, info.integer))
end
function JuMP.unset_binary(vref::MyVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           false, info.integer))
end
JuMP.is_integer(vref::MyVariableRef) = variable_info(vref).integer
function JuMP.set_integer(vref::MyVariableRef)
    @assert !JuMP.is_binary(vref)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, true))
end
function JuMP.unset_integer(vref::MyVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, false))
end

# Constraints
const MyConstraintRef = JuMP.ConstraintRef{MyModel, ConstraintIndex}
JuMP.constraint_type(::MyModel) = MyConstraintRef
function JuMP.add_constraint(model::MyModel, c::JuMP.AbstractConstraint,
                             name::String="")
    model.nextconidx += 1
    index = ConstraintIndex(model.nextconidx)
    cref = JuMP.ConstraintRef(model, index, JuMP.shape(c))
    model.constraints[index] = c
    JuMP.set_name(cref, name)
    return cref
end
function JuMP.delete(model::MyModel, constraint_ref::MyConstraintRef)
    @assert JuMP.is_valid(model, constraint_ref)
    delete!(model.constraints, constraint_ref.index)
    delete!(model.con_to_name, constraint_ref.index)
end
function JuMP.is_valid(model::MyModel, constraint_ref::MyConstraintRef)
    return (model === constraint_ref.model &&
            constraint_ref.index in keys(model.constraints))
end
function JuMP.constraint_object(cref::MyConstraintRef)
    return cref.model.constraints[cref.index]
end

# Objective
function JuMP.set_objective(m::MyModel, sense::MOI.OptimizationSense,
                            f::JuMP.AbstractJuMPScalar)
    m.objectivesense = sense
    m.objective_function = f
end
JuMP.objective_sense(model::MyModel) = model.objectivesense
function JuMP.set_objective_sense(model::MyModel, sense)
    model.objectivesense = sense
end
JuMP.objective_function_type(model::MyModel) = typeof(model.objective_function)
JuMP.objective_function(model::MyModel) = model.objective_function
function JuMP.objective_function(model::MyModel, FT::Type)
    # InexactError should be thrown, this is needed in `objective.jl`
    if !(model.objective_function isa FT)
        throw(InexactError(:objective_function, FT,
                           typeof(model.objective_function)))
    end
    model.objective_function
end

# Names
JuMP.name(vref::MyVariableRef) = vref.model.var_to_name[vref.idx]
function JuMP.set_name(vref::MyVariableRef, name::String)
    vref.model.var_to_name[vref.idx] = name
    vref.model.name_to_var[name] = vref.idx
end
function JuMP.variable_by_name(model::MyModel, name::String)
    return MyVariableRef(model, model.name_to_var[name])
end
JuMP.name(cref::MyConstraintRef) = cref.model.con_to_name[cref.index]
function JuMP.set_name(cref::MyConstraintRef, name::String)
    cref.model.con_to_name[cref.index] = name
end

end
