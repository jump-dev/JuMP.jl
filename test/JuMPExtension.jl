module JuMPExtension
# Simple example of JuMP extension used in the tests to check that JuMP works well with extensions
# The main difference between `JuMP.Model` and `JuMPExtension.MyModel` is the fact that in `addvariable` (resp. `addconstraint`),
# `JuMP.Model` applies the modification to its `moi_backend` field while
# `JuMPExtension.MyModel` stores the `AbstractVariable` (resp. `AbstractConstraint`) in a list.

using MathOptInterface
const MOI = MathOptInterface
using JuMP

mutable struct MyModel <: JuMP.AbstractModel
    nextvaridx::Int                                 # Next variable index is nextvaridx+1
    variables::Dict{Int, JuMP.AbstractVariable}     # Map varidx -> variable
    varnames::Dict{Int, String}                     # Map varidx -> name
    nextconidx::Int                                 # Next constraint index is nextconidx+1
    constraints::Dict{Int, JuMP.AbstractConstraint} # Map conidx -> variable
    connames::Dict{Int, String}                     # Map varidx -> name
    objectivesense::Symbol
    objective_function::JuMP.AbstractJuMPScalar
    obj_dict::Dict{Symbol, Any}                     # Same that JuMP.Model's field `obj_dict`
    function MyModel()
        new(0, Dict{Int, JuMP.AbstractVariable}(),   Dict{Int, String}(), # Variables
            0, Dict{Int, JuMP.AbstractConstraint}(), Dict{Int, String}(), # Constraints
            :Min, zero(JuMP.GenericAffExpr{Float64, MyVariableRef}),
            Dict{Symbol, Any}())
    end
end
if VERSION >= v"0.7-"
    Base.broadcastable(model::MyModel) = Ref(model)
end

JuMP.object_dictionary(model::MyModel) = model.obj_dict

# Variables
struct MyVariableRef <: JuMP.AbstractVariableRef
    model::MyModel # `model` owning the variable
    idx::Int       # Index in `model.variables`
end
Base.copy(v::MyVariableRef) = v
Base.:(==)(v::MyVariableRef, w::MyVariableRef) = v.model === w.model && v.idx == w.idx
if VERSION >= v"0.7-"
    Base.broadcastable(v::MyVariableRef) = Ref(v)
end
JuMP.owner_model(v::MyVariableRef) = v.model
JuMP.isequal_canonical(v::MyVariableRef, w::MyVariableRef) = v == w
JuMP.variable_type(::MyModel) = MyVariableRef
function JuMP.add_variable(m::MyModel, v::JuMP.AbstractVariable, name::String="")
    m.nextvaridx += 1
    vref = MyVariableRef(m, m.nextvaridx)
    m.variables[vref.idx] = v
    JuMP.set_name(vref, name)
    vref
end
function JuMP.delete(model::MyModel, variable_ref::MyVariableRef)
    @assert JuMP.is_valid(model, variable_ref)
    delete!(model.variables, variable_ref.idx)
    delete!(model.varnames, variable_ref.idx)
end
function JuMP.is_valid(model::MyModel, variable_ref::MyVariableRef)
    return (model === variable_ref.model &&
            variable_ref.idx in keys(model.variables))
end
JuMP.num_variables(m::MyModel) = length(m.variables)

JuMP.has_lower_bound(vref::MyVariableRef) = vref.model.variables[vref.idx].info.has_lb
function JuMP.lower_bound(vref::MyVariableRef)
    @assert !JuMP.is_fixed(vref)
    vref.model.variables[vref.idx].info.lower_bound
end
function JuMP.set_lower_bound(vref::MyVariableRef, lower)
    vref.model.variables[vref.idx].info.has_lb = true
    vref.model.variables[vref.idx].info.lower_bound = lower
end
function JuMP.delete_lower_bound(variable_ref::MyVariableRef)
    variable_ref.model.variables[variable_ref.idx].info.has_lb = false
end
JuMP.has_upper_bound(vref::MyVariableRef) = vref.model.variables[vref.idx].info.has_ub
function JuMP.upper_bound(vref::MyVariableRef)
    @assert !JuMP.is_fixed(vref)
    vref.model.variables[vref.idx].info.upper_bound
end
function JuMP.set_upper_bound(vref::MyVariableRef, upper)
    vref.model.variables[vref.idx].info.has_ub = true
    vref.model.variables[vref.idx].info.upper_bound = upper
end
function JuMP.delete_upper_bound(variable_ref::MyVariableRef)
    variable_ref.model.variables[variable_ref.idx].info.has_ub = false
end
JuMP.is_fixed(vref::MyVariableRef) = vref.model.variables[vref.idx].info.has_fix
JuMP.fix_value(vref::MyVariableRef) = vref.model.variables[vref.idx].info.fixed_value
function JuMP.fix(variable_ref::MyVariableRef, value)
    variable_ref.model.variables[variable_ref.idx].info.has_fix = true
    variable_ref.model.variables[variable_ref.idx].info.fixed_value = value
end
function JuMP.unfix(variable_ref::MyVariableRef)
    variable_ref.model.variables[variable_ref.idx].info.has_fix = false
end
JuMP.start_value(vref::MyVariableRef) = vref.model.variables[vref.idx].info.start
function JuMP.set_start_value(vref::MyVariableRef, start)
    vref.model.variables[vref.idx].info.start = start
end
JuMP.is_binary(vref::MyVariableRef) = vref.model.variables[vref.idx].info.binary
function JuMP.set_binary(vref::MyVariableRef)
    @assert !JuMP.is_integer(vref)
    vref.model.variables[vref.idx].info.binary = true
end
function JuMP.unset_binary(vref::MyVariableRef)
    vref.model.variables[vref.idx].info.binary = false
end
JuMP.is_integer(vref::MyVariableRef) = vref.model.variables[vref.idx].info.integer
function JuMP.set_integer(vref::MyVariableRef)
    @assert !JuMP.is_binary(vref)
    vref.model.variables[vref.idx].info.integer = true
end
function JuMP.unset_integer(vref::MyVariableRef)
    vref.model.variables[vref.idx].info.integer = false
end

# Constraints
struct MyConstraintRef
    model::MyModel # `model` owning the constraint
    idx::Int       # Index in `model.constraints`
end
JuMP.constraint_type(::MyModel) = MyConstraintRef
if VERSION >= v"0.7-"
    Base.broadcastable(cref::MyConstraintRef) = Ref(cref)
end
function JuMP.add_constraint(m::MyModel, c::JuMP.AbstractConstraint, name::String="")
    m.nextconidx += 1
    cref = MyConstraintRef(m, m.nextconidx)
    m.constraints[cref.idx] = c
    JuMP.set_name(cref, name)
    cref
end
function JuMP.delete(model::MyModel, constraint_ref::MyConstraintRef)
    @assert JuMP.is_valid(model, constraint_ref)
    delete!(model.constraints, constraint_ref.idx)
    delete!(model.connames, constraint_ref.idx)
end
function JuMP.is_valid(model::MyModel, constraint_ref::MyConstraintRef)
    return (model === constraint_ref.model &&
            constraint_ref.idx in keys(model.constraints))
end
function JuMP.constraint_object(cref::MyConstraintRef)
    return cref.model.constraints[cref.idx]
end

# Objective
function JuMP.set_objective(m::MyModel, sense::Symbol, f::JuMP.AbstractJuMPScalar)
    m.objectivesense = sense
    m.objective_function = f
end
JuMP.objective_sense(m::MyModel) = m.objectivesense
function JuMP.objective_function(m::MyModel, FT::Type)
    # InexactError should be thrown, this is needed in `objective.jl`
    if !(m.objective_function isa FT)
        if VERSION < v"0.7-"
            throw(InexactError())
        else
            throw(InexactError(:objective_function, FT, typeof(m.objective_function)))
        end
    end
    m.objective_function
end

# Names
JuMP.name(vref::MyVariableRef) = vref.model.varnames[vref.idx]
function JuMP.set_name(vref::MyVariableRef, name::String)
    vref.model.varnames[vref.idx] = name
end
JuMP.name(cref::MyConstraintRef) = cref.model.connames[cref.idx]
function JuMP.set_name(cref::MyConstraintRef, name::String)
    cref.model.connames[cref.idx] = name
end

end
