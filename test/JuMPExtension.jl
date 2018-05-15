module JuMPExtension
# Simple example of JuMP extension used in the tests to check that JuMP works well with extensions
# The main difference between `JuMP.Model` and `JuMPExtension.MyModel` is the fact that in `addvariable` (resp. `addconstraint`),
# `JuMP.Model` applies the modification to its `moibackend` field while
# `JuMPExtension.MyModel` stores the `AbstractVariable` (resp. `AbstractConstraint`) in a list.

using JuMP

mutable struct MyModel <: JuMP.AbstractModel
    variables::Vector{JuMP.AbstractVariable}
    varnames::Vector{String}
    constraints::Vector{JuMP.AbstractConstraint}
    connames::Vector{String}
    objdict::Dict{Symbol, Any}
    function MyModel()
        new(JuMP.AbstractVariable[], String[],   # Variables
            JuMP.AbstractConstraint[], String[], # Constraints
            Dict{Symbol, Any}())
    end
end

JuMP.object_dictionary(m::MyModel) = m.objdict

# Variables
struct MyVariableRef <: JuMP.AbstractVariableRef
    model::MyModel # `model` owning the variable
    idx::Int       # Index in `model.variables`
end
Base.copy(v::MyVariableRef) = v
JuMP.variabletype(::MyModel) = MyVariableRef
function JuMP.addvariable(m::MyModel, v::JuMP.AbstractVariable, name::String="")
    push!(m.variables, v)
    push!(m.varnames, name)
    MyVariableRef(m, length(m.variables))
end

# Constraints
struct MyConstraintRef
    model::MyModel # `model` owning the constraint
    idx::Int       # Index in `model.constraints`
end
function JuMP.addconstraint(m::MyModel, c::JuMP.AbstractConstraint, name::String="")
    push!(m.constraints, c)
    push!(m.connames, name)
    MyConstraintRef(m, length(m.constraints))
end
function JuMP.constraintobject(cref::MyConstraintRef, F::Type, S::Type)
    c = cref.model.constraints[cref.idx]
    # `TypeError` should be thrown is `F` and `S` are not correct
    # This is needed for the tests in `constraints.jl`
    c.func::F
    c.set::S
    c
end

# Names
JuMP.name(vref::MyVariableRef) = vref.model.varnames[vref.idx]
function JuMP.setname(vref::MyVariableRef, name::String)
    cref.model.varnames[vref.idx] = name
end
JuMP.name(cref::MyConstraintRef) = cref.model.connames[cref.idx]
function JuMP.setname(cref::MyConstraintRef, name::String)
    cref.model.connames[cref.idx] = name
end

end
