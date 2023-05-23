#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module JuMPExtension

# Simple example of JuMP extension used in the tests to check that JuMP works
# well with extensions.
#
# The main difference between `JuMP.Model` and `JuMPExtension.MyModel` is the
# fact that in `add_variable` (resp. `add_constraint`), `JuMP.Model` applies the
# modification to its `moi_backend` field while `JuMPExtension.MyModel` stores
# the `AbstractVariable` (resp. `AbstractConstraint`) in a list.

import JuMP
import Test
import OrderedCollections

struct ConstraintIndex
    value::Int
end

mutable struct MyModel <: JuMP.AbstractModel
    next_var_index::Int
    variables::OrderedCollections.OrderedDict{Int,JuMP.ScalarVariable}      # Map varidx -> variable
    var_to_name::Dict{Int,String}
    name_to_var::Union{Dict{String,Int},Nothing}
    next_con_index::Int
    constraints::Dict{ConstraintIndex,JuMP.AbstractConstraint}
    con_to_name::Dict{ConstraintIndex,String}
    name_to_con::Union{Dict{String,ConstraintIndex},Nothing}
    objective_sense::JuMP.MOI.OptimizationSense
    objective_function::Union{
        JuMP.AbstractJuMPScalar,
        Vector{<:JuMP.AbstractJuMPScalar},
    }
    obj_dict::Dict{Symbol,Any}
    function MyModel()
        return new(
            0,
            OrderedCollections.OrderedDict{Int,JuMP.ScalarVariable}(),
            Dict{Int,String}(),
            nothing,
            0,
            Dict{ConstraintIndex,JuMP.AbstractConstraint}(),
            Dict{ConstraintIndex,String}(),
            nothing,
            JuMP.MOI.FEASIBILITY_SENSE,
            zero(JuMP.GenericAffExpr{Float64,MyVariableRef}),
            Dict{Symbol,Any}(),
        )
    end
end

Base.broadcastable(model::MyModel) = Ref(model)

JuMP.object_dictionary(model::MyModel) = model.obj_dict

struct MyVariableRef <: JuMP.AbstractVariableRef
    model::MyModel
    idx::Int
end

JuMP.variable_ref_type(::Union{MyModel,Type{MyModel}}) = MyVariableRef
Base.copy(v::MyVariableRef) = v

function Base.:(==)(v::MyVariableRef, w::MyVariableRef)
    return v.model === w.model && v.idx == w.idx
end

Base.broadcastable(v::MyVariableRef) = Ref(v)

JuMP.isequal_canonical(v::MyVariableRef, w::MyVariableRef) = v == w

function JuMP.add_variable(
    m::MyModel,
    v::JuMP.AbstractVariable,
    name::String = "",
)
    m.next_var_index += 1
    vref = MyVariableRef(m, m.next_var_index)
    m.variables[vref.idx] = v
    JuMP.set_name(vref, name)
    return vref
end

function JuMP.add_variable(
    model::MyModel,
    variable::JuMP.VariableConstrainedOnCreation,
    name::String,
)
    var_ref = JuMP.add_variable(model, variable.scalar_variable, name)
    JuMP.add_constraint(model, JuMP.ScalarConstraint(var_ref, variable.set))
    return var_ref
end

function JuMP.add_variable(
    model::MyModel,
    variable::JuMP.VariablesConstrainedOnCreation,
    names,
)
    var_refs =
        JuMP.add_variable.(
            model,
            variable.scalar_variables,
            JuMP.vectorize(names, variable.shape),
        )
    if !isa(variable.set, JuMP.MOI.Reals)
        constraint = JuMP.VectorConstraint(var_refs, variable.set)
        JuMP.add_constraint(model, constraint)
    end
    return JuMP.reshape_vector(var_refs, variable.shape)
end

function JuMP.delete(model::MyModel, vref::MyVariableRef)
    @assert JuMP.is_valid(model, vref)
    delete!(model.variables, vref.idx)
    delete!(model.var_to_name, vref.idx)
    return
end

function JuMP.delete(model::MyModel, vref::Vector{MyVariableRef})
    JuMP.delete.(model, vref)
    return
end

function JuMP.is_valid(model::MyModel, vref::MyVariableRef)
    return model === vref.model && vref.idx in keys(model.variables)
end

function JuMP.all_variables(model::MyModel)
    return [MyVariableRef(model, idx) for idx in keys(model.variables)]
end

JuMP.num_variables(m::MyModel) = length(m.variables)

_variable_info(vref::MyVariableRef) = vref.model.variables[vref.idx].info

function _update_variable_info(vref::MyVariableRef, info::JuMP.VariableInfo)
    vref.model.variables[vref.idx] = JuMP.ScalarVariable(info)
    return
end

JuMP.has_lower_bound(vref::MyVariableRef) = _variable_info(vref).has_lb

function JuMP.lower_bound(vref::MyVariableRef)::Float64
    @assert !JuMP.is_fixed(vref)
    return _variable_info(vref).lower_bound
end

function JuMP.set_lower_bound(vref::MyVariableRef, lower)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            true,
            lower,
            info.has_ub,
            info.upper_bound,
            info.has_fix,
            info.fixed_value,
            info.has_start,
            info.start,
            info.binary,
            info.integer,
        ),
    )
    return
end

function JuMP.delete_lower_bound(vref::MyVariableRef)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            false,
            info.lower_bound,
            info.has_ub,
            info.upper_bound,
            info.has_fix,
            info.fixed_value,
            info.has_start,
            info.start,
            info.binary,
            info.integer,
        ),
    )
    return
end

JuMP.has_upper_bound(vref::MyVariableRef) = _variable_info(vref).has_ub

function JuMP.upper_bound(vref::MyVariableRef)::Float64
    @assert !JuMP.is_fixed(vref)
    return _variable_info(vref).upper_bound
end

function JuMP.set_upper_bound(vref::MyVariableRef, upper)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            info.has_lb,
            info.lower_bound,
            true,
            upper,
            info.has_fix,
            info.fixed_value,
            info.has_start,
            info.start,
            info.binary,
            info.integer,
        ),
    )
    return
end

function JuMP.delete_upper_bound(vref::MyVariableRef)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            info.has_lb,
            info.lower_bound,
            false,
            info.upper_bound,
            info.has_fix,
            info.fixed_value,
            info.has_start,
            info.start,
            info.binary,
            info.integer,
        ),
    )
    return
end

JuMP.is_fixed(vref::MyVariableRef) = _variable_info(vref).has_fix

function JuMP.fix_value(vref::MyVariableRef)::Float64
    return _variable_info(vref).fixed_value
end

function JuMP.fix(vref::MyVariableRef, value; force::Bool = false)
    info = _variable_info(vref)
    if !force && (info.has_lb || info.has_ub)
        error(
            "Unable to fix $(vref) to $(value) because it has existing bounds.",
        )
    end
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            false,
            info.lower_bound,
            false,
            info.upper_bound,
            true,
            value,
            info.has_start,
            info.start,
            info.binary,
            info.integer,
        ),
    )
    return
end

function JuMP.unfix(vref::MyVariableRef)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            info.has_lb,
            info.lower_bound,
            info.has_ub,
            info.upper_bound,
            false,
            info.fixed_value,
            info.has_start,
            info.start,
            info.binary,
            info.integer,
        ),
    )
    return
end

function JuMP.start_value(vref::MyVariableRef)::Union{Nothing,Float64}
    info = _variable_info(vref)
    return info.has_start ? info.start : nothing
end

function JuMP.set_start_value(vref::MyVariableRef, start)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            info.has_lb,
            info.lower_bound,
            info.has_ub,
            info.upper_bound,
            info.has_fix,
            info.fixed_value,
            true,
            start,
            info.binary,
            info.integer,
        ),
    )
    return
end

JuMP.is_binary(vref::MyVariableRef) = _variable_info(vref).binary

function JuMP.set_binary(vref::MyVariableRef)
    @assert !JuMP.is_integer(vref)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            info.has_lb,
            info.lower_bound,
            info.has_ub,
            info.upper_bound,
            info.has_fix,
            info.fixed_value,
            info.has_start,
            info.start,
            true,
            info.integer,
        ),
    )
    return
end

function JuMP.unset_binary(vref::MyVariableRef)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            info.has_lb,
            info.lower_bound,
            info.has_ub,
            info.upper_bound,
            info.has_fix,
            info.fixed_value,
            info.has_start,
            info.start,
            false,
            info.integer,
        ),
    )
    return
end

JuMP.is_integer(vref::MyVariableRef) = _variable_info(vref).integer

function JuMP.set_integer(vref::MyVariableRef)
    @assert !JuMP.is_binary(vref)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            info.has_lb,
            info.lower_bound,
            info.has_ub,
            info.upper_bound,
            info.has_fix,
            info.fixed_value,
            info.has_start,
            info.start,
            info.binary,
            true,
        ),
    )
    return
end

function JuMP.unset_integer(vref::MyVariableRef)
    info = _variable_info(vref)
    _update_variable_info(
        vref,
        JuMP.VariableInfo(
            info.has_lb,
            info.lower_bound,
            info.has_ub,
            info.upper_bound,
            info.has_fix,
            info.fixed_value,
            info.has_start,
            info.start,
            info.binary,
            false,
        ),
    )
    return
end

const MyConstraintRef = JuMP.ConstraintRef{MyModel,ConstraintIndex}

function JuMP.add_constraint(
    model::MyModel,
    c::JuMP.AbstractConstraint,
    name::String = "",
)
    model.next_con_index += 1
    index = ConstraintIndex(model.next_con_index)
    cref = JuMP.ConstraintRef(model, index, JuMP.shape(c))
    model.constraints[index] = c
    JuMP.set_name(cref, name)
    return cref
end

function JuMP.delete(model::MyModel, constraint_ref::MyConstraintRef)
    @assert JuMP.is_valid(model, constraint_ref)
    delete!(model.constraints, constraint_ref.index)
    delete!(model.con_to_name, constraint_ref.index)
    return
end

function JuMP.delete(model::MyModel, constraint_ref::Vector{<:MyConstraintRef})
    JuMP.delete.(model, constraint_ref)
    return
end

function JuMP.is_valid(model::MyModel, constraint_ref::MyConstraintRef)
    return model === constraint_ref.model &&
           haskey(model.constraints, constraint_ref.index)
end

function JuMP.all_constraints(
    model::MyModel,
    ::Type{F},
    ::Type{S},
) where {F,S<:JuMP.MOI.AbstractSet}
    constraints = JuMP.ConstraintRef[]
    for (index, c) in model.constraints
        if JuMP.jump_function(c) isa F && JuMP.moi_set(c) isa S
            push!(constraints, JuMP.ConstraintRef(model, index, JuMP.shape(c)))
        end
    end
    return constraints
end

function JuMP.constraint_object(cref::MyConstraintRef)
    return cref.model.constraints[cref.index]
end

function JuMP.num_constraints(
    model::MyModel,
    F::Type{<:JuMP.AbstractJuMPScalar},
    S::Type{<:JuMP.MOI.AbstractSet},
)
    return count(values(model.constraints)) do con
        return con isa JuMP.ScalarConstraint{F,S}
    end
end

function JuMP.num_constraints(
    model::MyModel,
    ::Type{<:Vector{F}},
    S::Type{<:JuMP.MOI.AbstractSet},
) where {F<:JuMP.AbstractJuMPScalar}
    return count(values(model.constraints)) do con
        return con isa JuMP.VectorConstraint{F,S}
    end
end

function JuMP.num_constraints(
    model::MyModel;
    count_variable_in_set_constraints::Bool,
)
    return count(values(model.constraints)) do con
        is_variable_set = JuMP.moi_set(con) isa MyVariableRef
        return count_variable_in_set_constraints && is_variable_set
    end
end

function JuMP.set_objective_function(
    m::MyModel,
    f::Union{JuMP.AbstractJuMPScalar,Vector{<:JuMP.AbstractJuMPScalar}},
)
    m.objective_function = f
    return
end

function JuMP.set_objective_function(m::MyModel, f::Real)
    m.objective_function = JuMP.GenericAffExpr{Float64,MyVariableRef}(f)
    return
end

JuMP.objective_sense(model::MyModel) = model.objective_sense

function JuMP.set_objective_sense(model::MyModel, sense)
    model.objective_sense = sense
    return
end

JuMP.objective_function_type(model::MyModel) = typeof(model.objective_function)

JuMP.objective_function(model::MyModel) = model.objective_function

function JuMP.objective_function(model::MyModel, FT::Type)
    F = typeof(model.objective_function)
    if F != FT
        throw(InexactError(:objective_function, FT, F))
    end
    return model.objective_function::FT
end

JuMP.name(vref::MyVariableRef) = vref.model.var_to_name[vref.idx]

function JuMP.set_name(vref::MyVariableRef, name::String)
    vref.model.var_to_name[vref.idx] = name
    vref.model.name_to_var = nothing
    return
end

function JuMP.variable_by_name(model::MyModel, name::String)
    if model.name_to_var === nothing
        # Inspired from MOI/src/Utilities/model.jl
        model.name_to_var = Dict{String,Int}()
        for (var, var_name) in model.var_to_name
            if haskey(model.name_to_var, var_name)
                # -1 is a special value that means this string does not map to
                # a unique variable name.
                model.name_to_var[var_name] = -1
            else
                model.name_to_var[var_name] = var
            end
        end
    end
    index = get(model.name_to_var, name, nothing)
    if index === nothing
        return nothing
    elseif index == -1
        error("Multiple variables have the name $name.")
    end
    return MyVariableRef(model, index)
end

JuMP.name(cref::MyConstraintRef) = cref.model.con_to_name[cref.index]

function JuMP.set_name(cref::MyConstraintRef, name::String)
    cref.model.con_to_name[cref.index] = name
    cref.model.name_to_con = nothing
    return
end

function JuMP.constraint_by_name(model::MyModel, name::String)
    if model.name_to_con === nothing
        # Inspired from MOI/src/Utilities/model.jl
        model.name_to_con = Dict{String,ConstraintIndex}()
        for (con, con_name) in model.con_to_name
            if haskey(model.name_to_con, con_name)
                # -1 is a special value that means this string does not map to
                # a unique constraint name.
                model.name_to_con[con_name] = ConstraintIndex(-1)
            else
                model.name_to_con[con_name] = con
            end
        end
    end
    index = get(model.name_to_con, name, nothing)
    if index === nothing
        return nothing
    elseif index.value == -1
        error("Multiple constraints have the name $name.")
    end
    # We have no information on whether this is a vector constraint
    # or a scalar constraint
    return JuMP.ConstraintRef(model, index, JuMP.ScalarShape())
end

JuMP.show_backend_summary(io::IO, model::MyModel) = nothing

function JuMP.show_objective_function_summary(io::IO, model::MyModel)
    T = JuMP.objective_function_type(model)
    println(io, "Objective function type: ", T)
    return
end

function JuMP.objective_function_string(print_mode, model::MyModel)
    return JuMP.function_string(print_mode, JuMP.objective_function(model))
end

_plural(n) = isone(n) ? "" : "s"

function JuMP.show_constraints_summary(io::IO, model::MyModel)
    n = length(model.constraints)
    print(io, "Constraint", _plural(n), ": ", n)
    return
end

function JuMP.constraints_string(print_mode, model::MyModel)
    return String[
        JuMP.constraint_string(print_mode, c) for
        (_, c) in sort(collect(model.constraints); by = c -> c.first.value)
    ]
end

###
### Tests
###

# Test printing of models of type `ModelType` for which the model is stored in
# its JuMP form, e.g., as `AbstractVariable`s and `AbstractConstraint`s.
# This is used by `JuMPExtension` but can also be used by external packages such
# as `StructJuMP`, see https://github.com/jump-dev/JuMP.jl/issues/1711

function test_model_extension_printing()
    repl(s) = JuMP._math_symbol(MIME("text/plain"), s)

    model_1 = MyModel()
    JuMP.@variable(model_1, a >= 1)
    JuMP.@variable(model_1, b <= 1)
    JuMP.@variable(model_1, -1 <= c <= 1)
    JuMP.@variable(model_1, a1 >= 1, Int)
    JuMP.@variable(model_1, b1 <= 1, Int)
    JuMP.@variable(model_1, -1 <= c1 <= 1, Int)
    JuMP.@variable(model_1, x, Bin)
    JuMP.@variable(model_1, y)
    JuMP.@variable(model_1, z, Int)
    JuMP.@variable(model_1, u[1:3], Bin)
    JuMP.@variable(model_1, fi == 9)
    JuMP.@objective(model_1, Max, a - b + 2a1 - 10x)
    JuMP.@constraint(model_1, a + b - 10c - 2x + c1 <= 1)
    JuMP.@constraint(model_1, a * b <= 2)
    JuMP.@constraint(model_1, [1 - a; u] in JuMP.SecondOrderCone())

    model_2 = MyModel()
    JuMP.@variable(model_2, x, Bin)
    JuMP.@variable(model_2, y, Int)
    JuMP.@constraint(model_2, x * y <= 1)

    model_3 = MyModel()
    JuMP.@variable(model_3, x)
    JuMP.@constraint(model_3, x <= 3)

    Test.@test JuMP.model_string(MIME("text/plain"), model_1) == """
Max a - b + 2 a1 - 10 x
Subject to
 a + b - 10 c - 2 x + c1 $(repl(:leq)) 1
 a*b $(repl(:leq)) 2
 [-a + 1, u[1], u[2], u[3]] $(repl(:in)) MathOptInterface.SecondOrderCone(4)
"""

    Test.@test JuMP.model_string(MIME("text/latex"), model_1) == """
\$\$ \\begin{aligned}
\\max\\quad & a - b + 2 a1 - 10 x\\\\
\\text{Subject to} \\quad & a + b - 10 c - 2 x + c1 \\leq 1\\\\
 & a\\times b \\leq 2\\\\
 & [-a + 1, u_{1}, u_{2}, u_{3}] \\in \\text{MathOptInterface.SecondOrderCone(4)}\\\\
\\end{aligned} \$\$"""

    Test.@test sprint(show, model_1) == """
An Abstract JuMP Model
Maximization problem with:
Variables: 13
Objective function type: $(JuMP.GenericAffExpr{Float64,MyVariableRef})
Constraints: 3
Names registered in the model: a, a1, b, b1, c, c1, fi, u, x, y, z"""

    Test.@test sprint(show, model_2) == """
An Abstract JuMP Model
Feasibility problem with:
Variables: 2
Constraint: 1
Names registered in the model: x, y"""

    Test.@test sprint(show, model_3) == """
An Abstract JuMP Model
Feasibility problem with:
Variable: 1
Constraint: 1
Names registered in the model: x"""
    return
end

end
