# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Multi-objective project planning

# The purpose of this tutorial is to demonstrate how to create and solve a
# multi-objective linear program.

# **This tutorial was originally contributed by Xavier Gandibleux. It was
# adapted from the EURO PhD School on MCDM: Case studies on MCDM with
# Mathematical Programming (Begoña Vitoriano). Complutense University, Madrid,
# Feb 17-28 2014.**

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import CSV
import DataFrames
import HiGHS
import MultiObjectiveAlgorithms as MOA
import Plots
import Test  #src

# ## Background

# In a business or industry is common that for tasks to be developed a certain
# total or partial order to carry them out is established. These tasks can
# belong to a project and what is sought is to determine the shortest possible
# time to finish the project. Sometimes, the production processes require the
# use of certain resources and it must be decided how to sequence them with
# different objectives or goals to achieve. What is common to both approaches is
# the existence of tasks and precedence relations between them, with or without
# resources. Project planning has its own entity and widespread use, so it has
# therefore been subject to detailed analysis and has led to specific
# techniques, although limited.

# The problem stated to illustrate this type of problems is building a house.
# For planning purposes the following activities are considered with their
# durations and people needed to perform it. Some tasks can be accelerated by a
# maximum of 2 days, except for M and R which only support 1 day of
# acceleration. Acceleration can be achieved by including 20% more people
# (rounding up). Tasks have a precedence order. Here's the data:

filename = joinpath(@__DIR__, "multi_objective_project_planning.csv")
tasks = CSV.read(filename, DataFrames.DataFrame)

# We can build a JuMP model to optimize the project plan.

model = Model();

# The first decision variable is when tasks start. We can compute an upper bound
# for starting any task by summing the total duration of all tasks.
start_ub = sum(tasks.duration)
@variable(model, 0 <= x_start_time[tasks.task] <= start_ub, Int)

# A second decision is to decide how many days to accelerate each task:
@variable(model, x_accelerate[tasks.task] >= 0, Int)

# We need a corresponding binary variable to track if a task is accelerated:
@variable(model, z_accelerate[tasks.task], Bin)

# and a constraint to link `x_accelerate` to `z_accelerate`:
@constraint(
    model,
    [r in eachrow(tasks)],
    x_accelerate[r.task] <= r.max_accelerate * z_accelerate[r.task],
)

# There are three objectives: makespan, cost, and staffing level.

# The makespan is the start time of the "Project complete" task:
@expression(model, obj_time, x_start_time["S"])

# The cost is the cost per day of accelerating a task:
@expression(
    model,
    obj_cost,
    sum(r.cost_accelerate * x_accelerate[r.task] for r in eachrow(tasks)),
)

# The number of people is the 20% increase in staff if a task is accelerated:
@expression(
    model,
    obj_people,
    sum(
        ceil(Int, 1.2 * r.staff) * z_accelerate[r.task] for r in eachrow(tasks)
    )
)

@objective(model, Min, [obj_time, obj_cost, obj_people])

# For the constraints, the start time of a task must be after any dependent
# tasks have finished.
for r in eachrow(tasks)
    for dependent_task in split(r.predecence, ","; keepempty = false)
        @constraint(
            model,
            x_start_time[dependent_task] >=
            x_start_time[r.task] + r.duration - x_accelerate[r.task]
        )
    end
end

set_optimizer(model, () -> MOA.Optimizer(HiGHS.Optimizer))
set_attribute(model, MOA.Algorithm(), MOA.KirlikSayin())
set_silent(model)
optimize!(model)
Test.@test result_count(model) == 25  #src
solution_summary(model)

# The algorithm found many solutions. Here's a way to plot their relative
# trade-offs in objective space:

plt = Plots.scatter(
    [value(obj_time; result) for result in 1:result_count(model)],
    [value(obj_people; result) for result in 1:result_count(model)];
    marker_z = [value(obj_cost; result) / 1_000 for result in 1:result_count(model)],
    xlabel = "Makespan",
    ylabel = "Extra staff required",
    colorbar_title = "Cost ('000 €)",
    label = nothing,
    markersize = 8,
    markerstrokewidth = 0,
)
for i in 1:result_count(model)
    obj = objective_value(model; result = i)
    Plots.annotate!(plt, obj[1]+0.5, obj[3]+0.5, (i, 6))
end
plt
# If we decrease the extra staff required, the makespan increases. Additionally,
# there are some solutions where we can keep the same makespan, but increase the
# number of staff for a reduction in cost. These solutions arise because we are
# accelerating cheaper tasks instead of more expensive ones.

# Here's a function to visualize the solution:

function plot_solution(; result::Int)
    shape(; x, y, w) = Plots.Shape(x .+ [0, w, w, 0], y .+ [0.1, 0.1, 0.9, 0.9])
    shapes = [
        shape(;
            x = value(x_start_time[r.task]; result),
            y = i,
            w = r.duration - value(x_accelerate[r.task]; result),
        ) for (i, r) in enumerate(eachrow(tasks))
    ]
    accelerate_color = Dict(0 => "#a0b1ec", 1 => "#4063d8", 2 => "#20326c")
    colors = map(eachrow(tasks)) do r
        return accelerate_color[round(Int, value(x_accelerate[r.task]; result))]
    end
    obj = round.(Int, objective_value(model; result))
    obj[2] /= 1_000
    return Plots.plot(
        shapes;
        color = permutedims(colors),
        linewidth = 0,
        legend = false,
        title = "makespan = $(obj[1]), cost = €$(obj[2])k, people =$(obj[3])",
        xlabel = "Time",
        ylabel = "Task",
        yticks = (1.5:19.5, 'A':'S'),
    )
end

# The solution with the shortest makespan is:

plot_solution(; result = 1)

# The solution with the longest makespan is:

plot_solution(; result = 25)

# There are a variety of solutions in between these extremes. One example is:

plot_solution(; result = 8)
