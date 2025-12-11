# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Multi-objective project planning

# The purpose of this tutorial is to demonstrate how to create and solve a
# multi-objective linear program.

# **This tutorial was originally contributed by Xavier Gandibleux. It was
# adapted from an example created by Begoña Vitoriano from the Universidad
# Complutense de Madrid and originally presented at the EURO PhD School on
# multicriteria decision making with mathematical programming, held February
# 17-28, 2014 in Madrid.**

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

# This tutorial investigates the problem of project planning. To complete the
# project, a number of tasks need to be completed in some precedence order.
# Each task requires a different amount of resources, for example time, money,
# and labor. The optimal schedule depends on the objective function. Do we want
# to minimize the makespan (the time to complete the project)? Or do we want to
# complete the project while minimizing the cost to complete it? Instead of
# choosing a single objective function, this tutorial solves a project planning
# problem as a multi-objective optimization problem in which we minimize the
# three objectives of time, cost, and labor units. Then, the decision maker can
# select a schedule from the set of possible solutions that best trades off the
# three objectives according to their personal preferences.

# ## Data

# We study how to optimize the project schedule for building a house. There are
# many tasks that need to be completed to build the house. Each task has a given
# duration and a number staff required to perform it. Some tasks can be
# accelerated (for a fixed cost per day) by hiring 20% more people (rounding
# up). Tasks have a precedence order. Here's the data we're going to use:

filename = joinpath(@__DIR__, "multi_objective_project_planning.csv")
tasks = CSV.read(filename, DataFrames.DataFrame)

# ## JuMP formulation

# We can build a JuMP model to optimize the project plan.

model = Model();

# The first decision variable is when tasks start. We can compute an upper bound
# for starting any task by summing the total duration of all tasks.

start_ub = sum(tasks.duration)
@variable(model, 0 <= x_start_time[tasks.task] <= start_ub, Int);

# A second decision is to decide how many days to accelerate each task:

@variable(model, x_accelerate[tasks.task] >= 0, Int);

# We need a corresponding binary variable to track if a task is accelerated:

@variable(model, z_accelerate[tasks.task], Bin);

# and a constraint to link `x_accelerate` to `z_accelerate`:

@constraint(
    model,
    [r in eachrow(tasks)],
    x_accelerate[r.task] <= r.max_accelerate * z_accelerate[r.task],
);

# There are three objectives: makespan, acceleration cost, and the number of
# extra staff required.

# The makespan is the start time of the "Project complete" task:

@expression(model, obj_time, x_start_time["S"]);

# The cost is the cost per day of accelerating a task:

@expression(
    model,
    obj_cost,
    sum(r.cost_accelerate * x_accelerate[r.task] for r in eachrow(tasks)),
);

# The number of extra staff is the 20% increase in staff if a task is
# accelerated:

@expression(
    model,
    obj_people,
    sum(
        ceil(Int, 0.2 * r.staff) * z_accelerate[r.task] for r in eachrow(tasks)
    )
);

# Our full objective is therefore:

@objective(model, Min, [obj_time, obj_cost, obj_people])

# For the constraints, the start time of a task must be after any required
# preceding tasks have finished:

for r in eachrow(tasks)
    for dependent_task in split(r.precedes, ","; keepempty = false)
        @constraint(
            model,
            x_start_time[dependent_task] >=
            x_start_time[r.task] + r.duration - x_accelerate[r.task]
        )
    end
end

# ## Solution

# Now we have built our model, we can solve it using MOA:

set_optimizer(model, () -> MOA.Optimizer(HiGHS.Optimizer))
set_attribute(model, MOA.Algorithm(), MOA.KirlikSayin())
set_silent(model)
optimize!(model)

# The algorithm found 21 non-dominated solutions:

Test.@test result_count(model) == 21  #src
solution_summary(model)

# The [`objective_bound`](@ref) is the ideal point, that is, a lower bound on
# each objective when solved independently.

objective_bound(model)

# In the ideal case, we cannot do better than a makespan of 77, paying \$0 in
# extra cost, and employing zero extra staff. In reality, we cannot achieve this
# ideal outcome, and we need to trade-off the three objectives against each
# other.

# Here's a way to plot their relative trade-offs in objective space. Each point
# is a separate solution, numbered by their `result` index.

plt = Plots.scatter(
    [value(obj_time; result) for result in 1:result_count(model)],
    [value(obj_people; result) for result in 1:result_count(model)];
    marker_z = [
        value(obj_cost; result) / 1_000 for result in 1:result_count(model)
    ],
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

# Here's a function to visualize the schedule of a given solution. It's not
# important to understand the details of the plotting function.

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

# The darkest bars are tasks that are accelerated by two days. The lightest bars
# are not accelerated.

# The solution with the longest makespan is:

plot_solution(; result = 21)

# No tasks are accelerated.

# There are a variety of solutions in between these extremes. One example is:

plot_solution(; result = 8)
