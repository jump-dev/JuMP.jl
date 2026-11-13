# Copyright (c) 2019 Arpit Bhatia and contributors                               #src
#                                                                                #src
# Permission is hereby granted, free of charge, to any person obtaining a copy   #src
# of this software and associated documentation files (the "Software"), to deal  #src
# in the Software without restriction, including without limitation the rights   #src
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      #src
# copies of the Software, and to permit persons to whom the Software is          #src
# furnished to do so, subject to the following conditions:                       #src
#                                                                                #src
# The above copyright notice and this permission notice shall be included in all #src
# copies or substantial portions of the Software.                                #src
#                                                                                #src
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #src
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       #src
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    #src
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         #src
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  #src
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  #src
# SOFTWARE.                                                                      #src

# # Solvers and Solutions

# **Originally Contributed by**: Arpit Bhatia

# The purpose of this part of the tutorial is to introduce you to solvers and
# how to use them with JuMP. We'll also learn what to do with a problem after
# the solver has finished optimizing it.

# ## What is a Solver?

# A solver is a software package that incorporates algorithms for finding
# solutions to one or more classes of problem. For example, GLPK, which we used
# in the previous tutorials is a solver for linear programming (LP) and mixed
# integer programming (MIP) problems. It incorporates algorithms such as the
# simplex method, interior-point method etc. JuMP currently supports a number of
# open-source and commercial solvers which can be viewed in the [Supported-solvers](@ref)
# table.

# ## What is MathOptInterface?

# Each mathematical optimization solver API has its own concepts and data
# structures for representing optimization models and obtaining results.
# However, it is often desirable to represent an instance of an optimization
# problem at a higher level so that it is easy to try using different solvers.

# MathOptInterface (MOI) is an abstraction layer designed to provide an
# interface to mathematical optimization solvers so that users do not need to
# understand multiple solver-specific APIs. MOI can be used directly, or through
# a higher-level modeling interface like JuMP.

# Note that JuMP re-exports MathOptInterface and you can use the shortcut `MOI`
# to refer to MathOptInterface in your code.

# ## Constructing a model

# JuMP models can be created in three different modes: `AUTOMATIC`, `MANUAL` and
# `DIRECT`. We'll use the following LP to illustrate them.

# ```math
# \begin{aligned}
# & \max_{x,y} & x + 2y \\
# & \;\;\text{s.t.} & x + y &\leq 1 \\
# & & 0\leq x, y &\leq 1 \\
# \end{aligned}
# ```

using JuMP
using GLPK

# ### `AUTOMATIC` Mode

# #### With Optimizer

# This is the easiest method to use a solver in JuMP. In order to do so, we
# simply set the solver inside the Model constructor.

model_auto = Model(GLPK.Optimizer)
@variable(model_auto, 0 <= x <= 1)
@variable(model_auto, 0 <= y <= 1)
@constraint(model_auto, x + y <= 1)
@objective(model_auto, Max, x + 2y)
optimize!(model_auto)
objective_value(model_auto)

# #### No Optimizer (at first)

# It is also possible to create a JuMP model with no optimizer attached. After
# the model object is initialized empty and all its variables, constraints and
# objective are set, then we can attach the solver at `optimize!` time.

model_auto_no = Model()
@variable(model_auto_no, 0 <= x <= 1)
@variable(model_auto_no, 0 <= y <= 1)
@constraint(model_auto_no, x + y <= 1)
@objective(model_auto_no, Max, x + 2y)
set_optimizer(model_auto_no, GLPK.Optimizer)
optimize!(model_auto_no)
objective_value(model_auto_no)

# Note that we can also enforce the automatic mode by passing
# `caching_mode = MOIU.AUTOMATIC` in the Model function call.

# ### `MANUAL` Mode

# This mode is similar to the `AUTOMATIC` mode, but there are less protections
# from the user getting errors from the solver API. On the other side, nothing
# happens silently, which might give the user more control. It requires
# attaching the solver before the solve step using the `MOIU.attach_optimizer()`
# function.

model_manual = Model(GLPK.Optimizer, caching_mode = MOIU.MANUAL)
@variable(model_manual, 0 <= x <= 1)
@variable(model_manual, 0 <= y <= 1)
@constraint(model_manual, x + y <= 1)
@objective(model_manual, Max, x + 2y)
MOIU.attach_optimizer(model_manual)
optimize!(model_manual)
objective_value(model_manual)

# ### `DIRECT` Mode

# Some solvers are able to handle the problem data directly. This is common for
# LP/MIP solver but not very common for open-source conic solvers. In this case
# we do not set a optimizer, we set a backend which is more generic and is able
# to hold data and not only solving a model.

model_direct = direct_model(GLPK.Optimizer())
@variable(model_direct, 0 <= x <= 1)
@variable(model_direct, 0 <= y <= 1)
@constraint(model_direct, x + y <= 1)
@objective(model_direct, Max, x + 2y)
optimize!(model_direct)
objective_value(model_direct)

# ### Solver Options

# Many of the solvers also allow options to be passed in. However, these options
# are solver-specific. To find out the various options available, please check
# out the individual solver packages. Some examples for the GLPK solver are
# given below.

using GLPK

# To turn off printing (i.e. silence the solver),

model = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0));

# To increase the maximum number of simplex iterations:

model = Model(optimizer_with_attributes(GLPK.Optimizer, "it_lim" => 10_000));

# To set the solution timeout limit (in milliseconds):

model = Model(optimizer_with_attributes(GLPK.Optimizer, "tm_lim" => 5_000));

# ## How to querying the solution

# So far we have seen all the elements and constructs related to writing a JuMP
# optimization model. In this section we reach the point of what to do with a
# solved problem. JuMP follows closely the concepts defined in MathOptInterface
# to answer user questions about a finished call to `optimize!(model)`. The
# three main steps in querying a solution are given below. We'll use the model
# we created in `AUTOMATIC` mode with an optimizer attached in this section.

# ### The termination status

# Termination statuses are meant to explain the reason why the optimizer stopped
# executing in the most recent call to `optimize!`.

termination_status(model_auto)

# You can view the different termination status codes by referring to the docs
# or though checking the possible types using the below command.

display(typeof(MOI.OPTIMAL))

# ### The primal and dual status

# These statuses indicate what kind of result is available to be queried from
# the model. It's possible that no result is available to be queried. We shall
# discuss more on the dual status and solutions in the Duality tutorial.

primal_status(model_auto)

#-

dual_status(model_auto)

# As we saw before, the result (solution) status codes can be viewed directly
# from Julia.

display(typeof(MOI.FEASIBLE_POINT))

# ### Getting the primal solution

# Provided the primal status is not `MOI.NO_SOLUTION`, we can inspect the
# solution values and optimal cost.

value(x)

#-

value(y)

#-

objective_value(model_auto)

# Since it is possible that no solution is available to be queried from the
# model, calls to [`value`](@ref) may throw errors. Hence, it is recommended to
# check for the presence of solutions.

model_no_solution = Model(GLPK.Optimizer)
@variable(model_no_solution, 0 <= x <= 1)
@variable(model_no_solution, 0 <= y <= 1)
@constraint(model_no_solution, x + y >= 3)
@objective(model_no_solution, Max, x + 2y)

optimize!(model_no_solution)

try #hide
if termination_status(model_no_solution) == MOI.OPTIMAL
    optimal_solution = value(x)
    optimal_objective = objective_value(model_no_solution)
elseif termination_status(model_no_solution) == MOI.TIME_LIMIT && has_values(model_no_solution)
    suboptimal_solution = value(x)
    suboptimal_objective = objective_value(model_no_solution)
else
    error("The model was not solved correctly.")
end
catch err; showerror(stderr, err); end  #hide
