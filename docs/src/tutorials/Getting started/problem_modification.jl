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

# # Problem Modification

# **Originally Contributed by**: Arpit Bhatia

# This tutorial deals with how to modify models after they have been created and
# solved. This functionality can be useful, for example, when we are solving
# many similar models in succession or generating the model dynamically.
# Additionally it is sometimes desirable for the solver to re-start from the
# last solution to reduce running times for successive solves. This is called
# warm-starting.

using JuMP

# ## Modifying variables

model = Model()
@variable(model, x);

# ### Variable bounds

# The [`set_lower_bound`](@ref) and [`set_upper_bound`](@ref) functions can be
# used to create as well as modify an existing variable bound.

set_lower_bound(x, 3)
lower_bound(x)

#-

set_lower_bound(x, 2)
lower_bound(x)

# We can delete variable bounds using the [`delete_lower_bound`](@ref) and
# [`delete_upper_bound`](@ref) functions.

delete_lower_bound(x)
has_lower_bound(x)

# We can assign a fixed value to a variable using [`fix`](@ref).

fix(x, 5)
fix_value(x)

# However, fixing a variable with existing bounds will throw an error.

@variable(model, y >= 0);

try #hide
fix(y, 2)
catch err; showerror(stderr, err); end  #hide

# As we can see in the error message above, we have to specify to JuMP that we
# wish to forcefuly remove the bound.

fix(y, 2; force = true)
fix_value(y)

# We can also call the [`unfix`](@ref) function to remove the fixed value.

unfix(x)
is_fixed(x)

# ### Deleting variables

# The [`all_variables`](@ref) function returns a list of all variables present
# in the model.

all_variables(model)

# To delete variables from a model, we can use the [`delete`](@ref) function.

delete(model, x)
all_variables(model)

# We can also check whether a variable is a valid JuMP variable in a model using
# the [`is_valid`](@ref) function.

is_valid(model, x)

# ## Modifying constraints

model = Model()
@variable(model, x);

# ### Modifying a variable xoefficient

# It is also possible to modify the scalar coefficients (but notably not the
# quadratic coefficients) using the [`set_normalized_coefficient`](@ref)
# function.

@constraint(model, con, 2x <= 1);

#-

set_normalized_coefficient(con, x, 3)

#-

con

# ### Deleting a constraint

# Just like for deleting variables, we can use the [`delete`](@ref) function for
# constraints as well.

delete(model, con)
is_valid(model, con)

# ## Modifying the objective

model = Model()
@variable(model, x)
@objective(model, Min, 7x + 4);

# The function [`objective_function`](@ref) is used to query the objective of a
# model.

objective_function(model)

# [`objective_sense`](@ref) is similarily used to query the objective sense of a
# model.

objective_sense(model)

# To easiest way to change the objective is to simply call [`@objective`](@ref)
# again, which will replace the objective function and objective sense.

@objective(model, Max, 8x + 3)

#-

objective_function(model)

#-

objective_sense(model)

# Another way is to change the objective is to use the low-level functions
# [`set_objective_function`](@ref) and [`set_objective_sense`](@ref).

set_objective_function(model, 5x + 11)
objective_function(model)

#-

set_objective_sense(model, MOI.MIN_SENSE)
objective_sense(model)

# Note that we can't use the `Min` and `Max` shortcuts here as
# [`set_objective_sense`](@ref) is a low level function.
