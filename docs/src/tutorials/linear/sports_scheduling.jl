# Copyright (c) 2026 Oscar Dowson and contributors                               #src
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

# # Sport scheduling

# **This tutorial was adapted from an [example written by Eli Towle](https://publicsectororcourse.wordpress.com/2016/05/09/optimization-in-julia/)
# in 2016 as part of an operations research course at the University of
# Wisconsin-Madison.**

# The purpose of this tutorial is to demonstrate a simple model for scheduling
# round-robin tournaments. As teams, it uses the [Big 10](https://en.wikipedia.org/wiki/Big_Ten_Conference).
# (You might notice that there are more than 10 teams. Our example was also
# written before the expansion of the Conference in 2024.)

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import HiGHS

# ## Basic model

# Here are the teams in our tournament:

M = [
    #!format:off
    "Ind", "UMD", "UMich", "MSU", "OSU", "Penn", "Rtgrs", "Ill", "Iowa", "UMN",
    "UNL", "NU", "Purd", "UW",
    #!format: on
];

# For each team to play each other exactly once, we need the number of teams - 1
# weeks:

T = length(M) - 1;

# Now we create a JuMP model to build our optimzation problem:

model = Model(HiGHS.Optimizer);

# Variable: a binary which is true for `x[m,n,t]` if team `m` is playing `n` at
# home in week `t`.

@variable(model, x[M, M, 1:T], Bin);

# Constraint: each team `m` can never play themselves.

@constraint(model, [m in M, t in 1:T], x[m, m, t] == 0);

# Constraint: each team `m` can play at most once per day

@constraint(model, [m in M, t in 1:T], sum(x[m, :, t]) + sum(x[:, m, t]) <= 1);

# Constraint: every team `m` plays at least half home games

@constraint(model, [m in M], div(T, 2) <= sum(x[m, :, :]) <= div(T, 2) + 1);

# Constraint: no more than two away games in any three-game window

@constraint(model, [m in M, t in 1:(T-2)], sum(x[:, m, t:(t+2)]) <= 2);

# Constraint: every team must play every other team exactly once

@constraint(
    model,
    [m in M, n in M; m != n],
    sum(x[m, n, :]) + sum(x[n, m, :]) == 1,
);

# Now we can solve our model:

set_silent(model)
optimize!(model)
assert_is_solved_and_feasible(model)

# and print the schedule:

function print_schedule(M::Vector{String}, T::Int, Y::AbstractArray{Bool})
    println("Week ", join(rpad.(M, 6), ' '))
    for t in 1:T
        print(rpad(t, 5))
        for m in M, n in M
            if Y[m, n, t]
                print(rpad(n, 7))
            elseif Y[n, m, t]
                print("@", rpad(n, 6))
            end
        end
        println()
    end
    return
end

Y = round.(Bool, value.(x))
print_schedule(M, T, Y)

# This schedule is okay, but it features a large number of back-to-back away
# games. Let's count them:

number_of_back_to_back_away_games =
    sum(round(Int, value(sum(x[:, m, (t-1):t]))) == 2 for m in M, t in 2:T)

# A better schedule would minimize this quantity.

# ## Minimizing the number of back-to-back away games

# To minimize the number of back-to-back away games, we modify our model.

# Variable: a binary which is true for `y[m,t]` if team `m` is playing
# back-to-back away games in week `t`.

@variable(model, y[M, 2:T], Bin);

# Objective: minimize the number of back-to-back away games

@objective(model, Min, sum(y));

# Constraint: count back-to-back away games

@constraint(model, [m in M, t in 2:T], y[m, t] >= sum(x[:, m, (t-1):t]) - 1);

# Now we can solve our model. However, this problem is actually very difficult
# to solve to optimality. Rather than wait a very long time, we set a time limit
# so that this documentation doesn't take too long to build:

set_time_limit_sec(model, 30.0)

# We're also going to set a start value based on the previous solution:

set_start_value.(x, Y)
for m in M, t in 2:T
    set_start_value(y[m, t], sum(Y[:, m, (t-1):t]) - 1)
end

# Now we can optimize:

optimize!(model)

# Because we hit a time limit, we can't use [`assert_is_solved_and_feasible`](@ref),
# but we still check that we found a feasible primal solution:

@assert termination_status(model) == TIME_LIMIT
@assert primal_status(model) == FEASIBLE_POINT

# This solution has fewer back-to-back away games:

number_of_back_to_back_away_games =
    sum(round(Int, value(sum(x[:, m, (t-1):t]))) == 2 for m in M, t in 2:T)

# And the final schedule is:

Y = round.(Bool, value.(x))
print_schedule(M, T, Y)

# Try running the model for longer. What is the smallest number of back-to-back
# away games you can find?
