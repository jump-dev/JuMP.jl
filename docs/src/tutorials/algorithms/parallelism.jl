# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Parallelism

# The purpose of this tutorial is to give a brief overview of parallelism in
# Julia as it pertains to JuMP, and to explain some of the things to be aware of
# when writing parallel algorithms involving JuMP models.

# ## Multi-threading and Distributed computing

# There are two main types of parallelism in Julia:
#
#  1. Multi-threading
#  2. Distributed computing

# In multi-threading, multiple tasks are run in a single Julia process and share
# the same memory. In distributed computing, tasks are run in multiple Julia
# processes with independent memory spaces. This can include processes across
# multiple physical machines, such as in a high-performance computing cluster.

# Choosing and understanding the type of parallelism you are using is important
# because the code you write for each type is different, and there are different
# limitations and benefits to each approach.

# ## Multi-threading

# To use multi-threading with Julia, you must either start Julia with the
# command line flag `--threads=N`, or you must set the `JULIA_NUM_THREADS`
# environment variable before launching Julia. For this documentation, we set
# the environment variable to:

ENV["JULIA_NUM_THREADS"]

# You can check how many threads are available using:

Threads.nthreads()

# The easiest way to use multi-threading in Julia is by placing the
# `Threads.@threads` macro in front of a `for`-loop:

ids = Int[]
@time begin
    Threads.@threads for i in 1:Threads.nthreads()
        push!(ids, Threads.threadid())
        sleep(1.0)
    end
end

# This for-loop sleeps for `1` second on each iteration. Thus, if it had
# executed sequentially, it should have taken the same number of seconds as
# there are threads availabe. Instead, it took only 1 second, showing that the
# iterations were executed simultaneously. We can verify this by checking the
# `Threads.threadid()` of the thread that executed each iteration:

ids

# !!! tip
#     For more information, read the Julia documentation
#     [Multi-Threading](https://docs.julialang.org/en/v1/manual/multi-threading/#man-multithreading).

# ## Distributed computing

# To use distributed computing with Julia, use the `Distributed` package:

import Distributed

# Like multi-threading, we need to tell Julia how many processes to add. We can
# do this either by starting Jlia with the `-p N` command line argument, or by
# using `Distributed.addprocs`:

import Pkg
project_path = Pkg.project().path
Distributed.addprocs(4; exeflags = "--project=$project_path")

# !!! warning
#     Not loading the environment with `--project` is a common mistake.

# The added processes are "worker" processes that we can use to do computation
# with. They are orchestrated by the process with the id `1`. You can check
# what process the code is currently running on using `Distributed.myid()`

Distributed.myid()

# As a general rule, to get maximum performance you should add as many processes
# as you have logical cores avaiable.

# Unlike the `for`-loop approach of multi-threading, distributed computing
# extends the Julia `map` function to a "parallel-map" function
# `Distributed.pmap`. For each element in the list of arguments to map over,
# Julia will copy the element to an idle worker process and evaluate the
# function passing the element as an input argument.

function hard_work(i::Int)
    sleep(1.0)
    return Distributed.myid()
end

try                         #hide
    @time ids = Distributed.pmap(hard_work, 1:4)
catch err                   #hide
    showerror(stderr, err)  #hide
end                         #hide

# Unforunately, if you try this code directly, you will get an error message
# that says `On worker 2: UndefVarError: hard_work not defined`. The error is
# thrown because, although process `1` knows what the `hard_work` function is,
# the worker processes do not know what it is.

# To fix the error, we need to use `Distributed.@everywhere`, which evaluates
# the code on every process:

Distributed.@everywhere begin
    function hard_work(i::Int)
        sleep(1.0)
        return Distributed.myid()
    end
end

# Now if we run `pmap`, we see that it took only 1 second instead of 4, and that
# it executed on each of the worker processes:

@time ids = Distributed.pmap(hard_work, 1:4)

# !!! tip
#     For more information, read the Julia documentation
#     [Distributed Computing](https://docs.julialang.org/en/v1/manual/distributed-computing/).

# ## Using parallelism the wrong way

# ### When building a JuMP model

# For large problems, building the model in JuMP can be a bottleneck, and you
# may consider trying to write code that builds the model in parallel, for
# example, by wrapping a `for`-loop that adds constraints with
# `Threads.@threads`.

# Unfortunately, you cannot build a JuMP model in parallel, and attempting to do
# so may error or produce incorrect results.

# In most cases, we find that the reason for the bottleneck is not JuMP, but in
# how you are constructing the problem data, and that with changes, it is
# possible to build a model in a way that is not the bottleneck in the solution
# process.

# !!! tip
#     Looking for help to make your code run faster? Ask for help on the
#     [community forum](https://discourse.julialang.org/c/domain/opt/13). Make
#     sure to include a reproducible example of your code.

# ### With a single JuMP model

# A common approach people try is to use parallelism with a single JuMP model.
# For example, to optimize a model over multiple right-hand side vectors, you
# may try:

# ```julia
# using JuMP
# import HiGHS
# model = Model(HiGHS.Optimizer)
# set_silent(model)
# @variable(model, x)
# @objective(model, Min, x)
# solutions = Float64[]
# Threads.@threads for i in 1:10
#     set_lower_bound(x, i)
#     optimize!(model)
#     push!(solutions, objective_value(model))
# end

# This will not work, and attempting to do so may error, crash Julia or produce
# incorrect results. The reason is because most solvers are written in C,

# ## Using parallelism the right way

# To use parallelism with JuMP, the simplest rule to remember is that each
# worker, both threads and distributed processes, must all have their own
# instance of a JuMP model.

# ### With multi-threading

# With multi-threading, create a new instance of `model` in each iteration of
# the for-loop:

using JuMP
import HiGHS
solutions = Float64[]
Threads.@threads for i in 1:10
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x)
    @objective(model, Min, x)
    set_lower_bound(x, i)
    optimize!(model)
    push!(solutions, objective_value(model))
end
solutions

# ### With distributed computing

# With distributed computing, remember to evaluate all of the code on all of the
# processes using `Distributed.@everywhere`, and then write a function which
# creates a new instance of the model on every evaluation:

Distributed.@everywhere begin
    using JuMP
    import HiGHS
end

Distributed.@everywhere begin
    function solve_model_with_right_hand_side(i)
        model = Model(HiGHS.Optimizer)
        set_silent(model)
        @variable(model, x)
        @objective(model, Min, x)
        set_lower_bound(x, i)
        optimize!(model)
        return objective_value(model)
    end
end

solutions = Distributed.pmap(solve_model_with_right_hand_side, 1:10)

