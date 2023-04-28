# Parallelism

The purpose of this tutorial is to give a brief overview of parallelism in
Julia as it pertains to JuMP, and to explain some of the things to be aware of
when writing parallel algorithms involving JuMP models.

## Multi-threading and distributed computing

There are two main types of parallelism in Julia:

 1. Multi-threading
 2. Distributed computing

In multi-threading, multiple tasks are run in a single Julia process and share
the same memory. In distributed computing, tasks are run in multiple Julia
processes with independent memory spaces. This can include processes across
multiple physical machines, such as in a high-performance computing cluster.

Choosing and understanding the type of parallelism you are using is important
because the code you write for each type is different, and there are different
limitations and benefits to each approach. However, the best choice is highly
problem dependent, so you may want to experiment with both approaches to
determine what works for your situation.

### Multi-threading

To use multi-threading with Julia, you must either start Julia with the
command line flag `--threads=N`, or you must set the `JULIA_NUM_THREADS`
environment variable before launching Julia. For this documentation, we set
the environment variable to:

````julia
julia> ENV["JULIA_NUM_THREADS"]
"4"
````

You can check how many threads are available using:

````julia
julia> Threads.nthreads()
4
````

The easiest way to use multi-threading in Julia is by placing the
`Threads.@threads` macro in front of a `for`-loop:

````julia
julia> @time begin
           ids = Int[]
           my_lock = Threads.ReentrantLock()
           Threads.@threads for i in 1:Threads.nthreads()
               global ids, my_lock
               Threads.lock(my_lock) do
                   push!(ids, Threads.threadid())
               end
               sleep(1.0)
           end
       end
  1.037087 seconds (31.32 k allocations: 1.836 MiB, 2.02% compilation time)
````

This for-loop sleeps for `1` second on each iteration. Thus, if it had
executed sequentially, it should have taken the same number of seconds as
there are threads available. Instead, it took only 1 second, showing that the
iterations were executed simultaneously. We can verify this by checking the
`Threads.threadid()` of the thread that executed each iteration:

````julia
julia> ids
4-element Vector{Int64}:
 2
 4
 1
 3
````

!!! warn
    When working with threads, you need to avoid race conditions, in which two
    threads attempt to write to the same variable at the same time. In the above
    example we avoided a race condition by using `ReentrantLock`. See the
    [Multi-threading](https://docs.julialang.org/en/v1/manual/multi-threading/)
    section of the Julia documentation for more details.

### Distributed computing

To use distributed computing with Julia, use the `Distributed` package:

````julia
julia> import Distributed
````

Like multi-threading, we need to tell Julia how many processes to add. We can
do this either by starting Julia with the `-p N` command line argument, or by
using `Distributed.addprocs`:

````julia
julia> import Pkg

julia> project = Pkg.project();

julia> workers = Distributed.addprocs(4; exeflags = "--project=$(project.path)")
4-element Vector{Int64}:
 2
 3
 4
 5
````

!!! warning
    Not loading the parent environment with `--project` is a common mistake.

The added processes are "worker" processes that we can use to do computation
with. They are orchestrated by the process with the id `1`. You can check
what process the code is currently running on using `Distributed.myid()`

````julia
julia> Distributed.myid()
1
````

As a general rule, to get maximum performance you should add as many processes
as you have logical cores available.

Unlike the `for`-loop approach of multi-threading, distributed computing
extends the Julia `map` function to a "parallel-map" function
`Distributed.pmap`. For each element in the list of arguments to map over,
Julia will copy the element to an idle worker process and evaluate the
function, passing the element as an input argument.

````julia
julia> function hard_work(i::Int)
           sleep(1.0)
           return Distributed.myid()
       end
hard_work (generic function with 1 method)

julia> Distributed.pmap(hard_work, 1:4)
ERROR: On worker 2:
UndefVarError: #hard_work not defined
Stacktrace:
[...]
````

Unfortunately, if you try this code directly, you will get an error message
that says `On worker 2: UndefVarError: hard_work not defined`. The error is
thrown because, although process `1` knows what the `hard_work` function is,
the worker processes do not.

To fix the error, we need to use `Distributed.@everywhere`, which evaluates
the code on every process:

````julia
julia> Distributed.@everywhere begin
           function hard_work(i::Int)
               sleep(1.0)
               return Distributed.myid()
           end
       end
````

Now if we run `pmap`, we see that it took only 1 second instead of 4, and that
it executed on each of the worker processes:

````julia
julia> @time ids = Distributed.pmap(hard_work, 1:4)
  1.202006 seconds (216.39 k allocations: 13.301 MiB, 4.07% compilation time)
4-element Vector{Int64}:
 2
 3
 5
 4
````

!!! tip
    For more information, read the Julia documentation
    [Distributed Computing](https://docs.julialang.org/en/v1/manual/distributed-computing/).

## Using parallelism the wrong way

### When building a JuMP model

For large problems, building the model in JuMP can be a bottleneck, and you
may consider trying to write code that builds the model in parallel, for
example, by wrapping a `for`-loop that adds constraints with
`Threads.@threads`.

Unfortunately, you cannot build a JuMP model in parallel, and attempting to do
so may error or produce incorrect results.

In most cases, we find that the reason for the bottleneck is not JuMP, but in
how you are constructing the problem data, and that with changes, it is
possible to build a model in a way that is not the bottleneck in the solution
process.

!!! tip
    Looking for help to make your code run faster? Ask for help on the
    [community forum](https://jump.dev/forum). Make sure to include a
    reproducible example of your code.

### With a single JuMP model

A common approach people try is to use parallelism with a single JuMP model.
For example, to optimize a model over multiple right-hand side vectors, you
may try:

```julia
using JuMP
import HiGHS
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
@objective(model, Min, x)
solutions = Pair{Int,Float64}[]
my_lock = Threads.ReentrantLock()
Threads.@threads for i in 1:10
    set_lower_bound(x, i)
    optimize!(model)
    Threads.lock(my_lock) do
        push!(solutions, i => objective_value(model))
    end
end
```

This will not work, and attempting to do so may error, crash Julia or produce
incorrect results.

## Using parallelism the right way

To use parallelism with JuMP, the simplest rule to remember is that each
worker must have its own instance of a JuMP model.

### With multi-threading

With multi-threading, create a new instance of `model` in each iteration of
the for-loop:

````julia
julia> using JuMP

julia> import HiGHS

julia> solutions = Pair{Int,Float64}[]

julia> my_lock = Threads.ReentrantLock();

julia> Threads.@threads for i in 1:10
           model = Model(HiGHS.Optimizer)
           set_silent(model)
           set_attribute(model, MOI.NumberOfThreads(), 1)
           @variable(model, x)
           @objective(model, Min, x)
           set_lower_bound(x, i)
           optimize!(model)
           Threads.lock(my_lock) do
               push!(solutions, i => objective_value(model))
           end
       end

julia> solutions
10-element Vector{Pair{Int64, Float64}}:
  7 => 7.0
  4 => 4.0
  1 => 1.0
  9 => 9.0
  5 => 5.0
  8 => 8.0
 10 => 10.0
  2 => 2.0
  6 => 6.0
  3 => 3.0
````

!!! warning
    For some solvers, it may be necessary to limit the number of threads used
    internally by the solver to 1 by setting the [`MOI.NumberOfThreads`](@ref)
    attribute.

### With distributed computing

With distributed computing, remember to evaluate all of the code on all of the
processes using `Distributed.@everywhere`, and then write a function which
creates a new instance of the model on every evaluation:

````julia
julia> Distributed.@everywhere begin
           using JuMP
           import HiGHS
       end

julia> Distributed.@everywhere begin
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

julia> solutions = Distributed.pmap(solve_model_with_right_hand_side, 1:10)
10-element Vector{Float64}:
  1.0
  2.0
  3.0
  4.0
  5.0
  6.0
  7.0
  8.0
  9.0
 10.0
````

## Other types of parallelism

### GPU

JuMP does not support GPU programming, and few solvers support execution on a
GPU.

### Parallelism within the solver

Many solvers use parallelism internally. For example, commercial solvers like
[Gurobi](https://github.com/jump-dev/Gurobi.jl) and [CPLEX](https://github.com/jump-dev/CPLEX.jl)
both parallelize the search in branch-and-bound. Solvers supporting internal
parallelism will typically support the [`MOI.NumberOfThreads`](@ref) attribute,
which you can set using [`set_attribute`](@ref).
