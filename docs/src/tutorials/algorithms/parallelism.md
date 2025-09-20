# Parallelism

The purpose of this tutorial is to give a brief overview of parallelism in
Julia as it pertains to JuMP, and to explain some of the things to be aware of
when writing parallel algorithms involving JuMP models.

## Overview

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

## Multi-threading

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

!!! danger
    The `Threads.threadid()` that a task runs on may change during execution.
    Therefore, it is not safe to use `Threads.threadid()` to index into, say, a
    vector of buffer or stateful objects. As an example, do not do:
    ```julia
    x = rand(Threads.nthreads())
    Threads.@threads for i in 1:Threads.nthreads()
        x[Threads.threadid()] *= 2  # Danger! This use of threadid is not safe
    end
    ```
    For more information, read
    [PSA: Thread-local state is no longer recommended](https://julialang.org/blog/2023/07/PSA-dont-use-threadid/).

### Data races

When working with threads, you must avoid data races. A data race occurs when
multiple threads access the same variable at the same time, at least one thread
modifies the variable, and the order of the reads and writes are not properly
coordinated.

Here's an example of a data race:
````julia
julia> begin
           x = Ref(0)
           Threads.@threads for i in 1:Threads.nthreads()
                for i in 1:1_000
                    x[] += 1
                end
           end
           x[]
       end
1106
````
The expected answer is `4_000` (because there are four threads each incrementing
1,000 times), but the actual result is much smaller. Moreover, the result is
non-deterministic; if we re-ran this code, we would get a different value each
time.

We got the wrong answer because multiple threads are reading and writing `x`
at the same time without coordination. For example, the following sequence
could occur:

 * Assume `x[]` currently has a value of `3`
 * Thread A reads `x[]` to get `3`
 * Thread B reads `x[]` to get `3`
 * Thread A writes `x[] = 3 + 1 = 4`
 * Thread B writes `x[] = 3 + 1 = 4`
 * The final value of `x[]` is `4`

The write from Thread A is overwritten, and so the value of `x[]` has increased
by `1` instead of `2`.

Similar to the earlier `ids` example, we can fix this data race using a
`ReentrantLock`. The lock ensures that only one thread can update (read and then
write) `x` at a time. Now we get the correct answer:
````julia
julia> begin
           x = Ref(0)
           l = ReentrantLock()
           Threads.@threads for i in 1:Threads.nthreads()
                for i in 1:1_000
                    lock(l) do
                        x[] += 1
                    end
                end
           end
           x[]
       end
4000
````

With the lock, the sequence of events goes something like:

 * Assume `x[]` currently has a value of `3`
 * Thread A acquires the lock
 * Thread A reads `x[]` to get `3`
 * Thread B asks for the lock, but is denied because A is currently using it
 * Thread A writes `x[] = 3 + 1 = 4`
 * Thread A releases the lock
 * Thread B acquires the lock
 * Thread B reads `x[]` to get `4`
 * Thread B writes `x[] = 4 + 1 = 5`
 * Thread B releases the lock
 * The final value of `x[]` is `5`

See the [Multi-threading](https://docs.julialang.org/en/v1/manual/multi-threading/)
section of the Julia documentation for more details.

### JuMP models are not thread-safe

An object is thread-safe if it can be modified by separate threads without
causing a data race. JuMP models are not thread-safe. Code that uses
multi-threading to simultaneously modify or optimize a single JuMP model
across threads may error, crash Julia, or silently produce incorrect results.

For example, the following incorrect use of multi-threading crashes Julia:
```julia
julia> using JuMP, HiGHS

julia> function an_incorrect_way_to_use_threading()
           model = Model(HiGHS.Optimizer)
           set_silent(model)
           @variable(model, x)
           Threads.@threads for i in 1:10
               optimize!(model)
           end
           return
       end
an_incorrect_way_to_use_threading (generic function with 1 method)

julia> an_incorrect_way_to_use_threading()
julia(76918,0x16c92f000) malloc: *** error for object 0x600003e52220: pointer being freed was not allocated
zsh: abort      julia -t 4
```

To avoid issues with thread safety, create a new instance of a JuMP model in
each iteration of the for-loop. In addition, you must avoid race conditions in
the rest of your Julia code, for example, by using a lock when pushing elements
to a shared vector.

### Thread safety and the closure capture bug

!!! danger
    This section is very important to understand. It is not specific to JuMP and
    it applies to all multithreaded Julia programs.

There is an upstream design issue in Julia ([julia#14948](https://github.com/JuliaLang/julia/issues/14948))
that can silently introduce race conditions to your code and violate thread
safety.

You can trigger this bug if you have a local variable inside the
`Threads.@threads` loop with the same name as a variable outside the loop.
Here's an example:

```julia
julia> using JuMP, HiGHS

julia> function _create_model(j)
           model = Model(HiGHS.Optimizer)
           @variable(model, x[1:j])
           return model
       end
_create_model (generic function with 1 method)

julia> function dont_run_segfault_likely()
           models = _create_model.(1:2)
           Threads.@threads for j in 1:2
               model = models[j]  # `model` is used inside the loop
               set_lower_bound.(model[:x], j)
               optimize!(model)
           end
           model = models[1]      # `model` is used outside the loop
       end
dont_run_segfault_likely (generic function with 1 method)
```

This code is problematic for the following reason. Because `model` appears
both inside `Threads.@threads` and outside, Julia's scoping rules treat it as
a single variable. Therefore, instead of creating a different `model` for each
thread, Julia creates a single mutable container called a `Core.Box` that is
used to store the value of `model` throughout the life of the function. Now
there is a race condition for reads and writes of `model`, and so a thread may
read the value of `model` only to find that its `model` was overwritten with the
value of `model` from another thread.

To diagnose this issue, use `@code_warntype`. If your code is problematic, you
will see a local variable with the type `::Core.Box`:
```julia
julia> @code_warntype dont_run_segfault_likely()
MethodInstance for dont_run_segfault_likely()
  from dont_run_segfault_likely() @ Main REPL[3]:1
Arguments
  ...
Locals
  ...
  model::Core.Box
  ...
```
If you see `Core.Box`, you must refactor your code to avoid re-using the same
variable name inside and outside the threading loop.

Alternatively, you can annotate the local variables inside the loop with `local`
to disambiguate them from the variables outside the loop.

```julia
julia> using JuMP, HiGHS

julia> function _create_model(j)
           model = Model(HiGHS.Optimizer)
           @variable(model, x[1:j])
           return model
       end
_create_model (generic function with 1 method)

julia> function safe_to_run()
           models = _create_model.(1:2)
           Threads.@threads for j in 1:2
               local model = models[j]  # annotated as `local`
               set_lower_bound.(model[:x], j)
               optimize!(model)
           end
           model = models[1]  # This `model` is not the inner `model`
       end
safe_to_run (generic function with 1 method)
```

### Example: parameter search with multi-threading

Here is an example of how to use multi-threading to solve a collection of JuMP
models in parallel.
````julia
julia> using JuMP, HiGHS

julia> function a_good_way_to_use_threading()
           solutions = Pair{Int,Float64}[]
           my_lock = Threads.ReentrantLock();
           Threads.@threads for i in 1:10
               model = Model(HiGHS.Optimizer)
               set_silent(model)
               set_attribute(model, MOI.NumberOfThreads(), 1)
               @variable(model, x >= i)
               @objective(model, Min, x)
               optimize!(model)
               assert_is_solved_and_feasible(model)
               Threads.lock(my_lock) do
                   push!(solutions, i => objective_value(model))
               end
           end
           return solutions
       end
a_good_way_to_use_threading (generic function with 1 method)

julia> a_good_way_to_use_threading()
10-element Vector{Pair{Int64, Float64}}:
  7 => 7.0
  9 => 9.0
  4 => 4.0
  1 => 1.0
  5 => 5.0
  2 => 2.0
  8 => 8.0
 10 => 10.0
  3 => 3.0
  6 => 6.0
````

!!! warning
    For some solvers, it may be necessary to limit the number of threads used
    internally by the solver to 1 by setting the [`MOI.NumberOfThreads`](@ref)
    attribute.

### Example: building data structures in parallel

For large problems, building the model in JuMP can be a bottleneck, and you
may consider trying to write code that builds the model in parallel, for
example, by wrapping a `for`-loop that adds constraints with
`Threads.@threads`. Here's an example:

```julia
julia> using JuMP

julia> function an_incorrect_way_to_build_with_multithreading()
           model = Model()
           @variable(model, x[1:10])
           Threads.@threads for i in 1:10
               @constraint(model, x[i] <= i)
           end
           return model
       end

julia> an_incorrect_way_to_build_with_multithreading()
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 10
├ num_constraints: 7
│ └ AffExpr in MOI.LessThan{Float64}: 7
└ Names registered in the model
  └ :x
```

Unfortunately, this model is wrong. It has only seven constraints instead of the
expected ten. This happens because JuMP models are not thread-safe. Code that
uses multi-threading to simultaneously modify or optimize a single JuMP model
across threads may error, crash Julia, or silently produce incorrect results.

The correct way to build a JuMP model with multi-threading is to build the
data structures in parallel, but add them to the JuMP model in a thread-safe
way:
```julia
julia> using JuMP

julia> function a_correct_way_to_build_with_multithreading()
           model = Model()
           @variable(model, x[1:10])
           my_lock = Threads.ReentrantLock()
           Threads.@threads for i in 1:10
               con = @build_constraint(x[i] <= i)
               Threads.lock(my_lock) do
                   add_constraint(model, con)
               end
           end
           return model
       end

julia> a_correct_way_to_build_with_multithreading()
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 10
├ num_constraints: 10
│ └ AffExpr in MOI.LessThan{Float64}: 10
└ Names registered in the model
  └ :x
```

!!! warning
    **Do not use multi-threading to build a JuMP model just because your original
    code is slow.** In most cases, we find that the reason for the bottleneck is
    not JuMP, but in how you are constructing the problem data, and that with
    changes, it is possible to build a model in a way that is not the bottleneck
    in the solution process. If you need help to make your code run faster, ask
    for help on the [community forum](https://jump.dev/forum). Make sure to
    include a reproducible example of your code.

## Distributed computing

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

### Example: parameter search with distributed computing

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
               assert_is_solved_and_feasible(sudoku)
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

## Parallelism within the solver

Many solvers use parallelism internally. For example, commercial solvers like
[Gurobi.jl](@ref) and [CPLEX.jl](@ref) both parallelize the search in
branch-and-bound.

Solvers supporting internal parallelism will typically support the
[`MOI.NumberOfThreads`](@ref) attribute, which you can set using
[`set_attribute`](@ref):

```julia
using JuMP, Gurobi
model = Model(Gurobi.Optimizer)
set_attribute(model, MOI.NumberOfThreads(), 4)
```

## GPU parallelism

JuMP does not support GPU programming, but some solvers support execution on a
GPU.

One example is [SCS.jl](@ref), which supports using a GPU to internally solve
a system of linear equations. If you are on `x86_64` Linux machine, do:
```julia
using JuMP, SCS, SCS_GPU_jll
model = Model(SCS.Optimizer)
set_attribute(model, "linear_solver", SCS.GpuIndirectSolver)
```
