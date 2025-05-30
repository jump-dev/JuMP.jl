# Copyright (c) 2021 Oscar Dowson and contributors                               #src
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

# # Design patterns for larger models

# JuMP makes it easy to build and solve optimization models. However, once you
# start to construct larger models, and especially ones that interact with
# external data sources or have customizable sets of variables and constraints
# based on client choices, you may find that your scripts become unwieldy. This
# tutorial demonstrates a variety of ways in which you can structure larger JuMP
# models to improve their readability and maintainability.

# !!! tip
#     This tutorial is more advanced than the other "Getting started" tutorials.
#     It's in the "Getting started" section to give you an early preview of how
#     JuMP makes it easy to structure larger models. However, if you are new to
#     JuMP you may want to briefly skim the tutorial, and come back to it once
#     you have written a few JuMP models.

# ## Overview

# This tutorial uses explanation-by-example. We're going to start with a simple
# [knapsack model](https://en.wikipedia.org/wiki/Knapsack_problem), and then
# expand it to add various features and structure.

# ## A simple script

# Your first prototype of a JuMP model is probably a script that uses a small
# set of hard-coded data.

using JuMP, HiGHS
profit = [5, 3, 2, 7, 4]
weight = [2, 8, 4, 2, 5]
capacity = 10
N = 5
model = Model(HiGHS.Optimizer)
@variable(model, x[1:N], Bin)
@objective(model, Max, sum(profit[i] * x[i] for i in 1:N))
@constraint(model, sum(weight[i] * x[i] for i in 1:N) <= capacity)
optimize!(model)
assert_is_solved_and_feasible(model)
value.(x)

# The benefits of this approach are:
#  * it is quick to code
#  * it is quick to make changes.

# The downsides include:
#  * all variables are global (read [Performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/))
#  * it is easy to introduce errors, for example, having `profit` and `weight` be
#    vectors of different lengths, or not match `N`
#  * the solution, `x[i]`, is hard to interpret without knowing the order in
#    which we provided the data.

# ## Wrap the model in a function

# A good next step is to wrap your model in a function. This is useful for a few
# reasons:
#  * it removes global variables
#  * it encapsulates the JuMP model and forces you to clarify your inputs and
#    outputs
#  * we can add some error checking.

function solve_knapsack_1(profit::Vector, weight::Vector, capacity::Real)
    if length(profit) != length(weight)
        throw(DimensionMismatch("profit and weight are different sizes"))
    end
    N = length(weight)
    model = Model(HiGHS.Optimizer)
    @variable(model, x[1:N], Bin)
    @objective(model, Max, sum(profit[i] * x[i] for i in 1:N))
    @constraint(model, sum(weight[i] * x[i] for i in 1:N) <= capacity)
    optimize!(model)
    assert_is_solved_and_feasible(model)
    return value.(x)
end

solve_knapsack_1([5, 3, 2, 7, 4], [2, 8, 4, 2, 5], 10)

# ## Create better data structures

# Although we can check for errors like mis-matched vector lengths, if you start
# to develop models with a lot of data, keeping track of vectors and lengths and
# indices is fragile and a common source of bugs. A good solution is to use
# Julia's type system to create an abstraction over your data.

# For example, we can create a `struct` that represents a single object, with a
# constructor that lets us validate assumptions on the input data:

struct KnapsackObject
    profit::Float64
    weight::Float64
    function KnapsackObject(profit::Float64, weight::Float64)
        if weight < 0
            throw(DomainError("Weight of object cannot be negative"))
        end
        return new(profit, weight)
    end
end

# as well as a `struct` that holds a dictionary of objects and the knapsack's
# capacity:

struct KnapsackData
    objects::Dict{String,KnapsackObject}
    capacity::Float64
end

# Here's what our data might look like now:

objects = Dict(
    "apple" => KnapsackObject(5.0, 2.0),
    "banana" => KnapsackObject(3.0, 8.0),
    "cherry" => KnapsackObject(2.0, 4.0),
    "date" => KnapsackObject(7.0, 2.0),
    "eggplant" => KnapsackObject(4.0, 5.0),
)
data = KnapsackData(objects, 10.0)

# If you want, you can add custom printing to make it easier to visualize:

function Base.show(io::IO, data::KnapsackData)
    println(io, "A knapsack with capacity $(data.capacity) and possible items:")
    for (k, v) in data.objects
        println(
            io,
            "  $(rpad(k, 8)) : profit = $(v.profit), weight = $(v.weight)",
        )
    end
    return
end

data

# Then, we can re-write our `solve_knapsack` function to take our `KnapsackData`
# as input:

function solve_knapsack_2(data::KnapsackData)
    model = Model(HiGHS.Optimizer)
    @variable(model, x[keys(data.objects)], Bin)
    @objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
    @constraint(
        model,
        sum(v.weight * x[k] for (k, v) in data.objects) <= data.capacity,
    )
    optimize!(model)
    assert_is_solved_and_feasible(model)
    return value.(x)
end

solve_knapsack_2(data)

# ## Read in data from files

# Having a data structure is a good step. But it is still annoying that we have
# to hard-code the data into Julia. A good next step is to separate the data
# into an external file format; JSON is a common choice.

json_data = """
{
    "objects": {
        "apple": {"profit": 5.0, "weight": 2.0},
        "banana": {"profit": 3.0, "weight": 8.0},
        "cherry": {"profit": 2.0, "weight": 4.0},
        "date": {"profit": 7.0, "weight": 2.0},
        "eggplant": {"profit": 4.0, "weight": 5.0}
    },
    "capacity": 10.0
}
"""
temp_dir = mktempdir()
knapsack_json_filename = joinpath(temp_dir, "knapsack.json")
## Instead of writing a new file here you could replace `knapsack_json_filename`
## with the path to a local file.
write(knapsack_json_filename, json_data);

# Now let's write a function that reads this file and builds a `KnapsackData`
# object:

import JSON

function read_data(filename)
    d = JSON.parsefile(filename)
    return KnapsackData(
        Dict(
            k => KnapsackObject(v["profit"], v["weight"]) for
            (k, v) in d["objects"]
        ),
        d["capacity"],
    )
end

data = read_data(knapsack_json_filename)

# ## Add options via if-else

# At this point, we have data in a file format which we can load and solve a
# single problem. For many users, this might be sufficient. However, at some
# point you may be asked to add features like "but what if we want to take more
# than one of a particular item?"

# If this is the first time that you've been asked to add a feature, adding
# options via `if-else` statements is a good approach. For example, we might
# write:

function solve_knapsack_3(data::KnapsackData; binary_knapsack::Bool)
    model = Model(HiGHS.Optimizer)
    if binary_knapsack
        @variable(model, x[keys(data.objects)], Bin)
    else
        @variable(model, x[keys(data.objects)] >= 0, Int)
    end
    @objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
    @constraint(
        model,
        sum(v.weight * x[k] for (k, v) in data.objects) <= data.capacity,
    )
    optimize!(model)
    assert_is_solved_and_feasible(model)
    return value.(x)
end

# Now we can solve the binary knapsack:

solve_knapsack_3(data; binary_knapsack = true)

# And an integer knapsack where we can take more than one copy of each item:

solve_knapsack_3(data; binary_knapsack = false)

# ## Add configuration options via dispatch

# If you get repeated requests to add different options, you'll quickly find
# yourself in a mess of different flags and `if-else` statements. It's hard to
# write, hard to read, and hard to ensure you haven't introduced any bugs.
# A good solution is to use Julia's type dispatch to control the configuration
# of the model. The easiest way to explain this is by example.

# First, start by defining a new abstract type, as well as new subtypes for each
# of our options. These types are going to control the configuration of the
# knapsack model.

abstract type AbstractConfiguration end

struct BinaryKnapsackConfig <: AbstractConfiguration end

struct IntegerKnapsackConfig <: AbstractConfiguration end

# Then, we rewrite our `solve_knapsack` function to take a `config` argument,
# and we introduce an `add_knapsack_variables` function to abstract the creation
# of our variables.

function solve_knapsack_4(data::KnapsackData, config::AbstractConfiguration)
    model = Model(HiGHS.Optimizer)
    x = add_knapsack_variables(model, data, config)
    @objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
    @constraint(
        model,
        sum(v.weight * x[k] for (k, v) in data.objects) <= data.capacity,
    )
    optimize!(model)
    assert_is_solved_and_feasible(model)
    return value.(x)
end

# For the binary knapsack problem, `add_knapsack_variables` looks like this:

function add_knapsack_variables(
    model::Model,
    data::KnapsackData,
    ::BinaryKnapsackConfig,
)
    return @variable(model, x[keys(data.objects)], Bin)
end

# For the integer knapsack problem, `add_knapsack_variables` looks like this:

function add_knapsack_variables(
    model::Model,
    data::KnapsackData,
    ::IntegerKnapsackConfig,
)
    return @variable(model, x[keys(data.objects)] >= 0, Int)
end

# Now we can solve the binary knapsack:

solve_knapsack_4(data, BinaryKnapsackConfig())

# and the integer knapsack problem:

solve_knapsack_4(data, IntegerKnapsackConfig())

# The main benefit of the dispatch approach is that you can quickly add new
# options without needing to modify the existing code. For example:

struct UpperBoundedKnapsackConfig <: AbstractConfiguration
    limit::Int
end

function add_knapsack_variables(
    model::Model,
    data::KnapsackData,
    config::UpperBoundedKnapsackConfig,
)
    return @variable(model, 0 <= x[keys(data.objects)] <= config.limit, Int)
end

solve_knapsack_4(data, UpperBoundedKnapsackConfig(3))

# ## Generalize constraints and objectives

# It's easy to extend the dispatch approach to constraints and objectives as
# well. The key points to notice in the next two functions are that:
#  * we can access registered variables via `model[:x]`
#  * we can define generic functions which accept any `AbstractConfiguration` as a
#    configuration argument. That means we can implement a single method and
#    have it apply to multiple configuration types.

function add_knapsack_constraints(
    model::Model,
    data::KnapsackData,
    ::AbstractConfiguration,
)
    x = model[:x]
    @constraint(
        model,
        capacity_constraint,
        sum(v.weight * x[k] for (k, v) in data.objects) <= data.capacity,
    )
    return
end

function add_knapsack_objective(
    model::Model,
    data::KnapsackData,
    ::AbstractConfiguration,
)
    x = model[:x]
    @objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
    return
end

function solve_knapsack_5(data::KnapsackData, config::AbstractConfiguration)
    model = Model(HiGHS.Optimizer)
    add_knapsack_variables(model, data, config)
    add_knapsack_constraints(model, data, config)
    add_knapsack_objective(model, data, config)
    optimize!(model)
    assert_is_solved_and_feasible(model)
    return value.(model[:x])
end

solve_knapsack_5(data, BinaryKnapsackConfig())

# ## Function barriers

# Querying a variable like `x = model[:x]` is not type stable. The lack of type
# stability means that Julia cannot statically prove what the type of `x` will
# be within a function. As a consequence, the functions
# `add_knapsack_constraints` and `add_knapsack_objective` are slower than they
# would be if Julia could prove that `x` was a [`Containers.DenseAxisArray`](@ref).

# One solution to the problem of type stability is to use a function barrier:

function add_knapsack_constraints(
    model::Model,
    data::KnapsackData,
    config::AbstractConfiguration,
)
    return add_knapsack_constraints_inner(model, data, config, model[:x])
end

function add_knapsack_constraints_inner(
    model::Model,
    data::KnapsackData,
    ::AbstractConfiguration,
    x,
)
    @constraint(
        model,
        capacity_constraint,
        sum(v.weight * x[k] for (k, v) in data.objects) <= data.capacity,
    )
    return
end

# Now, `add_knapsack_constraints_inner` is faster because Julia can compile a
# specific version that depends on the type of `x` passed as an input argument.

# The downside of function barriers is that they require more code. The upside is
# faster code and reduced memory allocations. In general, you should add a
# function barrier if you have a benchmark showing that the difference is
# meaningful for your code base. Function barriers are more likely to be useful
# if the `_inner` function does a lot of computational work, for example, it
# adds thousands of constraints or has a summation over thousands of items.

# ## Remove solver dependence, add error checks

# Compared to where we started, our knapsack model is now significantly
# different. We've wrapped it in a function, defined some data types, and
# introduced configuration options to control the variables and constraints that
# get added. There are a few other steps we can do to further improve things:
#  * remove the dependence on `HiGHS`
#  * add checks that we found an optimal solution
#  * add a helper function to avoid the need to explicitly construct the data.

function solve_knapsack_6(
    optimizer,
    data::KnapsackData,
    config::AbstractConfiguration,
)
    model = Model(optimizer)
    add_knapsack_variables(model, data, config)
    add_knapsack_constraints(model, data, config)
    add_knapsack_objective(model, data, config)
    optimize!(model)
    if !is_solved_and_feasible(model)
        @warn("Model not solved to optimality")
        return nothing
    end
    return value.(model[:x])
end

function solve_knapsack_6(
    optimizer,
    data::String,
    config::AbstractConfiguration,
)
    return solve_knapsack_6(optimizer, read_data(data), config)
end

solution = solve_knapsack_6(
    HiGHS.Optimizer,
    knapsack_json_filename,
    BinaryKnapsackConfig(),
)

# ## Create a module

# Now we're ready to expose our model to the wider world. That might be as part
# of a larger Julia project that we're contributing to, or as a stand-alone
# script that we can run on-demand. In either case, it's good practice to wrap
# everything in a module. This further encapsulates our code into a single
# namespace, and we can add documentation in the form of
# [docstrings](https://docs.julialang.org/en/v1/manual/documentation/).

# Some good rules to follow when creating a module are:
# * use `import` in a module instead of `using` to make it clear which functions
#   are from which packages
# * use `_` to start function and type names that are considered private
# * add docstrings to all public variables and functions.

module KnapsackModel

import JuMP
import JSON

struct _KnapsackObject
    profit::Float64
    weight::Float64
    function _KnapsackObject(profit::Float64, weight::Float64)
        if weight < 0
            throw(DomainError("Weight of object cannot be negative"))
        end
        return new(profit, weight)
    end
end

struct _KnapsackData
    objects::Dict{String,_KnapsackObject}
    capacity::Float64
end

function _read_data(filename)
    d = JSON.parsefile(filename)
    return _KnapsackData(
        Dict(
            k => _KnapsackObject(v["profit"], v["weight"]) for
            (k, v) in d["objects"]
        ),
        d["capacity"],
    )
end

abstract type _AbstractConfiguration end

"""
    BinaryKnapsackConfig()

Create a binary knapsack problem where each object can be taken 0 or 1 times.
"""
struct BinaryKnapsackConfig <: _AbstractConfiguration end

"""
    IntegerKnapsackConfig()

Create an integer knapsack problem where each object can be taken any number of
times.
"""
struct IntegerKnapsackConfig <: _AbstractConfiguration end

function _add_knapsack_variables(
    model::JuMP.Model,
    data::_KnapsackData,
    ::BinaryKnapsackConfig,
)
    return JuMP.@variable(model, x[keys(data.objects)], Bin)
end

function _add_knapsack_variables(
    model::JuMP.Model,
    data::_KnapsackData,
    ::IntegerKnapsackConfig,
)
    return JuMP.@variable(model, x[keys(data.objects)] >= 0, Int)
end

function _add_knapsack_constraints(
    model::JuMP.Model,
    data::_KnapsackData,
    config::_AbstractConfiguration,
)
    return _add_knapsack_constraints_inner(model, data, config, model[:x])
end

function _add_knapsack_constraints_inner(
    model::JuMP.Model,
    data::_KnapsackData,
    ::_AbstractConfiguration,
    x,
)
    JuMP.@constraint(
        model,
        capacity_constraint,
        sum(v.weight * x[k] for (k, v) in data.objects) <= data.capacity,
    )
    return
end

function _add_knapsack_objective(
    model::JuMP.Model,
    data::_KnapsackData,
    config::_AbstractConfiguration,
)
    return _add_knapsack_objective_inner(model, data, config, model[:x])
end

function _add_knapsack_objective_inner(
    model::JuMP.Model,
    data::_KnapsackData,
    ::_AbstractConfiguration,
    x,
)
    JuMP.@objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
    return
end

function _solve_knapsack(
    optimizer,
    data::_KnapsackData,
    config::_AbstractConfiguration,
)
    model = JuMP.Model(optimizer)
    _add_knapsack_variables(model, data, config)
    _add_knapsack_constraints(model, data, config)
    _add_knapsack_objective(model, data, config)
    JuMP.optimize!(model)
    if !JuMP.is_solved_and_feasible(model)
        @warn("Model not solved to optimality")
        return nothing
    end
    return JuMP.value.(model[:x])
end

"""
    solve_knapsack(
        optimizer,
        knapsack_json_filename::String,
        config::_AbstractConfiguration,
    )

Solve the knapsack problem and return the optimal primal solution

## Arguments

 * `optimizer` : an object that can be passed to `JuMP.Model` to construct a new
   JuMP model.
 * `knapsack_json_filename` : the filename of a JSON file containing the data for the
   problem.
 * `config` : an object to control the type of knapsack model constructed.
   Valid options are:
    * `BinaryKnapsackConfig()`
    * `IntegerKnapsackConfig()`

## Returns

 * If an optimal solution exists: a `JuMP.Containers.DenseAxisArray` that maps
   the `String` name of each object to the number of objects to pack into the
   knapsack.
 * Otherwise, `nothing`, indicating that the problem does not have an optimal
   solution.

## Example

```julia
solution = solve_knapsack(
    HiGHS.Optimizer,
    "path/to/data.json",
    BinaryKnapsackConfig(),
)
```

```julia
solution = solve_knapsack(
    MOI.OptimizerWithAttributes(HiGHS.Optimizer, "output_flag" => false),
    "path/to/data.json",
    IntegerKnapsackConfig(),
)
```
"""
function solve_knapsack(
    optimizer,
    knapsack_json_filename::String,
    config::_AbstractConfiguration,
)
    data = _read_data(knapsack_json_filename)
    return _solve_knapsack(optimizer, data, config)
end

end

# Finally, you can call your model:

import .KnapsackModel

KnapsackModel.solve_knapsack(
    HiGHS.Optimizer,
    knapsack_json_filename,
    KnapsackModel.BinaryKnapsackConfig(),
)

# !!! note
#     The `.` in `.KnapsackModel` denotes that it is a submodule and not a
#     separate package that we installed with `Pkg.add`. If you put the
#     `KnapsackModel` in a separate file, load it with:
#     ```julia
#     include("path/to/KnapsackModel.jl")
#     import .KnapsackModel
#     ```

# ## Add tests

# As a final step, you should add tests for your model. This often means testing
# on a small problem for which you can work out the optimal solution by hand.
# The Julia standard library `Test` has good unit-testing functionality.

import .KnapsackModel
using Test

@testset "KnapsackModel" begin
    @testset "feasible_binary_knapsack" begin
        x = KnapsackModel.solve_knapsack(
            HiGHS.Optimizer,
            knapsack_json_filename,
            KnapsackModel.BinaryKnapsackConfig(),
        )
        @test isapprox(x["apple"], 1, atol = 1e-5)
        @test isapprox(x["banana"], 0, atol = 1e-5)
        @test isapprox(x["cherry"], 0, atol = 1e-5)
        @test isapprox(x["date"], 1, atol = 1e-5)
        @test isapprox(x["eggplant"], 1, atol = 1e-5)
    end
    @testset "feasible_integer_knapsack" begin
        x = KnapsackModel.solve_knapsack(
            HiGHS.Optimizer,
            knapsack_json_filename,
            KnapsackModel.IntegerKnapsackConfig(),
        )
        @test isapprox(x["apple"], 0, atol = 1e-5)
        @test isapprox(x["banana"], 0, atol = 1e-5)
        @test isapprox(x["cherry"], 0, atol = 1e-5)
        @test isapprox(x["date"], 5, atol = 1e-5)
        @test isapprox(x["eggplant"], 0, atol = 1e-5)
    end
    @testset "infeasible_binary_knapsack" begin
        dir = mktempdir()
        infeasible_filename = joinpath(dir, "infeasible.json")
        write(
            infeasible_filename,
            """{
                "objects": {
                    "apple": {"profit": 5.0, "weight": 2.0},
                    "banana": {"profit": 3.0, "weight": 8.0},
                    "cherry": {"profit": 2.0, "weight": 4.0},
                    "date": {"profit": 7.0, "weight": 2.0},
                    "eggplant": {"profit": 4.0, "weight": 5.0}
                },
                "capacity": -10.0
            }""",
        )
        x = KnapsackModel.solve_knapsack(
            HiGHS.Optimizer,
            infeasible_filename,
            KnapsackModel.BinaryKnapsackConfig(),
        )
        @test x === nothing
    end
end

# !!! tip
#     Place these tests in a separate file `test_knapsack_model.jl` so that you
#     can run the tests by adding `include("test_knapsack_model.jl")` to any
#     file where needed.

# ## Next steps

# We've only briefly scratched the surface of ways to create and structure large
# JuMP models, so consider this tutorial a starting point, rather than a
# comprehensive list of all the possible ways to structure JuMP models.  If you
# are embarking on a large project that uses JuMP, a good next step is to
# look at ways people have written large JuMP projects "in the wild."

# Here are some good examples (all co-incidentally related to energy):
# * AnyMOD.jl
#   * [JuMP-dev 2021 talk](https://www.youtube.com/watch?v=QE_tNDER0F4)
#   * [source code](https://github.com/leonardgoeke/AnyMOD.jl)
# * PowerModels.jl
#   * [JuMP-dev 2021 talk](https://www.youtube.com/watch?v=POOt1FCA8LI)
#   * [source code](https://github.com/lanl-ansi/PowerModels.jl)
# * PowerSimulations.jl
#    * [JuliaCon 2021 talk](https://www.youtube.com/watch?v=-ZoO3npjwYU)
#    * [source code](https://github.com/NREL-SIIP/PowerSimulations.jl)
# * UnitCommitment.jl
#   * [JuMP-dev 2021 talk](https://www.youtube.com/watch?v=rYUZK9kYeIY)
#   * [source code](https://github.com/ANL-CEEESA/UnitCommitment.jl)
