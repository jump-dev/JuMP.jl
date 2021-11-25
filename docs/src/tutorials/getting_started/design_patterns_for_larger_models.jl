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

# # Design partterns for larger models

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
# knapsack model, and then expand it to add various features and structure.

# ## A simple script

# Your first prototype of a JuMP model is probably a script that uses a small
# set of hard-coded data.

using JuMP, GLPK
profit = [5, 3, 2, 7, 4]
weight = [2, 8, 4, 2, 5]
capacity = 10
N = 5
model = Model(GLPK.Optimizer)
@variable(model, x[1:N], Bin)
@objective(model, Max, sum(profit[i] * x[i] for i in 1:N))
@constraint(model, sum(weight[i] * x[i] for i in 1:N) <= capacity)
optimize!(model)
value.(x)

# The benefit of this approach are:
#  * it is quick to code
#  * it is quick to make changes.

# The downsides include:
#  * all variables are global (read [Performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/))
#  * it is easy to introduce errors, e.g., having `profit` and `weight` be
#    vectors of different lengths, or not match `N`
#  * the solution, `x[i]`, is hard to interpret without knowing the order in
#    which we provided the data.

# ## Wrap the model in a fuction

# A good next step is to wrap your model in a function. This is useful for a few
# reasons:
#  * it removes global variables
#  * it encapsulates the JuMP model and forces you to clarify your inputs and
#    outputs
#  * we can add some error checking.

function solve_knapsack_1(profit::Vector, weight::Vector, capacity::Real)
    @assert length(profit) == length(weight)
    N = length(weight)
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:N], Bin)
    @objective(model, Max, sum(profit[i] * x[i] for i in 1:N))
    @constraint(model, sum(weight[i] * x[i] for i in 1:N) <= capacity)
    optimize!(model)
    return value.(x)
end

solve_knapsack_1([5, 3, 2, 7, 4], [2, 8, 4, 2, 5], 10)

# ## Create better data structures

# Although we can check for errors like mis-matched vector lengths, if you start
# to develop models with a lot of data, keeping track of vectors and lengths and
# indices is fragile and a common source of bugs. A good solution is to use
# Julia's type system to create an abstraction over your data.

# For example, we can create a struct that represents a single object:

struct KnapsackObject
    profit::Float64
    weight::Float64
end

# as well as a struct that holds a dictionary of objects and the knapsack's
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
        println(io, "  $(rpad(k, 8)) : profit = $(v.profit), weight = $(v.weight)")
    end
    return
end

data

# Then, we can re-write our `solve_knapsack` function to take our `KnapsackData`
# as input:

function solve_knapsack_2(data::KnapsackData)
    model = Model(GLPK.Optimizer)
    @variable(model, x[keys(data.objects)], Bin)
    @objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
    @constraint(
        model,
        sum(v.weight * x[k] for (k, v) in data.objects) <= data.capacity,
    )
    optimize!(model)
    return value.(x)
end

solve_knapsack_2(data)

# ## Read in data from files

# Having a datastructure is a good step. But it is still annoying that we have
# to hard-code the data into Julia. A good next step is to separate the data
# into an external file format--JSON is a common choice.

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

data = read_data(joinpath(@__DIR__, "data", "knapsack.json"))

# ## Add options via if-else

# At this point, we have data in a file format which we can load and solve a
# single problem. For many users, this might be sufficient. However, at some
# point you may be asked to add features like "but what if I want to take more
# than one of a particular item?"

# If this is the first time that you've been asked to add a feature, adding
# options via `if-else` statements is a good approach. For example, we might
# write:

function solve_knapsack_3(data::KnapsackData; binary_knapsack::Bool)
    model = Model(GLPK.Optimizer)
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
    return value.(x)
end

# Now we can solve the binary knapsack:

solve_knapsack_3(data; binary_knapsack = true)

# And an integer knapsack where we can take more than one copy of each item:

solve_knapsack_3(data; binary_knapsack = false)

# ## Add options via dispatch

# If you get repeated requests to add different options, you'll quickly find
# yourself in a mess of different flags and `if-else` statements. It's hard to
# write, hard to read, and hard to ensure you haven't introduced any bugs.
# A good solution is to use Julia's type dispatch. The easiest way to explain
# this is by example.

# First, start by defining a new abstract type, as well as new subtypes for each
# of our options.

abstract type AbstractKnapsack end

struct BinaryKnapsack <: AbstractKnapsack end

struct IntegerKnapsack <: AbstractKnapsack end

# Then, we rewrite our `solve_knapsack` function to take a `config` argument,
# and we introduce an `add_knapsack_variables` function to abstract the creation
# of our variables.

function solve_knapsack_4(data::KnapsackData, config::AbstractKnapsack)
    model = Model(GLPK.Optimizer)
    x = add_knapsack_variables(model, data, config)
    @objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
    @constraint(
        model,
        sum(v.weight * x[k] for (k, v) in data.objects) <= data.capacity,
    )
    optimize!(model)
    return value.(x)
end

# For the binary knapsack problem, `add_knapsack_variables` looks like this:

function add_knapsack_variables(
    model::Model,
    data::KnapsackData,
    ::BinaryKnapsack,
)
    return @variable(model, x[keys(data.objects)], Bin)
end

# For the integer knapsack problem, `add_knapsack_variables` looks like this:

function add_knapsack_variables(
    model::Model,
    data::KnapsackData,
    ::IntegerKnapsack,
)
    return @variable(model, x[keys(data.objects)] >= 0, Int)
end

# Now we can solve the binary knapsack:

solve_knapsack_4(data, BinaryKnapsack())

# and the integer knapsack problem:

solve_knapsack_4(data, IntegerKnapsack())

# The main benefit of the dispatch approach is that you can quickly add new
# options without needing to modify the existing code. For example:

struct UpperBoundedKnapsack <: AbstractKnapsack
    limit::Int
end

function add_knapsack_variables(
    model::Model,
    data::KnapsackData,
    config::UpperBoundedKnapsack,
)
    return @variable(model, 0 <= x[keys(data.objects)] <= config.limit, Int)
end

solve_knapsack_4(data, UpperBoundedKnapsack(3))

# ## Generalize constraints and objectives

# It's easy to extend the dispatch approach to constraints and objectives as
# well. The key points to notice in the next two functions are that:
#  * we can access registered variables via `model[:x]`
#  * we can define generic functions which accept any `AbstractKnapsack` as a
#    configuration argument. That means we can implement a single method and
#    have it apply to multiple configuration types.

function add_knapsack_constraints(
    model::Model,
    data::KnapsackData,
    ::AbstractKnapsack,
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
    ::AbstractKnapsack,
)
    x = model[:x]
    @objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
    return
end

function solve_knapsack_5(data::KnapsackData, config::AbstractKnapsack)
    model = Model(GLPK.Optimizer)
    add_knapsack_variables(model, data, config)
    add_knapsack_constraints(model, data, config)
    add_knapsack_objective(model, data, config)
    optimize!(model)
    return value.(model[:x])
end

solve_knapsack_5(data, BinaryKnapsack())

# ## Remove solver dependence, add error checks

# Compared to where we started, our knapsack model is now significantly
# different. We've wrapped it in a function, defined some data types, and
# introduced configuration options to control the variables and constraints that
# get added. There's a few other steps we can do to further improve things. They
# are:
#  * remove the dependence on `GLPK`
#  * add checks that we found an optimal solution
#  * add a helper function to avoid the need to explicitly construct the data.

function solve_knapsack_6(
    optimizer,
    data::KnapsackData,
    config::AbstractKnapsack,
)
    model = Model(optimizer)
    _add_knapsack_variables(model, data, config)
    _add_knapsack_constraints(model, data, config)
    _add_knapsack_objective(model, data, config)
    optimize!(model)
    if termination_status(model) != OPTIMAL
        @warn("Model not solved to optimality")
        return nothing
    end
    return value.(model[:x])
end

function solve_knapsack_6(optimizer, data::String, config::AbstractKnapsack)
    return solve_knapsack_6(optimizer, read_data(data), config)
end

data_filename = joinpath(@__DIR__, "data", "knapsack.json")
solution = solve_knapsack_6(GLPK.Optimizer, data_filename, BinaryKnapsack())

# ## Create a module

# Now we're ready to expose our model to the wider-world. That might be as part
# of a larger Julia project that we're contributing to, or as a stand-alone
# script that we can run on-demand. In either case, it's good practice to wrap
# everything in a module. This further encapsulates our code into a single
# namespace, and we can add documentation in the form of docstrings.

# Some good rules to follow whe creating a module are:
# * Use `import` in a module instead of `using` to make it clear which functions
#   are from which packages
# * Use `_` to start function and type names that are considered private
# * Add docstrings to all public variables and functions.

module KnapsackModel

import JuMP
import JSON

struct _KnapsackObject
    profit::Float64
    weight::Float64
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

abstract type _AbstractKnapsack end

"""
    BinaryKnapsack()

Create a binary knapsack problem where each object can be taken 0 or 1 times.
"""
struct BinaryKnapsack <: _AbstractKnapsack end

"""
    IntegerKnapsack()

Create an integer knapsack problem where each object can be taken any amount of
times.
"""
struct IntegerKnapsack <: _AbstractKnapsack end

function _add_knapsack_variables(
    model::JuMP.Model,
    data::_KnapsackData,
    ::BinaryKnapsack,
)
    return JuMP.@variable(model, x[keys(data.objects)], Bin)
end

function _add_knapsack_variables(
    model::JuMP.Model,
    data::_KnapsackData,
    ::IntegerKnapsack,
)
    return JuMP.@variable(model, x[keys(data.objects)] >= 0, Int)
end

function _add_knapsack_constraints(
    model::JuMP.Model,
    data::_KnapsackData,
    ::_AbstractKnapsack,
)
    x = model[:x]
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
    ::_AbstractKnapsack,
)
    x = model[:x]
    JuMP.@objective(model, Max, sum(v.profit * x[k] for (k, v) in data.objects))
    return
end

function _solve_knapsack(
    optimizer,
    data::_KnapsackData,
    config::_AbstractKnapsack,
)
    model = JuMP.Model(optimizer)
    _add_knapsack_variables(model, data, config)
    _add_knapsack_constraints(model, data, config)
    _add_knapsack_objective(model, data, config)
    JuMP.optimize!(model)
    if JuMP.termination_status(model) != JuMP.OPTIMAL
        @warn("Model not solved to optimality")
        return nothing
    end
    return JuMP.value.(model[:x])
end

"""
    solve_knapsack(optimizer, data_filename::String, config::_AbstractKnapsack)

Solve the knapsack problem and return the optimal primal solution

## Arguments

 * `optimizer` : an object that can be passed to `JuMP.Model` to construct a new
   JuMP model.
 * `data_filename` : the filename of a JSON file containing the data for the
   problem.
 * `config` : an object to control the type of knapsack model constructed.
   Valid options are:
    * `BinaryKnapsack()`
    * `IntegerKnapsack()`

## Returns

 * If an optimal solution exists: a `JuMP.DenseAxisArray` that maps the `String`
   name of each object to the number of objects to pack into the knapsack.
 * Otherwise, `nothing`, indicating that the problem does not have an optimal
   solution.

## Examples

```julia
solution = solve_knapsack(GLPK.Optimizer, "path/to/data.json", BinaryKnapsack())
```

```julia
solution = solve_knapsack(
    MOI.OptimizerWithAttributes(GLPK.Optimizer, "msg_lev" => 0),
    "path/to/data.json",
    IntegerKnapsack(),
)
```
"""
function solve_knapsack(
    optimizer,
    data_filename::String,
    config::_AbstractKnapsack,
)
    return _solve_knapsack(optimizer, _read_data(data_filename), config)
end

end

# Finally, you can call your model:

import .KnapsackModel

KnapsackModel.solve_knapsack(
    GLPK.Optimizer,
    joinpath(@__DIR__, "data", "knapsack.json"),
    KnapsackModel.BinaryKnapsack(),
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
            GLPK.Optimizer,
            joinpath(@__DIR__, "data", "knapsack.json"),
            KnapsackModel.BinaryKnapsack(),
        )
        @test isapprox(x["apple"], 1, atol=1e-5)
        @test isapprox(x["banana"], 0, atol=1e-5)
        @test isapprox(x["cherry"], 0, atol=1e-5)
        @test isapprox(x["date"], 1, atol=1e-5)
        @test isapprox(x["eggplant"], 1, atol=1e-5)
    end
    @testset "feasible_integer_knapsack" begin
        x = KnapsackModel.solve_knapsack(
            GLPK.Optimizer,
            joinpath(@__DIR__, "data", "knapsack.json"),
            KnapsackModel.IntegerKnapsack(),
        )
        @test isapprox(x["apple"], 0, atol=1e-5)
        @test isapprox(x["banana"], 0, atol=1e-5)
        @test isapprox(x["cherry"], 0, atol=1e-5)
        @test isapprox(x["date"], 5, atol=1e-5)
        @test isapprox(x["eggplant"], 0, atol=1e-5)
    end
    @testset "feasible_binary_knapsack" begin
        x = KnapsackModel.solve_knapsack(
            GLPK.Optimizer,
            joinpath(@__DIR__, "data", "knapsack_infeasible.json"),
            KnapsackModel.BinaryKnapsack(),
        )
        @test x === nothing
    end
end

# !!! tip
#     Place these tests in a separate file `test_knapsack_model.jl` so that you
#     can run the tests by `include`ing the file.
