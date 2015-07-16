# Introduction to modeling with JuMP

In this chapter we will present the basics of modeling linear and quadratic optimization problems with JuMP. We will assume you are familiar with the language of constrained optimization, but may not have used another modeling language like AMPL, PuLP, or similar before.

## Loading JuMP

To begin using JuMP in a script or IJulia notebook, you will need to install it the first time you use it:

```julia
Pkg.add("JuMP")
```

and then load it (every time you use it):

```julia
using JuMP
```

This makes various JuMP functions and types available to use. New features and bug fixes are added to JuMP periodically, so to make sure you have the latest stable version you should run `Pkg.update()` periodically.


## The `Model`

At the core of a JuMP optimization model is the `Model` type. The following code constructs a new empty model:

```julia
m = Model()
```

The variable name `m` can be whatever you like. You'll be using it a lot, so many users prefer shorter names. Just be careful you don't accidently use `m` for something else, like `for m = 1:10`, later on.

By default JuMP will automagically load and use one of the solvers you have installed based on the model you construct. If you'd like to set the solver yourself, or set options of the solver, you should pass a `solver` keyword argument to the `Model` constructor:

```julia
using Gurobi
milp_model = Model(solver=GurobiSolver(OutputFlag=0))
using Ipopt
nlp_model = Model(solver=IpoptSolver(file_print_level=3))
```

There are two things to take note of here. First, we need to load the solver's 