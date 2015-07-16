# Introduction to modeling with JuMP

In this chapter we will present the basics of modeling linear and quadratic optimization problems with JuMP. We will assume you are familiar with the language of constrained optimization, but may not have used another modeling language like AMPL, PuLP, or similar before.

## Loading JuMP

To begin using JuMP in a script or IJulia notebook, you will need to install it the first time you use it:

    Pkg.add("JuMP")

and then load it (every time you use it):

    using JuMP

This makes various JuMP functions and types available to use. New features and bug fixes are added to JuMP periodically, so to make sure you have the latest stable version you should run `Pkg.update()` periodically.


## The `Model`

At the core of a JuMP optimization model is the `Model` type. The following code constructs a new empty model:

    m = Model()

The variable name `m` can be whatever you like. You'll be using it a lot, so many users prefer shorter names. Just be careful you don't accidently use `m` for something else, like `for m = 1:10`, later on.

By default JuMP will automagically load and use one of the solvers you have installed based on the model you construct. If you'd like to set the solver yourself, or set options of the solver, you should pass a `solver` keyword argument to the `Model` constructor:

    using Gurobi
    milp_model = Model(solver=GurobiSolver(OutputFlag=0))
    using Ipopt
    nlp_model = Model(solver=IpoptSolver(file_print_level=3))

There are two things to take note of here.
First, we need to load the respective solver's package (e.g. `Gurobi.jl`) so the solver object (e.g. `GurobiSolver`) is available.
Second, the options to the solver are passed as keyword arguments to the solver object. The options are solver-specific: to find out what they are, and the values they can take, refer to the resepctive solver's documentation. If you'd like to change the solver later, you can do that with `setSolver`, e.g.:

    using Cbc
    setSolver(milp_solver, CbcSolver())


## Variables

Now we have an empty model, we can build up the problem by adding variables. We need to distinguish between two different senses of the word "variable":

* *Decision variables* are the mathematical variables of our problem. In a linear optimization context, they are sometimes referred to as columns.
* *Julia variables* are variables in the programming sense. They are a "shorthand" for values.

In JuMP, we add decision variables to our model and JuMP gives us a Julia variable by which we can refer to these decision variables. Decision variables are internally defined by a simple index that starts at 1 and increases; JuMP provides an abstraction over the decision variables that lets us refer to these variables in a more intuitive and maintainable way.

### Single variables

To create a decision variable we use the JuMP macro `@defVar`. The simplest case is just a single variable, for example:

    m = Model()
    @defVar(m, 0 <= stock_order <= 100)
    @defVar(m, num_elephants >= 0, Int)
    @defVar(m, open_facility, Bin)
    @defVar(m, final_cost <= 9999)

Each of these lines creates a single *decision variable* that we can then refer to by the name we gave it (which is a *Julia variable*). We can set variable lower and upper bounds, as well as variable types like integer `Int` and `Bin`.

!!! warning "Warning:"
    If you create a new decision variable with the same Julia variable name, you will no longer be able to refer to that decision variable. For example, the following code defines two decision variables but `x` refers only to the second:

        m = Model()
        @defVar(m, x >= 0)  # This gets 'lost'
        @defVar(m, x >= 0)

    If you want to keep using the same name for some reason, but still need to keep a reference to the original variable, consider storing the variables in a list, e.g.

        m = Model()
        @defVar(m, x >= 0)
        vars = [x]
        @defVar(m, x >= 0)
        push!(vars, x)
        @show length(x)  # length(x) => 2

For a comprehensive description of the arguments to `@defVar`, see the Variables chapter.

### Groups of variables

It is common to not just have single variables, but groups of variables. We can define these with `@defVar` as well, for example:

    m = Model()
    items = [:sock,:sandal,:boot,:highheel]
    @defVar(m, stock_levels[items,size=1:16] >= 0)
    @defVar(m, 0 <= pixel_intensity[1:255,1:255] <= 1)
    @defVar(m, bit_pattern[0:8:63], Bin)
    @defVar(m, trucks_dispatched[i=1:N,j=1:N] >= 0, Int)

All these commands make a single Julia variable that refers to multiple decision variables. We can refer to the individual decision variables using indexing:

    pixel_intensity[2,10]
    stock_levels[:boot,12]

There are lot of interesting things going on here:

  * We can use an iterable set to index on. This might be a simple range like `1:10`, more complicated ranges like `-10:2:10`, or general vectors like `[1,4,sqrt(2)]` and `["my","index"]`