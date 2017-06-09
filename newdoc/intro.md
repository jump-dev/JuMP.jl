# Introduction to modeling with JuMP

In this chapter we will present the basics of modeling linear and quadratic optimization problems with JuMP. We will assume you are familiar with the language of constrained optimization, but may not have used another modeling language like AMPL, PuLP, or similar before.

!!! note "Note:"
    Nonlinear modeling has some differences. We'll discuss those later in the Nonlinear Modeling chapter.

## Loading JuMP

To begin using JuMP in a script or IJulia notebook, you will need to install it the first time you use it:

    Pkg.add("JuMP")

and then load it (every time you use it):

    using JuMP

This makes JuMP's functionality available to use in your program.
New features and bug fixes are added to JuMP periodically, so to make sure you have the latest stable version you should run `Pkg.update()` periodically.


## The Model

At the core of a JuMP optimization model is the `Model` type. A model describes an optimization problem: the decision variables, the constraints, and the objective. The following code constructs a new empty model:

    m = Model()

The name `m` can be any valid Julia variable name, like `portfolio`, `mymodel`, `contprob`, etc. You'll be using it a lot, so many users prefer shorter names.

!!! warning "Warning:"
    Be careful you don't accidently use `m` for something else later on, like `for m = 1:10`, later on. Also beware of `mod`: that is the name of the modulo division function in Julia, so you might run into trouble later on.

By default JuMP will *automagically* load and use one of the solvers you have installed based on the model you construct. If you'd like to set the solver yourself, or set options of the solver, you should pass a `solver` keyword argument to the `Model` constructor:

    using Gurobi
    milp_model = Model(solver=GurobiSolver(OutputFlag=0))
    using Ipopt
    nlp_model = Model(solver=IpoptSolver(file_print_level=3))

There are two things to take note of here.
First, we need to load the respective solver's package (e.g. `Gurobi.jl`) so the solver object (e.g. `GurobiSolver`) is available.
Second, the options to the solver are passed as keyword arguments to the solver object. The options are solver-specific: to find out what they are, and the values they can take, refer to the resepctive solver's documentation. If you'd like to change the solver later, you can do that with `setSolver`, e.g.:

    using Cbc
    setSolver(milp_solver, CbcSolver())

!!! note "Note:"

    JuMP supports many solvers. The full list of solvers is available at the [JuliaOpt website](http://juliaopt.org).


## Variables

Now we have an empty model, we can build up the problem by adding variables. We need to distinguish between two different senses of the word "variable":

* *Decision variables* are the mathematical variables of our problem. In a linear optimization context, they are sometimes referred to as columns.
* *Julia variables* are variables in the programming sense. They are a "shorthand" for values.

In JuMP, we add decision variables to our model and JuMP gives us a Julia variable by which we can refer to these decision variables. Decision variables are internally defined by a simple index that starts at 1 and increases; JuMP provides an abstraction over the decision variables that lets us refer to these variables in a more intuitive and maintainable way. 
For a more comprehensive description of the ways to create variables, including column generation, warm starts, and conditions, please see the Variables chapter.

### Single variables

To create a decision variable we use the JuMP macro `@defVar`. The simplest case is just a single variable, for example:

    m = Model()
    @defVar(m, 0 <= stock_order <= 100)
    @defVar(m, num_elephants >= 0, Int)
    @defVar(m, open_facility, Bin)
    @defVar(m, final_cost <= 9999)

Each of these lines creates a single *decision variable* that we can then refer to by the name we gave it (which is a *Julia variable*). We can set variable lower and upper bounds, as well as variable types like integer `Int` and `Bin`.

!!! note "Note:"
    If you wish to set just a lower bound, use `x >= 0`. The syntax `0 <= x` is not valid.

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

### Groups of variables

It is common to not just have single variables, but groups of variables. We can define these with `@defVar` as well, for example:

    m = Model()
    @defVar(m, boringvariable[1:9,1:9,1:9])
    @defVar(m, 0 <= pixel_intensity[1:255,1:255] <= 1)
    @defVar(m, bit_pattern[0:8:63], Bin)
    N = 5, M = 10
    @defVar(m, trucks_dispatched[i=1:N,j=1:M] >= 0, Int)
    items = [:sock,:sandal,:boot]
    max_stock = [:sock => 10, :sandal => 13, :boot => 5]
    @defVar(m, 0 <= stock_levels[item=items] <= max_stock[item])

All these commands make a single *Julia variable* that refers to multiple *decision variables*. We can refer to the individual decision variables using indexing, for example:

    pixel_intensity[2,10]
    stock_levels[:boot,12]

We can also use complex indexing schemes: in fact, we can index on almost any iterable set. This might be a simple range like `1:9`, more complicated ranges like `-10:2:10` (step size of 2), or general vectors of items like `[1,sqrt(2),4]` and `["my","index"]`.

You don't need to, but if you'd like you can name the indices, e.g. `i=1:N` and `item=items`. If you do name the indices, you can use those to give bounds on a per-variable basis. In the `stock_levels` example above we used the `item` index name to set a different maximum for each shoe type.


## Objective and Constraints

Now we have our model and our decision variables, we can add an objective function and constraints. Setting the objective function is done with `@setObjective`:

    @defVar(m, x)
    @defVar(m, y[1:10])
    @defVar(m, z[1:10])
    @setObjective(m, Max, x + y[1] + y[2] + y[3])
    @setObjective(m, Min, sum(z) - sum(y))
    @setObjective(m, Min, dot(y,z))

The first argument to `@setObjective` is the model (in this case, `m`). The second argument is the *sense*: JuMP supports `Min` (for minimization) and `Max` (for maximization). The third argument is an *expression*. JuMP expressions are combinations of numbers, data, and decision variables. The following are all valid expressions:

    # Linear expressions
    x - z[1]
    4x + 5y[1] + (x + 2*z[3])*5
    sin(123)*x + abs(mydata)*z[2]
    dot(rand(10), z)

    # Quadratic expressions
    x*x + (y-2)^2

!!! note "Note:"
    Some nonlinear expressions are valid too, but we will defer that to the Nonlinear Modeling chapter.

### Advanced sums

In the examples above that we used both the Julia `sum` and `dot` functions. But sometimes we need more complicated sums of our variables: for example, fixing one dimension and summing over another, or doing calculations for each index.