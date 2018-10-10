var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#![JuMP-logo](assets/jump-logo-with-text.svg)-1",
    "page": "Introduction",
    "title": "(Image: JuMP logo)",
    "category": "section",
    "text": "(Image: Powered by NumFOCUS)# These comments do not display in the HTML output.\n# See https://github.com/JuliaDocs/Documenter.jl/issues/674.\n\n# Style conventions for the JuMP documentation:\n# - Respect the 80-character line limit whenever possible.\n# - Be concise.\n# - Use lists instead of long sentences.\n#   - Use numbered lists when describing a sequence, e.g., (1) do X, (2) then Y.\n#   - Use bullet points when the items are not ordered.\n# - Example code should be covered by doctests.\n#   - But it\'s unclear what to do if the code depends on a solver, see\n#     https://github.com/JuliaOpt/JuMP.jl/issues/1175.warning: Warning\nThis documentation is for the development branch of JuMP. JuMP is undergoing a major transition to MathOptInterface, and the documentation is in the process of being rewritten. The development version is alpha-quality with breaking changes still in progress. Please provide feedback and file issues if you use this branch.JuMP is a domain-specific modeling language for mathematical optimization embedded in Julia. It currently supports a number of open-source and commercial solvers (see below) for a variety of problem classes, including linear programming, mixed-integer programming, second-order conic programming, semidefinite programming, and nonlinear programming. JuMP\'s features include:User friendliness\nSyntax that mimics natural mathematical expressions.\nComplete documentation (WIP!)\nSpeed\nBenchmarking has shown that JuMP can create problems at similar speeds   to special-purpose modeling languages such as   AMPL.\nJuMP communicates with most solvers in memory, avoiding the need to   write intermediary files.\nSolver independence\nJuMP uses a generic solver-independent interface provided by the   MathOptInterface   package, making it easy to change between a number of open-source and   commercial optimization software packages (\"solvers\").\nCurrently supported solvers include   Artelys Knitro,   Bonmin,   Cbc,   Clp,   Couenne,   CPLEX,   ECOS,   FICO Xpress,   GLPK,   Gurobi,   Ipopt,   MOSEK,   NLopt, and   SCS.\nAccess to advanced algorithmic techniques\nIncluding efficient LP re-solves which previously required using   solver-specific and/or low-level C++ libraries.\nEase of embedding\nJuMP itself is written purely in Julia. Solvers are the only binary   dependencies.\nBeing embedded in a general-purpose programming language makes it easy   to solve optimization problems as part of a larger workflow (e.g.,   inside a simulation, behind a web server, or as a subproblem in a   decomposition algorithm).\nAs a trade-off, JuMP\'s syntax is constrained by the syntax available   in Julia.\nJuMP is MPL licensed, meaning that   it can be embedded in commercial software that complies with the terms   of the license.While neither Julia nor JuMP have reached version 1.0 yet, the releases are stable enough for everyday use and are being used in a number of research projects and neat applications by a growing community of users who are early adopters. JuMP remains under active development, and we welcome your feedback, suggestions, and bug reports."
},

{
    "location": "index.html#Contents-1",
    "page": "Introduction",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"installation.md\",\n    \"quickstart.md\",\n    \"concepts.md\",\n    \"variables.md\",\n    \"expressions.md\",\n    \"constraints.md\",\n    \"containers.md\",\n    \"names.md\",\n    \"solvers.md\",\n    \"nlp.md\",\n    \"style.md\",\n    \"extensions.md\",\n    \"updating.md\",\n    \"howdoi.md\"]\nDepth = 2"
},

{
    "location": "index.html#Citing-JuMP-1",
    "page": "Introduction",
    "title": "Citing JuMP",
    "category": "section",
    "text": "If you find JuMP useful in your work, we kindly request that you cite the following paper (pdf):@article{DunningHuchetteLubin2017,\nauthor = {Iain Dunning and Joey Huchette and Miles Lubin},\ntitle = {JuMP: A Modeling Language for Mathematical Optimization},\njournal = {SIAM Review},\nvolume = {59},\nnumber = {2},\npages = {295-320},\nyear = {2017},\ndoi = {10.1137/15M1020575},\n}For an earlier work where we presented a prototype implementation of JuMP, see here:@article{LubinDunningIJOC,\nauthor = {Miles Lubin and Iain Dunning},\ntitle = {Computing in Operations Research Using Julia},\njournal = {INFORMS Journal on Computing},\nvolume = {27},\nnumber = {2},\npages = {238-248},\nyear = {2015},\ndoi = {10.1287/ijoc.2014.0623},\n}A preprint of this paper is freely available.(Image: NumFOCUS logo)JuMP is a fiscally sponsored project of NumFOCUS, a nonprofit dedicated to supporting the open source scientific computing community."
},

{
    "location": "installation.html#",
    "page": "Installation Guide",
    "title": "Installation Guide",
    "category": "page",
    "text": ""
},

{
    "location": "installation.html#Installation-Guide-1",
    "page": "Installation Guide",
    "title": "Installation Guide",
    "category": "section",
    "text": "JuMP is a package for Julia. To use JuMP, first download and install Julia or open up a remote notebook at JuliaBox or similar services.JuMP is compatible with both Julia 1.0 and 0.6. The compatibility with 0.6 is intended to facilitate upgrading from previous JuMP and Julia releases. The following instructions assume Julia 1.0.From Julia, JuMP is installed by using the built-in package manager:import Pkg\nPkg.add(\"JuMP\")note: Note\nThe installation instructions above assume that JuMP 0.19 has already been released. Until that time, see the JuMP README for instructions on installing a development release that\'s compatible with Julia 1.0."
},

{
    "location": "installation.html#Getting-Solvers-1",
    "page": "Installation Guide",
    "title": "Getting Solvers",
    "category": "section",
    "text": "JuMP depends on solvers to solve optimization problems. Most solvers are not written in Julia, and some require commercial licenses to use, so installation is often more complex. We list below the currently available solvers.note: Note\nThis list is open for new contributions. See also Interacting with solvers and the MathOptInterface docs for more details on how JuMP interacts with solvers. Please get in touch with any questions about connecting new solvers with JuMP.Solver Julia Package License Supports\nCbc Cbc.jl EPL MILP\nClp Clp.jl EPL LP\nCPLEX CPLEX.jl Comm. LP, MILP, SOCP, MISOCP\nCSDP CSDP.jl EPL LP, SDP\nECOS ECOS.jl GPL LP, SOCP\nFICO Xpress Xpress.jl Comm. LP, MILP, SOCP, MISOCP\nGLPK GLPK GPL LP, MILP\nGurobi Gurobi.jl Comm. LP, MILP, SOCP, MISOCP\nIpopt Ipopt.jl EPL LP, QP, NLP\nMOSEK MathOptInterfaceMosek.jl Comm. LP, MILP, SOCP, MISOCP, SDP\nOSQP OSQP.jl Apache LP, QP\nSCS SCS.jl MIT LP, SOCP, SDPWhere:LP = Linear programming\nQP = Quadratic programming\nSOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)\nMILP = Mixed-integer linear programming\nNLP = Nonlinear programming\nMINLP = Mixed-integer nonlinear programming\nSDP = Semidefinite programming\nMISDP = Mixed-integer semidefinite programmingTo install Gurobi, for example, and use it with a JuMP model model, run:import Pkg\nPkg.add(\"Gurobi\")\nusing JuMP\nusing Gurobi\nmodel = Model(with_optimizer(Gurobi.Optimizer))Most packages follow the ModuleName.Optimizer naming convention, but exceptions may exist. See the corresponding Julia package README for more details on how to use the solver.TODO: Discuss setting solver options.The following solvers were compatible with JuMP up to release 0.18 but are not yet compatible with the latest version because they do not implement the new MathOptInterface API:Artelys Knitro\nBARON\nBonmin and Couenne via AmplNLWriter.jl\nCDD\nNLopt\nPavito\nPajarito\nSCIPSolver-specific notes follow below."
},

{
    "location": "installation.html#Artelys-Knitro-1",
    "page": "Installation Guide",
    "title": "Artelys Knitro",
    "category": "section",
    "text": "Requires a license. The KNITRO.jl interface currently supports only nonlinear problems."
},

{
    "location": "installation.html#BARON-1",
    "page": "Installation Guide",
    "title": "BARON",
    "category": "section",
    "text": "Requires a license. A trial version is available for small problem instances."
},

{
    "location": "installation.html#COIN-OR-Cbc-1",
    "page": "Installation Guide",
    "title": "COIN-OR Cbc",
    "category": "section",
    "text": "Cbc supports \"SOS\" constraints."
},

{
    "location": "installation.html#CPLEX-1",
    "page": "Installation Guide",
    "title": "CPLEX",
    "category": "section",
    "text": "Requires a working installation of CPLEX with a license (free for faculty members and graduate teaching assistants). The interface requires using CPLEX as a shared library, which is unsupported by the CPLEX developers. Special installation steps are required on Mac OS. CPLEX supports \"SOS\" constraints."
},

{
    "location": "installation.html#ECOS-1",
    "page": "Installation Guide",
    "title": "ECOS",
    "category": "section",
    "text": "ECOS can be used by JuMP to solve LPs and SOCPs. ECOS does not support general quadratic objectives or constraints, only second-order conic constraints specified by using the SecondOrderCone set."
},

{
    "location": "installation.html#Gurobi-1",
    "page": "Installation Guide",
    "title": "Gurobi",
    "category": "section",
    "text": "￼Requires a working installation of Gurobi with an activated license (free for academic use). Gurobi supports \"SOS\" constraints."
},

{
    "location": "installation.html#FICO-Xpress-1",
    "page": "Installation Guide",
    "title": "FICO Xpress",
    "category": "section",
    "text": "Requires a working installation of Xpress with an active license (it is possible to get a license for academic use, see FICO Academic Partner Program). Supports SOCP and \"SOS\" constraints."
},

{
    "location": "installation.html#MOSEK-1",
    "page": "Installation Guide",
    "title": "MOSEK",
    "category": "section",
    "text": "￼Requires a license (free for academic use). The Mosek interface is maintained by the Mosek team. (Thanks!)"
},

{
    "location": "installation.html#SCS-1",
    "page": "Installation Guide",
    "title": "SCS",
    "category": "section",
    "text": "SCS can be used by JuMP to solve LPs and SOCPs, and SDPs. SCS is a first order solver and has low accuracy (10^4) by default; see the SCS.jl documentation for more information."
},

{
    "location": "quickstart.html#",
    "page": "Quick Start Guide",
    "title": "Quick Start Guide",
    "category": "page",
    "text": ""
},

{
    "location": "quickstart.html#Quick-Start-Guide-1",
    "page": "Quick Start Guide",
    "title": "Quick Start Guide",
    "category": "section",
    "text": "This quick start guide will introduce the main concepts of JuMP. If you are familiar with another modeling language embedded in a high-level language such as PuLP (Python) or a solver-specific interface you will find most of this familiar. If you are coming from an AMPL or similar background, you may find some of the concepts novel but the general appearance will still be familiar.The example in this guide is deliberately kept simple. There are more complex examples in the JuMP/examples/ folder.Once JuMP is installed, to use JuMP in your programs, you just need to say:julia> using JuMPYou also need to include a Julia package which provides an appropriate solver. One such solver is GLPK.Optimizer, which is provided by the GLPK.jl package.julia> using GLPKSee Installation Guide for a list of other solvers you can use.Models are created with the Model() function. The with_optimizer syntax is used to specify the optimizer to be used:julia> model = Model(with_optimizer(GLPK.Optimizer))\nA JuMP ModelDocTestSetup = quote\n    # Using a caching optimizer removes the need to # load a solver such as GLPK\n    # for building the documentation.\n    const MOI = JuMP.MathOptInterface\n    model = Model(with_optimizer(MOI.Utilities.MockOptimizer,\n                                 JuMP.JuMPMOIModel{Float64}(),\n                                 eval_objective_value = false,\n                                 eval_variable_constraint_dual = false))\nendnote: Note\nYour model doesn\'t have to be called model - it\'s just a name.There are a few options for defining a variable, depending on whether you want to have lower bounds, upper bounds, both bounds, or even no bounds. The following commands will create two variables, x and y, with both lower and upper bounds. Note the first argument is our model variable model. These variables are associated with this model and cannot be used in another model.julia> @variable(model, 0 <= x <= 2)\nx\n\njulia> @variable(model, 0 <= y <= 30)\nySee the Variables section for more information on creating variables.DocTestSetup = nothingNext we\'ll set our objective. Note again the model, so we know which model\'s objective we are setting! The objective sense, Max or Min, should be provided as the second argument. Note also that we don\'t have a multiplication * symbol between 5 and our variable x - Julia is smart enough to not need it! Feel free to stick with * if it makes you feel more comfortable, as we have done with 3 * y. (We have been intentionally inconsistent here to demonstrate different syntax; however, it is good practice to pick one way or the other consistently in your code.)julia> @objective(model, Max, 5x + 3 * y)Adding constraints is a lot like setting the objective. Here we create a less-than-or-equal-to constraint using <=, but we can also create equality constraints using == and greater-than-or-equal-to constraints with >=:julia> @constraint(model, con, 1x + 5y <= 3)\ncon : x + 5 y <= 3.0Note that in a similar manner to the @variable macro, we have named the constraint con. This will bind the constraint to the Julia variable con for later analysis.Models are solved with the JuMP.optimize! function:julia> JuMP.optimize!(model)DocTestSetup = quote\n    # Now we load in the solution. Using a caching optimizer removes the need to\n    # load a solver such as GLPK for building the documentation.\n    mock = JuMP.caching_optimizer(model).optimizer\n    MOI.set(mock, MOI.TerminationStatus(), MOI.Success)\n    MOI.set(mock, MOI.PrimalStatus(), MOI.FeasiblePoint)\n    MOI.set(mock, MOI.DualStatus(), MOI.FeasiblePoint)\n    MOI.set(mock, MOI.ResultCount(), 1)\n    MOI.set(mock, MOI.ObjectiveValue(), 10.6)\n    MOI.set(mock, MOI.VariablePrimal(), JuMP.optimizer_index(x), 2.0)\n    MOI.set(mock, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.2)\n    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(con), -0.6)\n    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(JuMP.UpperBoundRef(x)), -4.4)\n    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(JuMP.LowerBoundRef(y)), 0.0)\nendAfter the call to JuMP.optimize! has finished, we need to understand why the optimizer stopped. This can be for a number of reasons. First, the solver might have found the optimal solution, or proved that the problem is infeasible. However, it might also have run into numerical difficulties, or terminated due to a setting such as a time limit. We can ask the solver why it stopped using the JuMP.termination_status function:julia> JuMP.termination_status(model)\nSuccess::TerminationStatusCode = 0In this case, GLPK returned Success. This does not mean that it has found the optimal solution. Instead, it indicates that GLPK has finished running and did not encounter any errors or user-provided termination limits.DocTestSetup = nothingTo understand the reason for termination in more detail, we need to query JuMP.primalstatus:julia> JuMP.primal_status(model)\nFeasiblePoint::ResultStatusCode = 1This indicates that GLPK has found a FeasiblePoint to the primal problem. Coupled with the Success from JuMP.termination_status, we can infer that GLPK has indeed found the optimal solution. We can also query JuMP.dual_status:julia> JuMP.dual_status(model)\nFeasiblePoint::ResultStatusCode = 1Like the primal_status, GLPK indicates that it has found a FeasiblePoint to the dual problem.Finally, we can query the result of the optimization. First, we can query the objective value:julia> JuMP.objective_value(model)\n10.6We can also query the primal result values of the x and y variables:julia> JuMP.result_value(x)\n2.0\n\njulia> JuMP.result_value(y)\n0.2We can also query the value of the dual variable associated with the constraint con (which we bound to a Julia variable when defining the constraint):julia> JuMP.result_dual(con)\n-0.6To query the dual variables associated with the variable bounds, things are a little trickier as we first need to obtain a reference to the constraint:julia> x_upper = JuMP.UpperBoundRef(x)\nx <= 2.0\n\njulia> JuMP.result_dual(x_upper)\n-4.4A similar process can be followed to obtain the dual of the lower bound constraint on y:julia> y_lower = JuMP.LowerBoundRef(y)\ny >= 0.0\n\njulia> JuMP.result_dual(y_lower)\n0.0"
},

{
    "location": "concepts.html#",
    "page": "Concepts",
    "title": "Concepts",
    "category": "page",
    "text": ""
},

{
    "location": "concepts.html#Concepts-and-Definitions-1",
    "page": "Concepts",
    "title": "Concepts and Definitions",
    "category": "section",
    "text": "TODO: Use this section to define mathematical concepts used across JuMP. (?)MOI: MathOptInterface"
},

{
    "location": "variables.html#",
    "page": "Variables",
    "title": "Variables",
    "category": "page",
    "text": "CurrentModule = JuMP\nDocTestSetup = quote\n    using JuMP\nend"
},

{
    "location": "variables.html#Variables-1",
    "page": "Variables",
    "title": "Variables",
    "category": "section",
    "text": ""
},

{
    "location": "variables.html#What-is-a-JuMP-variable?-1",
    "page": "Variables",
    "title": "What is a JuMP variable?",
    "category": "section",
    "text": "The term variable in mathematical optimization has many meanings. Here, we distinguish between the following three types of variables:optimization variables, which are the mathematical x in the problem maxf_0(x)  f_i(x) in S_i.\nJulia variables, which are bindings between a name and a value, for example x = 1. (See here for the Julia docs.)\nJuMP variables, which are instances of the JuMP.VariableRef struct defined by JuMP that contains a reference to an optimization variable in a model. (Extra for experts: the VariableRef struct is a thin wrapper around a MOI.VariableIndex, and also contains a reference to the JuMP model.)To illustrate these three types of variables, consider the following JuMP code (the full syntax is explained below):julia> model = Model()\nA JuMP Model\n\njulia> @variable(model, x[1:2])\n2-element Array{VariableRef,1}:\n x[1]\n x[2]This code does three things:it adds two optimization variables to model\nit creates two JuMP variables that act as references to those optimization variables\nit binds those JuMP variables as a vector with two elements to the Julia variable x.To reduce confusion, we will attempt, where possible, to always refer to variables with their corresponding prefix.warn: Warn\nCreating two JuMP variables with the same name results in an error at runtime.JuMP variables can have attributes, such as names or an initial primal start value. We illustrate the name attribute in the following example:julia> @variable(model, y, base_name=\"decision variable\")\ndecision variableThis code does four things:it adds one optimization variable to model\nit creates one JuMP variable that acts as a reference to that optimization variable\nit binds the JuMP variable to the Julia variable y\nit tells JuMP that the name attribute of this JuMP variable is \"decisionvariable\". JuMP uses the value of base_name when it has to print the variable as a string.For example, when we print y at the REPL we get:julia> y\ndecision variableBecause y is a Julia variable, we can bind it to a different value. For example, if we write:julia> y = 1\n1y is no longer a binding to a JuMP variable. This does not mean that the JuMP variable has been destroyed. It still exists and is still a reference to the same optimization variable. The binding can be reset by querying the model for the symbol as it was written in the @variable macro. For example:julia> model[:y]\ndecision variableThis act of looking up the JuMP variable by using the symbol is most useful when composing JuMP models across multiple functions, as illustrated by the following example:function add_component_to_model(model::JuMP.Model)\n    x = model[:x]\n    # ... code that uses x\nend\nfunction build_model()\n    model = Model()\n    @variable(model, x)\n    add_component_to_model(model)\nend\n# TODO(@odow): add a section on looking up by stringNow that we understand the difference between optimization, JuMP, and Julia variables, we can introduce more of the functionality of the @variable macro."
},

{
    "location": "variables.html#Variable-bounds-1",
    "page": "Variables",
    "title": "Variable bounds",
    "category": "section",
    "text": "We have already seen the basic usage of the @variable macro. The next extension is to add lower- and upper-bounds to each optimization variable. This can be done as follows:julia> @variable(model, x_free)\nx_free\n\njulia> @variable(model, x_lower >= 0)\nx_lower\n\njulia> @variable(model, x_upper <= 1)\nx_upper\n\njulia> @variable(model, 2 <= x_interval <= 3)\nx_interval\n\njulia> @variable(model, x_fixed == 4)\nx_fixedIn the above examples, x_free represents an unbounded optimization variable, x_lower represents an optimization variable with a lower bound and so forth.note: Note\nWhen creating a variable with only a lower-bound or an upper-bound, and the value of the bound is not a numeric literal, the name must appear on the left-hand side. Putting the name on the right-hand side will result in an error. For example:@variable(model, 1 <= x)  # works\na = 1\n@variable(model, a <= x)  # errorsWe can query whether an optimization variable has a lower- or upper-bound via the JuMP.has_lower_bound and JuMP.has_upper_bound functions. For example:julia> JuMP.has_lower_bound(x_free)\nfalse\n\njulia> JuMP.has_upper_bound(x_upper)\ntrueIf a variable has a lower or upper bound, we can query the value of it via the JuMP.lower_bound and JuMP.upper_bound functions. For example:julia> JuMP.lower_bound(x_interval)\n2.0\n\njulia> JuMP.upper_bound(x_interval)\n3.0Querying the value of a bound that does not exist will result in an error.Instead of using the <= and >= syntax, we can also use the lower_bound and upper_bound keyword arguments. For example:julia> @variable(model, x, lower_bound=1, upper_bound=2)\nx\n\njulia> JuMP.lower_bound(x)\n1.0Another option is to use the JuMP.set_lower_bound and JuMP.set_upper_bound functions. These can also be used to modify an existing variable bound. For example:julia> @variable(model, x >= 1)\nx\n\njulia> JuMP.lower_bound(x)\n1.0\n\njulia> JuMP.set_lower_bound(x, 2)\n\njulia> JuMP.lower_bound(x)\n2.0Finally, we can delete variable bounds using JuMP.delete_lower_bound and JuMP.delete_upper_bound:julia> @variable(model, 1 <= x <= 2)\nx\n\njulia> JuMP.lower_bound(x)\n1.0\n\njulia> JuMP.delete_lower_bound(x)\n\njulia> JuMP.has_lower_bound(x)\nfalse\n\njulia> JuMP.upper_bound(x)\n2.0\n\njulia> JuMP.delete_upper_bound(x)\n\njulia> JuMP.has_upper_bound(x)\nfalse"
},

{
    "location": "variables.html#Variable-containers-1",
    "page": "Variables",
    "title": "Variable containers",
    "category": "section",
    "text": "In the examples above, we have mostly created scalar variables. By scalar, we mean that the Julia variable is bound to exactly one JuMP variable. However,  it is often useful to create collections of JuMP variables inside more complicated datastructures.JuMP provides a mechanism for creating three types of these datastructures, which we refer to as containers. The three types are Arrays, JuMPArrays, and Dictionaries. We explain each of these in the following."
},

{
    "location": "variables.html#Arrays-1",
    "page": "Variables",
    "title": "Arrays",
    "category": "section",
    "text": "We have already seen the creation of an array of JuMP variables with the x[1:2] syntax. This can naturally be extended to create multi-dimensional arrays of JuMP variables. For example:julia> @variable(model, x[1:2, 1:2])\n2×2 Array{VariableRef,2}:\n x[1,1]  x[1,2]\n x[2,1]  x[2,2]Arrays of JuMP variables can be indexed and sliced as follows:julia> x[1, 2]\nx[1,2]\n\njulia> x[2, :]\n2-element Array{VariableRef,1}:\n x[2,1]\n x[2,2]We can also name each index, and variable bounds can depend upon the indices:julia> @variable(model, x[i=1:2, j=1:2] >= 2i + j)\n2×2 Array{VariableRef,2}:\n x[1,1]  x[1,2]\n x[2,1]  x[2,2]\n\njulia> JuMP.lower_bound.(x)\n2×2 Array{Float64,2}:\n 3.0  4.0\n 5.0  6.0JuMP will form an Array of JuMP variables when it can determine at compile time that the indices are one-based integer ranges. Therefore x[1:b] will work, but x[a:b] will throw an error."
},

{
    "location": "variables.html#JuMPArrays-1",
    "page": "Variables",
    "title": "JuMPArrays",
    "category": "section",
    "text": "We often want to create arrays where the indices are not one-based integer ranges. For example, we may want to create a variable indexed by the name of a product or a location. The syntax is the same as that above, except with an arbitrary vector as an index as opposed to a one-based range. The biggest difference is that instead of returning an Array of JuMP variables, JuMP will return a JuMPArray. For example:julia> @variable(model, x[1:2, [:A,:B]])\n2-dimensional JuMPArray{VariableRef,2,...} with index sets:\n    Dimension 1, 1:2\n    Dimension 2, Symbol[:A, :B]\nAnd data, a 2×2 Array{VariableRef,2}:\n x[1,A]  x[1,B]\n x[2,A]  x[2,B]JuMPArray\'s can be indexed and sliced as follows:julia> x[1, :A]\nx[1,A]\n\njulia> x[2, :]\n1-dimensional JuMPArray{VariableRef,1,...} with index sets:\n    Dimension 1, Symbol[:A, :B]\nAnd data, a 2-element Array{VariableRef,1}:\n x[2,A]\n x[2,B]Similarly to the Array case, the indices in a JuMPArray can be named, and the bounds can depend upon these names. For example:julia> @variable(model, x[i=2:3, j=1:2:3] >= 0.5i + j)\n2-dimensional JuMPArray{VariableRef,2,...} with index sets:\n    Dimension 1, 2:3\n    Dimension 2, 1:2:3\nAnd data, a 2×2 Array{VariableRef,2}:\n x[2,1]  x[2,3]\n x[3,1]  x[3,3]\n\njulia> JuMP.lower_bound.(x)\n2-dimensional JuMPArray{Float64,2,...} with index sets:\n    Dimension 1, 2:3\n    Dimension 2, 1:2:3\nAnd data, a 2×2 Array{Float64,2}:\n 2.0  4.0\n 2.5  4.5"
},

{
    "location": "variables.html#Dictionaries-1",
    "page": "Variables",
    "title": "Dictionaries",
    "category": "section",
    "text": "The third datatype that JuMP supports the efficient creation of are dictionaries. These dictionaries are created when the indices do not form a rectangular set. One example is when indices have a dependence upon previous indices (called triangular indexing). JuMP supports this as follows:julia> @variable(model, x[i=1:2, j=i:2])\nDict{Any,VariableRef} with 3 entries:\n  (1, 2) => x[1,2]\n  (2, 2) => x[2,2]\n  (1, 1) => x[1,1]x is a standard Julia dictionary. Therefore, slicing cannot be performed.We can also conditionally create variables via a JuMP-specific syntax. This sytax appends a comparison check that depends upon the named indices and is separated from the indices by a semi-colon (;). For example:julia> @variable(model, x[i=1:4; mod(i, 2)==0])\nDict{Any,VariableRef} with 2 entries:\n  4 => x[4]\n  2 => x[2]"
},

{
    "location": "variables.html#Forcing-the-container-type-1",
    "page": "Variables",
    "title": "Forcing the container type",
    "category": "section",
    "text": "When creating a container of JuMP variables, JuMP will attempt to choose the tightest container type that can store the JuMP variables. Thus, it will prefer to create an Array before a JuMPArray, and a JuMPArray before a dictionary. However, because this happens at compile time, it does not always make the best choice. To illustrate this, consider the following example:julia> A = 1:2\n1:2\n\njulia> @variable(model, x[A])\n1-dimensional JuMPArray{VariableRef,1,...} with index sets:\n    Dimension 1, 1:2\nAnd data, a 2-element Array{VariableRef,1}:\n x[1]\n x[2]Since the value (and type) of A is unknown at compile time, JuMP is unable to infer that A is a one-based integer range. Therefore, JuMP creates a JuMPArray, even though it could store these two variables in a standard one-dimensional Array.We can share our knowledge that it is possible to store these JuMP variables as an array by setting the container keyword:julia> @variable(model, y[A], container=Array)\n2-element Array{VariableRef,1}:\n y[1]\n y[2]JuMP now creates a vector of JuMP variables, instead of a JuMPArray. Note that choosing an invalid container type will throw an error."
},

{
    "location": "variables.html#Integrality-shortcuts-1",
    "page": "Variables",
    "title": "Integrality shortcuts",
    "category": "section",
    "text": "Adding integrality constraints to a model such as @constraint(model, x in MOI.ZeroOne()) and @constraint(model, x in MOI.Integer()) is a common operation. Therefore, JuMP supports two shortcuts for adding such constraints."
},

{
    "location": "variables.html#Binary-(ZeroOne)-constraints-1",
    "page": "Variables",
    "title": "Binary (ZeroOne) constraints",
    "category": "section",
    "text": "Binary optimization variables are constrained to the set x in 0 1. (The MOI.ZeroOne set in MathOptInterface.) Binary optimization variables can be created in JuMP by passing Bin as an optional positional argument:julia> @variable(model, x, Bin)\nxWe can check if an optimization variable is binary by calling JuMP.is_binary on the JuMP variable:julia> JuMP.is_binary(x)\ntrueBinary optimization variables can also be created by setting the binary keyword to true.julia> @variable(model, x, binary=true)\nx"
},

{
    "location": "variables.html#Integer-constraints-1",
    "page": "Variables",
    "title": "Integer constraints",
    "category": "section",
    "text": "Integer optimization variables are constrained to the set x in mathbbZ. (The MOI.Integer set in MathOptInterface.) Integer optimization variables can be created in JuMP by passing Int as an optional positional argument:julia> @variable(model, x, Int)\nxInteger optimization variables can also be created by setting the integer keyword to true.julia> @variable(model, x, integer=true)\nxWe can check if an optimization variable is integer by calling JuMP.is_integer on the JuMP variable:julia> JuMP.is_integer(x)\ntrue"
},

{
    "location": "variables.html#Semidefinite-variables-1",
    "page": "Variables",
    "title": "Semidefinite variables",
    "category": "section",
    "text": "JuMP also supports modeling with semidefinite variables. A square symmetric matrix X is positive semidefinite if all eigenvalues are nonnegative. We can declare a matrix of JuMP variables to be positive semidefinite as follows:julia> @variable(model, x[1:2, 1:2], PSD)\n2×2 LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}:\n x[1,1]  x[1,2]\n x[1,2]  x[2,2]Note that x must be a square 2-dimensional Array of JuMP variables; it cannot be a JuMP array or a dictionary. (See Variable containers, above, for more on this.)You can also impose a slightly weaker constraint that the square matrix is only symmetric (instead of positive semidefinite) as follows:julia> @variable(model, x[1:2, 1:2], Symmetric)\n2×2 LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}:\n x[1,1]  x[1,2]\n x[1,2]  x[2,2]"
},

{
    "location": "variables.html#Anonymous-JuMP-variables-1",
    "page": "Variables",
    "title": "Anonymous JuMP variables",
    "category": "section",
    "text": "In all of the above examples, we have created named JuMP variables. However, it is also possible to create so called anonymous JuMP variables. To create an anonymous JuMP variable, we drop the name of the variable from the macro call. This means dropping the second positional argument if the JuMP variable is a scalar, or dropping the name before the square bracket ([) if a container is being created. For example:julia> x = @variable(model)\nnonameThis shows how (model, x) is really short for:julia> x = model[:x] = @variable(model, base_name=\"x\")\nxAn Array of anonymous JuMP variables can be created as follows:julia> y = @variable(model, [i=1:2])\n2-element Array{VariableRef,1}:\n noname\n nonameIf necessary, you can store x in model as follows:julia> model[:x] = xThe <= and >= short-hand cannot be used to set bounds on anonymous JuMP variables. Instead, you should use the lower_bound and upper_bound keywords.Passing the Bin and Int variable types are also invalid. Instead, you should use the binary and integer keywords.Thus, the anonymous variant of @variable(model, x[i=1:2] >= i, Int) is:julia> x = @variable(model, [i=1:2], base_name=\"x\", lower_bound=i, integer=true)\n2-element Array{VariableRef,1}:\n x[1]\n x[2]"
},

{
    "location": "variables.html#User-defined-containers-1",
    "page": "Variables",
    "title": "User-defined containers",
    "category": "section",
    "text": "In the section Variable containers, we explained how JuMP supports the efficient creation of collections of JuMP variables in three types of containers. However, users are also free to create collections of JuMP variables in their own datastructures. For example, the following code creates a dictionary with symmetric matrices as the values:julia> variables = Dict{Symbol, Array{VariableRef,2}}()\nDict{Symbol,Array{VariableRef,2}} with 0 entries\n\njulia> for key in [:A, :B]\n           global variables[key] = @variable(model, [1:2, 1:2])\n       end\n\njulia> variables\nDict{Symbol,Array{VariableRef,2}} with 2 entries:\n  :A => VariableRef[noname noname; noname noname]\n  :B => VariableRef[noname noname; noname noname]"
},

{
    "location": "variables.html#Deleting-variables-1",
    "page": "Variables",
    "title": "Deleting variables",
    "category": "section",
    "text": "JuMP supports the deletion of optimization variables.  To delete variables, we can use the JuMP.delete method. We can also check whether x is a valid JuMP variable in model using the JuMP.is_valid method:julia> @variable(model, x)\nx\n\njulia> JuMP.is_valid(model, x)\ntrue\n\njulia> JuMP.delete(model, x)\n\njulia> JuMP.is_valid(model, x)\nfalse"
},

{
    "location": "variables.html#JuMP.@variable",
    "page": "Variables",
    "title": "JuMP.@variable",
    "category": "macro",
    "text": "@variable(model, kwargs...)\n\nAdd an anonymous (see Names) variable to the model model described by the keyword arguments kwargs and returns the variable.\n\n@variable(model, expr, args..., kwargs...)\n\nAdd a variable to the model model described by the expression expr, the positional arguments args and the keyword arguments kwargs. The expression expr can either be (note that in the following the symbol <= can be used instead of ≤ and the symbol >=can be used instead of ≥)\n\nof the form varexpr creating variables described by varexpr;\nof the form varexpr ≤ ub (resp. varexpr ≥ lb) creating variables described by varexpr with upper bounds given by ub (resp. lower bounds given by lb);\nof the form varexpr == value creating variables described by varexpr with fixed values given by value; or\nof the form lb ≤ varexpr ≤ ub or ub ≥ varexpr ≥ lb creating variables described by varexpr with lower bounds given by lb and upper bounds given by ub.\n\nThe expression varexpr can either be\n\nof the form varname creating a scalar real variable of name varname;\nof the form varname[...] or [...] creating a container of variables (see Containers in macro.\n\nThe recognized positional arguments in args are the following:\n\nBin: Sets the variable to be binary, i.e. either 0 or 1.\nInt: Sets the variable to be integer, i.e. one of ..., -2, -1, 0, 1, 2, ...\nSymmetric: Only available when creating a square matrix of variables, i.e. when varexpr is of the form varname[1:n,1:n] or varname[i=1:n,j=1:n]. It creates a symmetric matrix of variable, that is, it only creates a new variable for varname[i,j] with i ≤ j and sets varname[j,i] to the same variable as varname[i,j].\nPSD: The square matrix of variable is both Symmetric and constrained to be positive semidefinite.\n\nThe recognized keyword arguments in kwargs are the following:\n\nbase_name: Sets the name prefix used to generate variable names. It corresponds to the variable name for scalar variable, otherwise, the variable names are set to base_name[...] for each index ... of the axes axes.\nlower_bound: Sets the value of the variable lower bound.\nupper_bound: Sets the value of the variable upper bound.\nstart: Sets the variable starting value used as initial guess in optimization.\nbinary: Sets whether the variable is binary or not.\ninteger: Sets whether the variable is integer or not.\nvariable_type: See the \"Note for extending the variable macro\" section below.\ncontainer: Specify the container type, see Containers in macro.\n\nExamples\n\nThe following are equivalent ways of creating a variable x of name x with lower bound 0:\n\n# Specify everything in `expr`\n@variable(model, x >= 0)\n# Specify the lower bound using a keyword argument\n@variable(model, x, lower_bound=0)\n# Specify everything in `kwargs`\nx = @variable(model, base_name=\"x\", lower_bound=0)\n\nThe following are equivalent ways of creating a JuMPArray of index set [:a, :b] and with respective upper bounds 2 and 3 and names x[a] and `x[b].\n\nub = Dict(:a => 2, :b => 3)\n# Specify everything in `expr`\n@variable(model, x[i=keys(ub)] <= ub[i])\n# Specify the upper bound using a keyword argument\n@variable(model, x[i=keys(ub)], upper_bound=ub[i])\n\nNote for extending the variable macro\n\nThe single scalar variable or each scalar variable of the container are created using add_variable(model, build_variable(_error, info, extra_args...; extra_kwargs...)) where\n\nmodel is the model passed to the @variable macro;\n_error is an error function with a single String argument showing the @variable call in addition to the error message given as argument;\ninfo is the VariableInfo struct containing the information gathered in expr, the recognized keyword arguments (except base_name and variable_type) and the recognized positional arguments (except Symmetric and PSD);\nextra_args are the unrecognized positional arguments of args plus the value of the variable_type keyword argument if present. The variable_type keyword argument allows the user to pass a position argument to build_variable without the need to give a positional argument to @variable. In particular, this allows the user to give a positional argument to the build_variable call when using the anonymous single variable syntax @variable(model, kwargs...); and\nextra_kwargs are the unrecognized keyword argument of kwargs.\n\nExamples\n\nThe following creates a variable x of name x with lower_bound 0 as with the first example above but does it without using the @variable macro\n\ninfo = VariableInfo(true, 0, false, NaN, false, NaN, false, NaN, false, false)\nJuMP.add_variable(model, JuMP.build_variable(error, info), \"x\")\n\nThe following creates a JuMPArray of index set [:a, :b] and with respective upper bounds 2 and 3 and names x[a] and x[b] as with the second example above but does it without using the @variable macro\n\n# Without the `@variable` macro\ndata = Vector{JuMP.variable_type(model)}(undef, length(keys(ub)))\nx = JuMPArray(data, keys(ub))\nfor i in keys(ub)\n    info = VariableInfo(false, NaN, true, ub[i], false, NaN, false, NaN, false, false)\n    x[i] = JuMP.add_variable(model, JuMP.build_variable(error, info), \"x[$i]\")\nend\n\nThe following are equivalent ways of creating a Matrix of size N x N with variables custom variables created with a JuMP extension using the Poly(X) positional argument to specify its variables:\n\n# Using the `@variable` macro\n@variable(model, x[1:N,1:N], Symmetric, Poly(X))\n# Without the `@variable` macro\nx = Matrix{JuMP.variable_type(model, Poly(X))}(N, N)\ninfo = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)\nfor i in 1:N, j in i:N\n    x[i,j] = x[j,i] = JuMP.add_variable(model, build_variable(error, info, Poly(X)), \"x[$i,$j]\")\nend\n\n\n\n\n\n"
},

{
    "location": "variables.html#JuMP.owner_model",
    "page": "Variables",
    "title": "JuMP.owner_model",
    "category": "function",
    "text": "owner_model(s::AbstractJuMPScalar)\n\nReturn the model owning the scalar s.\n\n\n\n\n\n"
},

{
    "location": "variables.html#Reference-1",
    "page": "Variables",
    "title": "Reference",
    "category": "section",
    "text": "@variable\nowner_model"
},

{
    "location": "expressions.html#",
    "page": "Expressions",
    "title": "Expressions",
    "category": "page",
    "text": "DocTestSetup = quote\n    using JuMP\nend"
},

{
    "location": "expressions.html#Expressions-1",
    "page": "Expressions",
    "title": "Expressions",
    "category": "section",
    "text": "DRAFT: JuMP has multiple types of expressions: affine, quadratic, and nonlinear. Just talk about affine and quadratic here (see other section for nonlinear). Describe the basic data structures and ways to construct expressions (these are: constructors, operators, add_to_expression!, and macros).Example code with doc tests:m = Model()\n@variable(m, x)\n@variable(m, y)\nAffExpr(-1.0, x => 2.0, y => 1.0)\n\n# output\n\n2 x + y - 1"
},

{
    "location": "expressions.html#Objective-functions-1",
    "page": "Expressions",
    "title": "Objective functions",
    "category": "section",
    "text": "TODO: Describe how JuMP expressions relate to MOI functions. How to set, query, and modify an objective function."
},

{
    "location": "constraints.html#",
    "page": "Constraints",
    "title": "Constraints",
    "category": "page",
    "text": ""
},

{
    "location": "constraints.html#JuMP.@constraint",
    "page": "Constraints",
    "title": "JuMP.@constraint",
    "category": "macro",
    "text": "@constraint(m::Model, expr)\n\nAdd a constraint described by the expression expr.\n\n@constraint(m::Model, ref[i=..., j=..., ...], expr)\n\nAdd a group of constraints described by the expression expr parametrized by i, j, ...\n\nThe expression expr can either be\n\nof the form func in set constraining the function func to belong to the set set, e.g. @constraint(m, [1, x-1, y-2] in MOI.SecondOrderCone(3)) constrains the norm of [x-1, y-2] be less than 1;\nof the form a sign b, where sign is one of ==, ≥, >=, ≤ and <= building the single constraint enforcing the comparison to hold for the expression a and b, e.g. @constraint(m, x^2 + y^2 == 1) constrains x and y to lie on the unit circle;\nof the form a ≤ b ≤ c or a ≥ b ≥ c (where ≤ and <= (resp. ≥ and >=) can be used interchangeably) constraining the paired the expression b to lie between a and c;\nof the forms @constraint(m, a .sign b) or @constraint(m, a .sign b .sign c) which broadcast the constraint creation to each element of the vectors.\n\nNote for extending the constraint macro\n\nEach constraint will be created using add_constraint(m, build_constraint(_error, func, set)) where\n\n_error is an error function showing the constraint call in addition to the error message given as argument,\nfunc is the expression that is constrained\nand set is the set in which it is constrained to belong.\n\nFor expr of the first type (i.e. @constraint(m, func in set)), func and set are passed unchanged to build_constraint but for the other types, they are determined from the expressions and signs. For instance, @constraint(m, x^2 + y^2 == 1) is transformed into add_constraint(m, build_constraint(_error, x^2 + y^2, MOI.EqualTo(1.0))).\n\nTo extend JuMP to accept new constraints of this form, it is necessary to add the corresponding methods to build_constraint. Note that this will likely mean that either func or set will be some custom type, rather than e.g. a Symbol, since we will likely want to dispatch on the type of the function or set appearing in the constraint.\n\n\n\n\n\n"
},

{
    "location": "constraints.html#JuMP.@SDconstraint",
    "page": "Constraints",
    "title": "JuMP.@SDconstraint",
    "category": "macro",
    "text": "@SDconstraint(m::Model, expr)\n\nAdd a semidefinite constraint described by the expression expr.\n\n@SDconstraint(m::Model, ref[i=..., j=..., ...], expr)\n\nAdd a group of semidefinite constraints described by the expression expr parametrized by i, j, ...\n\nThe expression expr needs to be of the form a sign b where sign is ⪰, ≥, >=, ⪯, ≤ or <= and a and b are square matrices. It constrains that a - b (or b - a if the sign is ⪯, ≤ or <=) is positive semidefinite.\n\n\n\n\n\n"
},

{
    "location": "constraints.html#Constraints-1",
    "page": "Constraints",
    "title": "Constraints",
    "category": "section",
    "text": "DRAFT: Describe how constraints are represented (link to MOI docs). Constraints are very similar to variables in (1) how names work (2) how attributes work, and (3) the macro syntax for constructing them. They\'re a bit different because they\'re parameterized by function-set type. Describe constraints vs. ConstraintRefs. Describe JuMP.constraint_object. How to delete constraints. How to modify constraints by setting attributes and MOI.modifyconstraint!. Describe semidefinite constraints and symmetry handling. Refer to NLP docs for nonlinear constraints.@constraint\n@SDconstraint"
},

{
    "location": "containers.html#",
    "page": "Containers",
    "title": "Containers",
    "category": "page",
    "text": ""
},

{
    "location": "containers.html#JuMP.generatecontainer",
    "page": "Containers",
    "title": "JuMP.generatecontainer",
    "category": "function",
    "text": "generatecontainer(T, indexvars, indexsets, requestedtype)\n\nReturn a tuple, the first element of which is code that generates a container for objects of type T given the index variables, index sets, and requestedtype. requestedtype may be one of :Array, :JuMPArray, :Dict, or :Auto. Return error-producing code if requested type is incompatible. For the case of :Auto, the following rules are used to determine the appropriate container:\n\nIf all index sets are either explicit 1:B objects for any B or symbols which refer to objects of type Base.OneTo, then an Array is generated of the appropriate size. Types of symbols/expressions are not known at compile time, so we defer to type-safe functions to check the Base.OneTo condition.\nIf condition (1) does not hold, and the index sets are independent (the index variable for one set does not appear in the definition of another), then an JuMPArray is generated of the appropriate size.\nOtherwise, generate an empty Dict{Any,T}.\n\nThe second element of the return tuple is a Bool, true if the container type automatically checks for duplicate terms in the index sets and false otherwise.\n\nExamples\n\ngeneratecontainer(VariableRef, [:i,:j], [:(1:N), :(1:T)], :Auto)\n# Returns code equivalent to:\n# :(Array{VariableRef}(length(1:N), length(1:T))\n\ngeneratecontainer(VariableRef, [:i,:j], [:(1:N), :(2:T)], :Auto)\n# Returns code equivalent to:\n# :(JuMPArray(Array{VariableRef}(length(1:N), length(2:T)), $indexvars...))\n\ngeneratecontainer(VariableRef, [:i,:j], [:(1:N), :(S)], :Auto)\n# Returns code that generates an Array if S is of type Base.OneTo,\n# otherwise an JuMPArray.\n\ngeneratecontainer(VariableRef, [:i,:j], [:(1:N), :(1:j)], :Auto)\n# Returns code equivalent to:\n# :(Dict{Any,VariableRef}())\n\n\n\n\n\n"
},

{
    "location": "containers.html#Containers-1",
    "page": "Containers",
    "title": "Containers",
    "category": "section",
    "text": "Containers can be created using the generatecontainer functionJuMP.generatecontainer"
},

{
    "location": "containers.html#Containers-in-macro-1",
    "page": "Containers",
    "title": "Containers in macro",
    "category": "section",
    "text": "In the @variable (resp. @constraint) macro, containers of variables (resp. constraints) can be created the following syntaxname[index_set_1,index_set_2,...,index_set_n] creating an n-dimensional container of name name; or\n[index_set_1,index_set_2,...,index_set_n] creating an anonymous (see Names) n-dimensional container.Each expression index_set_i can either beof the form index_set specifying that the ith index set of the container is index_set; or\nof the form index_name=index_set specifying that the ith index set of the container is index_set and allowing values used in the macro expression and keyword arguments to be expressions depending on the index_name.The macro then creates the container using the JuMP.generatecontainer function with the following arguments:VariableRef for the @variable macro and ConstraintRef for the @constraint macro.\nThe index variables and arbitrary symbols for dimensions for which no variable index is specified.\nThe index sets specified.\nThe value of the keyword argument if given or :Auto."
},

{
    "location": "names.html#",
    "page": "Names",
    "title": "Names",
    "category": "page",
    "text": ""
},

{
    "location": "names.html#Names-1",
    "page": "Names",
    "title": "Names",
    "category": "section",
    "text": "TODO: This discussion is redundant with and inconsistent with the discussion of names for variables. Also mention that expressions can be anonymous or named.There a two different naming aspects that need to be distinguished when creating variables/contraints (resp. a container of variables/constraints):The name of the local variable created (if any) holding the reference (resp. the container of references) which corresponds to the name that can be used to retrieve it using m[:name].\nThe name of the variable/constraint (resp. each variable/constraint in the container) used for printing. This corresponds to the MOI.VariableName/MOI.ConstraintName attribute.When creating a variable using the syntax @variable(m; kwargs...), creating a constraint using the syntax @constraint(m, expr) or when creating a container with the syntax [...] in a macro, we say that the variable or constraint is anonymous. For anonymous variables/constraints, no local variable is created holding the reference or container of references and it is not stored in the model, i.e. it is not possible to retrieve it using m[:name].Otherwise, when it is not anonymous, the name used both for the local variable created and the key for retrieving the reference or container of references in the model are determined from the macro expression. For instance, when creating a container with the syntax name[...] or when creating a constraint with @constraint(m, name, expr), the name used is name.The name of the variable/constraint used for printing is based on the base name which is specified by the base_name keyword argument. When the base_name keyword argument is not specified, the name depends on whether the variable is anonymous:if the variable/constraint is anonymous, then the MOI.VariableName/MOI.ConstraintName attribute is not set and the name used for printing is noname,\notherwise, the base name is set to the name used for the local variable created.The name of the variable/constraint set to the MOI.VariableName/MOI.ConstraintName attribute and used for printing is then base_name for single variable/constraint and base_name[i1,i2,...,in] for the reference at indices i1, i2, ..., in in a container."
},

{
    "location": "solvers.html#",
    "page": "Solvers",
    "title": "Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers.html#Interacting-with-solvers-1",
    "page": "Solvers",
    "title": "Interacting with solvers",
    "category": "section",
    "text": "A JuMP model keeps a MathOptInterface (MOI) backend of type MOI.ModelLike internally that stores the optimization problem and acts as the optimization solver. We call it an MOI backend and not optimizer as it can also be a wrapper around an optimization file format such as MPS that writes the JuMP model in a file. JuMP can be viewed as a lightweight user-friendly layer on top of the MOI backend:JuMP does not maintain any copy of the model outside this MOI backend.\nJuMP variable (resp. constraint) references are simple structures containing both a reference to the JuMP model and the MOI index of the variable (resp. constraint).\nJuMP gives the constraints to the MOI backend in the form provided by the user without doing any automatic reformulation.\nvariables additions, constraints additions/modifications and objective modifications are directly applied to the MOI backend thus expecting the backend to support such modifications.While this allows JuMP to be a thin wrapper on top of the solver API, as mentioned in the last point above, this seems rather demanding on the solver. Indeed, while some solvers support incremental building of the model and modifications before and after solve, other solvers only support the model being copied at once before solve. Moreover it seems to require all solvers to implement all possible reformulations independently which seems both very ambitious and might generate a lot of duplicated code.These apparent limitations are in fact addressed at the MOI level in a manner that is completely transparent to JuMP. While the MOI API may seem very demanding, it allows MOI models to be a succession of lightweight MOI layers that fill the gap between JuMP requirements and the solver capabilities.JuMP models can be created in three different modes: Automatic, Manual and Direct."
},

{
    "location": "solvers.html#JuMP.with_optimizer",
    "page": "Solvers",
    "title": "JuMP.with_optimizer",
    "category": "function",
    "text": "with_optimizer(constructor, args...; kwargs...)\n\nReturn an OptimizerFactory that creates optimizers using the constructor constructor with positional arguments args and keyword arguments kwargs.\n\nExamples\n\nThe following returns an optimizer factory that creates IpoptOptimizers using the constructor call IpoptOptimizer(print_level=0):\n\nwith_optimizer(IpoptOptimizer, print_level=0)\n\n\n\n\n\n"
},

{
    "location": "solvers.html#JuMP.optimize!",
    "page": "Solvers",
    "title": "JuMP.optimize!",
    "category": "function",
    "text": "function optimize!(model::Model,\n                   optimizer_factory::Union{Nothing, OptimizerFactory}=nothing;\n                   ignore_optimize_hook=(model.optimize_hook === nothing))\n\nOptimize the model. If optimizer_factory is not nothing, it first set the optimizer to a new one created using the optimizer factory.\n\n\n\n\n\n"
},

{
    "location": "solvers.html#JuMP.Model-Tuple{}",
    "page": "Solvers",
    "title": "JuMP.Model",
    "category": "method",
    "text": "Model(; caching_mode::MOIU.CachingOptimizerMode=MOIU.Automatic,\n        bridge_constraints::Bool=true)\n\nReturn a new JuMP model without any optimizer; the model is stored the model in a cache. The mode of the CachingOptimizer storing this cache is caching_mode. The optimizer can be set later in the JuMP.optimize! call. If bridge_constraints is true, constraints that are not supported by the optimizer are automatically bridged to equivalent supported constraints when an appropriate is defined in the MathOptInterface.Bridges module or is defined in another module and is explicitely added.\n\n\n\n\n\n"
},

{
    "location": "solvers.html#JuMP.Model-Tuple{JuMP.OptimizerFactory}",
    "page": "Solvers",
    "title": "JuMP.Model",
    "category": "method",
    "text": "Model(optimizer_factory::OptimizerFactory;\n      caching_mode::MOIU.CachingOptimizerMode=MOIU.Automatic,\n      bridge_constraints::Bool=true)\n\nReturn a new JuMP model using the optimizer factory optimizer_factory to create the optimizer. The optimizer factory can be created by the with_optimizer function.\n\nExamples\n\nThe following creates a model using the optimizer IpoptOptimizer(print_level=0):\n\nmodel = JuMP.Model(with_optimizer(IpoptOptimizer, print_level=0))\n\n\n\n\n\n"
},

{
    "location": "solvers.html#Automatic-and-Manual-modes-1",
    "page": "Solvers",
    "title": "Automatic and Manual modes",
    "category": "section",
    "text": "In Automatic and Manual modes, two MOI layers are automatically applied to the optimizer:CachingOptimizer: maintains a cache of the model so that when the optimizer does not support an incremental change to the model, the optimizer\'s internal model can be discarded and restored from the cache just before optimization. The CachingOptimizer has two different modes: Automatic and Manual corresponding to the two JuMP modes with the same names.\nLazyBridgeOptimizer (this can be disabled using the bridge_constraints keyword argument to Model constructor): when a constraint added is not supported by the optimizer, it tries transform the constraint into an equivalent form, possibly adding new variables and constraints that are supported by the optimizer. The applied transformations are selected among known recipes which are called bridges. A few default bridges are defined in MOI but new ones can be defined and added to the LazyBridgeOptimizer used by JuMP.See the MOI documentation for more details on these two MOI layers.To attach an optimizer to a JuMP model, JuMP needs to create a new empty optimizer instance. New optimizer instances can be obtained using an OptimizerFactory that can be created using the with_optimizer function:with_optimizerThe factory can be provided either at model construction time or at JuMP.optimize! time:JuMP.optimize!New JuMP models are created using the Model constructor:Model()\nModel(::JuMP.OptimizerFactory)TODO: how to control the caching optimizer states"
},

{
    "location": "solvers.html#JuMP.direct_model",
    "page": "Solvers",
    "title": "JuMP.direct_model",
    "category": "function",
    "text": "direct_model(backend::MOI.ModelLike)\n\nReturn a new JuMP model using backend to store the model and solve it. As opposed to the Model constructor, no cache of the model is stored outside of backend and no bridges are automatically applied to backend. The absence of cache reduces the memory footprint but it is important to bear in mind the following implications of creating models using this direct mode:\n\nWhen backend does not support an operation such as adding variables/constraints after solver or modifying constraints, an error is thrown. With models created using the Model constructor, such situations can be dealt with by storing the modifications in a cache and loading them into the optimizer when JuMP.optimize! is called.\nNo constraint bridging is supported by default.\nThe optimizer used cannot be changed the model is constructed.\nThe model created cannot be copied.\n\n\n\n\n\n"
},

{
    "location": "solvers.html#Direct-mode-1",
    "page": "Solvers",
    "title": "Direct mode",
    "category": "section",
    "text": "JuMP models can be created in Direct mode using the JuMP.direct_model function.JuMP.direct_modelTODO: How to set parameters (solver specific and generic). Status codes. Accessing the result. How to accurately measure the solve time."
},

{
    "location": "nlp.html#",
    "page": "Nonlinear Modeling",
    "title": "Nonlinear Modeling",
    "category": "page",
    "text": ""
},

{
    "location": "nlp.html#Nonlinear-Modeling-1",
    "page": "Nonlinear Modeling",
    "title": "Nonlinear Modeling",
    "category": "section",
    "text": "CurrentModule = JuMP\nDocTestSetup = quote\n    using JuMP\nendJuMP has support for general smooth nonlinear (convex and nonconvex) optimization problems. JuMP is able to provide exact, sparse second-order derivatives to solvers. This information can improve solver accuracy and performance.Nonlinear objectives and constraints are specified by using the @NLobjective and @NLconstraint macros. The familiar sum() syntax is supported within these macros, as well as prod() which analogously represents the product of the terms within. Note that the @objective and @constraint macros (and corresponding functions) do not currently support nonlinear expressions. However, a model can contain a mix of linear, quadratic, and nonlinear contraints or objective functions. Starting points may be provided by using the start keyword argument to @variable.  For nonconvex problems, the returned solution is only guaranteed to be locally optimal. Convexity detection is not currently provided.TODO(issue #1460): Describe how starting points are computed if none are provided.For example, we can solve the classical Rosenbrock problem (with a twist) as follows:using Ipopt\nmodel = Model(with_optimizer(Ipopt.Optimizer))\n@variable(model, x, start = 0.0)\n@variable(model, y, start = 0.0)\n\n@NLobjective(model, Min, (1 - x) ^ 2 + 100 * (y - x ^ 2) ^ 2)\n\nJuMP.optimize!(model)\nprintln(\"x = \", JuMP.result_value(x), \" y = \", JuMP.result_value(y))\n\n# adding a (linear) constraint\n@constraint(model, x + y == 10)\nJuMP.optimize!(model)\nprintln(\"x = \", JuMP.result_value(x), \" y = \", JuMP.result_value(y))TODO: Add links to NLP examples after they are updated.The NLP solver tests contain additional examples."
},

{
    "location": "nlp.html#Syntax-notes-1",
    "page": "Nonlinear Modeling",
    "title": "Syntax notes",
    "category": "section",
    "text": "The syntax accepted in nonlinear expressions is more restricted than the syntax for linear and quadratic expressions. We note some important points below.All expressions must be simple scalar operations. You cannot use dot, matrix-vector products, vector slices, etc. Translate vector operations into explicit sum() operations or use the AffExpr plus auxiliary variable trick described below.\nThere is no operator overloading provided to build up nonlinear expressions. For example, if x is a JuMP variable, the code 3x will return an AffExpr object that can be used inside of future expressions and linear constraints. However, the code sin(x) is an error. All nonlinear expressions must be inside of macros.\nUser-defined Functions may be used within nonlinear expressions only after they are registered. For example, the follow code results in an error because JuMP.register() must be called first to register my_function.model = Model()\nmy_function(a, b) = exp(a) * b\n@variable(model, x)\n@variable(model, y)\n@NLobjective(model, Min, my_function(x, y))\n\n# output\n\nERROR: Unrecognized function \"my_function\" used in nonlinear expression.AffExpr and QuadExpr objects cannot currently be used inside nonlinear expressions. Instead, introduce auxiliary variables, e.g.:    my_expr = dot(c, x) + 3y # where x and y are variables\n    @variable(model, aux)\n    @constraint(model, aux == my_expr)\n    @NLobjective(model, Min, sin(aux))You can declare embeddable nonlinear expressions with @NLexpression. For example:    @NLexpression(model, my_expr[i = 1:n], sin(x[i]))\n    @NLconstraint(model, my_constr[i = 1:n], my_expr[i] <= 0.5)Anonymous syntax is supported in @NLexpression and @NLconstraint:    my_expr = @NLexpression(model, [i = 1:n], sin(x[i]))\n    my_constr = @NLconstraint(model, [i = 1:n], my_expr[i] <= 0.5)"
},

{
    "location": "nlp.html#JuMP.@NLparameter",
    "page": "Nonlinear Modeling",
    "title": "JuMP.@NLparameter",
    "category": "macro",
    "text": "@NLparameter(model, param == value)\n\nCreate and return a nonlinear parameter param attached to the model model with initial value set to value. Nonlinear parameters may be used only in nonlinear expressions.\n\nExample\n\nmodel = Model()\n@NLparameter(model, x == 10)\nJuMP.value(x)\n\n# output\n10.0\n\n@NLparameter(model, param_collection[...] == value_expr)\n\nCreate and return a collection of nonlinear parameters param_collection attached to the model model with initial value set to value_expr (may depend on index sets). Uses the same syntax for specifying index sets as @variable.\n\nExample\n\nmodel = Model()\n@NLparameter(model, y[i = 1:10] == 2 * i)\nJuMP.value(y[9])\n\n# output\n18.0\n\n\n\n\n\n"
},

{
    "location": "nlp.html#JuMP.value-Tuple{JuMP.NonlinearParameter}",
    "page": "Nonlinear Modeling",
    "title": "JuMP.value",
    "category": "method",
    "text": "value(p::NonlinearParameter)\n\nReturn the current value stored in the nonlinear parameter p.\n\nExample\n\nmodel = Model()\n@NLparameter(model, p == 10)\nJuMP.value(p)\n\n# output\n10.0\n\n\n\n\n\n"
},

{
    "location": "nlp.html#JuMP.set_value-Tuple{JuMP.NonlinearParameter,Number}",
    "page": "Nonlinear Modeling",
    "title": "JuMP.set_value",
    "category": "method",
    "text": "set_value(p::NonlinearParameter, v::Number)\n\nStore the value v in the nonlinear parameter p.\n\nExample\n\nmodel = Model()\n@NLparameter(model, p == 0)\nJuMP.set_value(p, 5)\nJuMP.value(p)\n\n# output\n5.0\n\n\n\n\n\n"
},

{
    "location": "nlp.html#Nonlinear-Parameters-1",
    "page": "Nonlinear Modeling",
    "title": "Nonlinear Parameters",
    "category": "section",
    "text": "For nonlinear models only, JuMP offers a syntax for explicit \"parameter\" objects which can be used to modify a model in-place just by updating the value of the parameter. Nonlinear parameters are declared by using the @NLparameter macro and may be indexed by arbitrary sets analogously to JuMP variables and expressions. The initial value of the parameter must be provided on the right-hand side of the == sign. There is no anonymous syntax for creating parameters.@NLparameterYou may use value and set_value to query or update the value of a parameter.value(::JuMP.NonlinearParameter)\nset_value(::JuMP.NonlinearParameter, ::Number)Nonlinear parameters can be used within nonlinear expressions only:@NLparameter(model, x == 10)\n@variable(model, z)\n@objective(model, Max, x * z)               # Error: x is a nonlinear parameter.\n@NLobjective(model, Max, x * z)             # Ok.\n@expression(model, my_expr, x * z ^ 2)      # Error: x is a nonlinear parameter.\n@NLexpression(model, my_nl_expr, x * z ^ 2) # Ok.Nonlinear parameters are useful when solving nonlinear models in a sequence:using Ipopt\nmodel = Model(with_optimizer(Ipopt.Optimizer))\n@variable(model, z)\n@NLparameter(model, x == 1.0)\n@NLobjective(model, Min, (z - x) ^ 2)\nJuMP.optimize!(model)\nJuMP.result_value(z) # Equals 1.0.\n\n# Now, update the value of x to solve a different problem.\nJuMP.set_value(x, 5.0)\nJuMP.optimize!(model)\nJuMP.result_value(z) # Equals 5.0Using nonlinear parameters can be faster than creating a new model from scratch with updated data because JuMP is able to avoid repeating a number of steps in processing the model before handing it off to the solver."
},

{
    "location": "nlp.html#User-defined-Functions-1",
    "page": "Nonlinear Modeling",
    "title": "User-defined Functions",
    "category": "section",
    "text": "JuMP\'s library of recognized univariate functions is derived from the Calculus.jl package. If you encounter a standard special function not currently supported by JuMP, consider contributing to the list of derivative rules there. In addition to this built-in list of functions, it is possible to register custom (user-defined) nonlinear functions to use within nonlinear expressions. JuMP does not support black-box optimization, so all user-defined functions must provide derivatives in some form. Fortunately, JuMP supports automatic differentiation of user-defined functions, a feature to our knowledge not available in any comparable modeling systems.Automatic differentiation is not finite differencing. JuMP\'s automatically computed derivatives are not subject to approximation error.JuMP uses ForwardDiff.jl to perform automatic differentiation; see the ForwardDiff.jl documentation for a description of how to write a function suitable for automatic differentiation. The general guideline is to write code that is generic with respect to the number type; don\'t assume that the input to the function is Float64. To register a user-defined function with derivatives computed by automatic differentiation, use the JuMP.register method as in the following example:my_square(x) = x ^ 2\nmy_f(x,y) = (x - 1) ^ 2 + (y - 2) ^ 2\n\nmodel = Model()\n\nJuMP.register(model, :my_f, 2, my_f, autodiff=true)\nJuMP.register(model, :my_square, 1, my_square, autodiff=true)\n\n@variable(model, x[1:2] >= 0.5)\n@NLobjective(model, Min, my_f(x[1], my_square(x[2])))The above code creates a JuMP model with the objective function (x[1] - 1) ^ 2 + (x[2] ^ 2 - 2) ^ 2. The first argument to JuMP.register the model for which the functions are registered. The second argument is a Julia symbol object which serves as the name of the user-defined function in JuMP expressions; the JuMP name need not be the same as the name of the corresponding Julia method. The third argument specifies how many arguments the function takes. The fourth argument is the name of the Julia method which computes the function, and autodiff=true instructs JuMP to compute exact gradients automatically.note: Note\nAll arguments to user-defined functions are scalars, not vectors. To define a function which takes a large number of arguments, you may use the splatting syntax f(x...) = ....Forward-mode automatic differentiation as implemented by ForwardDiff.jl has a computational cost that scales linearly with the number of input dimensions. As such, it is not the most efficient way to compute gradients of user-defined functions if the number of input arguments is large. In this case, users may want to provide their own routines for evaluating gradients. The more general syntax for JuMP.register which accepts user-provided derivative evaluation routines is:JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function,\n              ∇f::Function, ∇²f::Function)The input differs for functions which take a single input argument and functions which take more than one. For univariate functions, the derivative evaluation routines should return a number which represents the first and second-order derivatives respectively. For multivariate functions, the derivative evaluation routines will be passed a gradient vector which they must explicitly fill. Second-order derivatives of multivariate functions are not currently supported; this argument should be omitted. The following example sets up the same optimization problem as before, but now we explicitly provide evaluation routines for the user-defined functions:my_square(x) = x ^ 2\nmy_square_prime(x) = 2x\nmy_square_prime_prime(x) = 2\n\nmy_f(x, y) = (x - 1) ^ 2 + (y - 2) ^ 2\nfunction ∇f(g, x, y)\n    g[1] = 2 * (x - 1)\n    g[2] = 2 * (y - 2)\nend\n\nmodel = Model()\n\nJuMP.register(model, :my_f, 2, my_f, ∇f)\nJuMP.register(model, :my_square, 1, my_square, my_square_prime,\n              my_square_prime_prime)\n\n@variable(model, x[1:2] >= 0.5)\n@NLobjective(model, Min, my_f(x[1], my_square(x[2])))"
},

{
    "location": "nlp.html#Factors-affecting-solution-time-1",
    "page": "Nonlinear Modeling",
    "title": "Factors affecting solution time",
    "category": "section",
    "text": "The execution time when solving a nonlinear programming problem can be divided into two parts, the time spent in the optimization algorithm (the solver) and the time spent evaluating the nonlinear functions and corresponding derivatives. Ipopt explicitly displays these two timings in its output, for example:Total CPU secs in IPOPT (w/o function evaluations)   =      7.412\nTotal CPU secs in NLP function evaluations           =      2.083For Ipopt in particular, one can improve the performance by installing advanced sparse linear algebra packages, see Installation Guide. For other solvers, see their respective documentation for performance tips.The function evaluation time, on the other hand, is the responsibility of the modeling language. JuMP computes derivatives by using reverse-mode automatic differentiation with graph coloring methods for exploiting sparsity of the Hessian matrix [1]. As a conservative bound, JuMP\'s performance here currently may be expected to be within a factor of 5 of AMPL\'s."
},

{
    "location": "nlp.html#Querying-derivatives-from-a-JuMP-model-1",
    "page": "Nonlinear Modeling",
    "title": "Querying derivatives from a JuMP model",
    "category": "section",
    "text": "TODO: This section is out of date and hasn\'t been updated yet for MOI.For some advanced use cases, one may want to directly query the derivatives of a JuMP model instead of handing the problem off to a solver. Internally, JuMP implements the AbstractNLPEvaluator interface from MathOptInterface. To obtain an NLP evaluator object from a JuMP model, use JuMP.NLPEvaluator. JuMP.index returns the MOI.VariableIndex corresponding to a JuMP variable. MOI.VariableIndex itself is a type-safe wrapper for Int64 (stored in the value field.)For example:MOI = JuMP.MathOptInterface\nraw_index(v::MOI.VariableIndex) = v.value\nmodel = Model()\n@variable(model, x)\n@variable(model, y)\n@NLobjective(model, Min, sin(x) + sin(y))\nvalues = zeros(2)\nx_index = raw_index(JuMP.index(x))\ny_index = raw_index(JuMP.index(y))\nvalues[x_index] = 2.0\nvalues[y_index] = 3.0\nd = JuMP.NLPEvaluator(model)\nMOI.initialize(d, [:Grad])\nMOI.eval_objective(d, values) # == sin(2.0) + sin(3.0)\n\n# output\n1.0504174348855488∇f = zeros(2)\nMOI.eval_objective_gradient(d, ∇f, values)\n(∇f[x_index], ∇f[y_index]) # == (cos(2.0), cos(3.0))\n\n# output\n(-0.4161468365471424, -0.9899924966004454)Only nonlinear constraints (those added with @NLconstraint), and nonlinear objectives (added with @NLobjective) exist in the scope of the NLPEvaluator. The NLPEvaluator does not evaluate derivatives of linear or quadratic constraints or objectives. The index method applied to a nonlinear constraint reference object returns its index as a NonlinearConstraintIndex. The value field of NonlinearConstraintIndex stores the raw integer index. For example:julia> model = Model();\n\njulia> @variable(model, x);\n\njulia> @NLconstraint(model, cons1, sin(x) <= 1);\n\njulia> @NLconstraint(model, cons2, x + 5 == 10);\n\njulia> typeof(cons1)\nConstraintRef{Model,JuMP.NonlinearConstraintIndex,JuMP.ScalarShape}\n\njulia> JuMP.index(cons1)\nJuMP.NonlinearConstraintIndex(1)\n\njulia> JuMP.index(cons2)\nJuMP.NonlinearConstraintIndex(2)TODO: Provide a link for how to access the linear and quadratic parts of the model.Note that for one-sided nonlinear constraints, JuMP subtracts any values on the right-hand side when computing expressions. In other words, one-sided nonlinear constraints are always transformed to have a right-hand side of zero.This method of querying derivatives directly from a JuMP model is convenient for interacting with the model in a structured way, e.g., for accessing derivatives of specific variables. For example, in statistical maximum likelihood estimation problems, one is often interested in the Hessian matrix at the optimal solution, which can be queried using the JuMP.NLPEvaluator."
},

{
    "location": "nlp.html#Raw-expression-input-1",
    "page": "Nonlinear Modeling",
    "title": "Raw expression input",
    "category": "section",
    "text": "In addition to the @NLobjective and @NLconstraint macros, it is also possible to provide Julia Expr objects directly by using JuMP.set_NL_objective and JuMP.add_NL_constraint. This input form may be useful if the expressions are generated programmatically. JuMP variables should be spliced into the expression object. For example:@variable(model, 1 <= x[i = 1:4] <= 5)\nJuMP.set_NL_objective(model, :Min, :($(x[1])*$(x[4])*($(x[1])+$(x[2])+$(x[3])) + $(x[3])))\nJuMP.add_NL_constraint(model, :($(x[1])*$(x[2])*$(x[3])*$(x[4]) >= 25))\n\n# Equivalent form using traditional JuMP macros:\n@NLobjective(model, Min, x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3])\n@NLconstraint(model, x[1] * x[2] * x[3] * x[4] >= 25)See the Julia documentation for more examples and description of Julia expressions.[1]: Dunning, Huchette, and Lubin, \"JuMP: A Modeling Language for Mathematical Optimization\", arXiv."
},

{
    "location": "style.html#",
    "page": "Style Guide",
    "title": "Style Guide",
    "category": "page",
    "text": ""
},

{
    "location": "style.html#Style-guide-and-design-principles-1",
    "page": "Style Guide",
    "title": "Style guide and design principles",
    "category": "section",
    "text": ""
},

{
    "location": "style.html#Style-guide-1",
    "page": "Style Guide",
    "title": "Style guide",
    "category": "section",
    "text": "This section describes the coding style rules that apply to JuMP code and that we recommend for JuMP models and surrounding Julia code. The motivations for a style guide include:conveying best practices for writing readable and maintainable code\nreducing the amount of time spent on bike-shedding by establishing basic naming and formatting conventions\nlowering the barrier for new contributors by codifying the existing practices (e.g., you can be more confident your code will pass review if you follow the style guide)In some cases, the JuMP style guide diverges from the Julia style guide. All such cases will be explicitly noted and justified.The JuMP style guide adopts many recommendations from the Google style guides.info: Info\nThe style guide is always a work in progress, and not all JuMP code follows the rules. When modifying JuMP, please fix the style violations of the surrounding code (i.e., leave the code tidier than when you started). If large changes are needed, consider separating them into another PR."
},

{
    "location": "style.html#Formatting-1",
    "page": "Style Guide",
    "title": "Formatting",
    "category": "section",
    "text": "Julia unfortunately does not have an autoformatting tool like gofmt. Until a reliable autoformatting tool is available, we adopt the following conventions."
},

{
    "location": "style.html#Whitespace-1",
    "page": "Style Guide",
    "title": "Whitespace",
    "category": "section",
    "text": "Julia is mostly insensitive to whitespace characters within lines. For consistency:Use spaces between binary operators\nUse a single space after commas and semicolons\nDo not use extra spaces for unary operators, parentheses, or braces\nIndent within new blocks (except module) using 4 spacesGood:f(x, y) = [3 * dot(x, y); x\']Bad:f(x,y) = [ 3*dot(x,y) ; x\' ]Good:module Foo\n\nfunction f(x)\n    return x + 1\nend\n\nend # module Foo"
},

{
    "location": "style.html#TODO:-Line-breaks-1",
    "page": "Style Guide",
    "title": "TODO: Line breaks",
    "category": "section",
    "text": ""
},

{
    "location": "style.html#Syntax-1",
    "page": "Style Guide",
    "title": "Syntax",
    "category": "section",
    "text": "Julia sometimes provides equivalent syntax to express the same basic operation. We discuss these cases below."
},

{
    "location": "style.html#for-loops-1",
    "page": "Style Guide",
    "title": "for loops",
    "category": "section",
    "text": "Julia allows both for x = 1:N and for x in 1:N. Always prefer to use in over =, because in generalizes better to other index sets like for x in eachindex(A)."
},

{
    "location": "style.html#Empty-vectors-1",
    "page": "Style Guide",
    "title": "Empty vectors",
    "category": "section",
    "text": "For a type T, T[] and Vector{T}() are equivalent ways to create an empty vector with element type T. Prefer T[] because it is more concise."
},

{
    "location": "style.html#Trailing-periods-in-floating-point-constants-1",
    "page": "Style Guide",
    "title": "Trailing periods in floating-point constants",
    "category": "section",
    "text": "Both 1.0 and 1. create a Float64 with value 1.0. Prefer 1.0 over 1. because it is more easily distinguished from the integer constant 1."
},

{
    "location": "style.html#Naming-1",
    "page": "Style Guide",
    "title": "Naming",
    "category": "section",
    "text": "module SomeModule end\nfunction some_function end\nconst SOME_CONSTANT = ...\nstruct SomeStruct end\nsome_local_variable = ...\nsome_file.jl # Except for ModuleName.jl."
},

{
    "location": "style.html#Use-of-underscores-1",
    "page": "Style Guide",
    "title": "Use of underscores",
    "category": "section",
    "text": "The Julia style guide recommends avoiding underscores \"when readable\", for example, haskey, isequal, remotecall, and remotecall_fetch. This convention creates the potential for unnecessary bikeshedding and also forces the user to recall the presence/absence of an underscore, e.g., \"was that argument named basename or base_name?\". For consistency, always use underscores in variable names and function names to separate words."
},

{
    "location": "style.html#Use-of-!-1",
    "page": "Style Guide",
    "title": "Use of !",
    "category": "section",
    "text": "Julia has a convention of appending ! to a function name if the function modifies its arguments. We recommend to:Omit ! when the name itself makes it clear that modification is taking place, e.g., add_constraint and set_name. We depart from the Julia style guide because ! does not provide a reader with any additional information in this case, and adherence to this convention is not uniform even in base Julia itself (consider Base.println and Base.finalize).\nUse ! in all other cases. In particular it can be used to distinguish between modifying and non-modifying variants of the same function like scale and scale!.Note that ! is not a self-documenting feature because it is still ambiguous which arguments are modified when multiple arguments are present. Be sure to document which arguments are modified in the method\'s docstring.See also the Julia style guide recommendations for ordering of function arguments."
},

{
    "location": "style.html#Abbreviations-1",
    "page": "Style Guide",
    "title": "Abbreviations",
    "category": "section",
    "text": "Abbreviate names to make the code more readable, not to save typing. Don\'t arbitrarily delete letters from a word to abbreviate it (e.g., indx). Use abbreviations consistently within a body of code (e.g., do not mix con and constr, idx and indx).Common abbreviations:num for numberTODO: add more"
},

{
    "location": "style.html#Miscellaneous-1",
    "page": "Style Guide",
    "title": "Miscellaneous",
    "category": "section",
    "text": "(TODO: Rethink categories.)"
},

{
    "location": "style.html#User-facing-MethodError-1",
    "page": "Style Guide",
    "title": "User-facing MethodError",
    "category": "section",
    "text": "Specifying argument types for methods is mostly optional in Julia, which means that it\'s possible to find out that you are working with unexpected types deep in the call chain. Avoid this situation or handle it with a helpful error message. A user should see a MethodError only for methods that they called directly.Bad:internal_function(x::Integer) = x + 1\n# The user sees a MethodError for internal_function when calling\n# public_function(\"a string\"). This is not very helpful.\npublic_function(x) = internal_function(x)Good:internal_function(x::Integer) = x + 1\n# The user sees a MethodError for public_function when calling\n# public_function(\"a string\"). This is easy to understand.\npublic_function(x::Integer) = internal_function(x)If it is hard to provide an error message at the top of the call chain, then the following pattern is also ok:internal_function(x::Integer) = x + 1\nfunction internal_function(x)\n    error(\"Internal error. This probably means that you called \" *\n          \"public_function() with the wrong type.\")\nend\npublic_function(x) = internal_function(x)"
},

{
    "location": "style.html#Design-principles-1",
    "page": "Style Guide",
    "title": "Design principles",
    "category": "section",
    "text": "TODO: How to structure and test large JuMP models, libraries that use JuMP.For how to write a solver, see MOI."
},

{
    "location": "extensions.html#",
    "page": "Extensions",
    "title": "Extensions",
    "category": "page",
    "text": ""
},

{
    "location": "extensions.html#Extending-JuMP-1",
    "page": "Extensions",
    "title": "Extending JuMP",
    "category": "section",
    "text": "TODO: How to extend JuMP: discussion on different ways to build on top of JuMP. How to extend JuMP\'s macros and how to avoid doing this."
},

{
    "location": "updating.html#",
    "page": "Updating Guide",
    "title": "Updating Guide",
    "category": "page",
    "text": ""
},

{
    "location": "updating.html#Updating-Guide-1",
    "page": "Updating Guide",
    "title": "Updating Guide",
    "category": "section",
    "text": "DRAFT: See also NEWS.md for updates between releases."
},

{
    "location": "updating.html#Updating-from-JuMP-0.18-to-JuMP-0.19-1",
    "page": "Updating Guide",
    "title": "Updating from JuMP 0.18 to JuMP 0.19",
    "category": "section",
    "text": "TODO: XX% of JuMP\'s source code changed between JuMP 0.18 and JuMP 0.19. Switched from MPB to MOI. Explain what broke and how to update."
},

{
    "location": "howdoi.html#",
    "page": "How do I ...? (FAQ)",
    "title": "How do I ...? (FAQ)",
    "category": "page",
    "text": ""
},

{
    "location": "howdoi.html#How-do-I-...?-(FAQ)-1",
    "page": "How do I ...? (FAQ)",
    "title": "How do I ...? (FAQ)",
    "category": "section",
    "text": "Q: I\'m using a solver that supports warm-starts for integer problems. How do I communicate the initial solution to the solver?TODO: use start= in @variable or the VariablePrimalStart() attribute.Q: How can I suppress output in a solver-independent way?TODO: Update answer for JuMP 0.19.A: JuMP does not currently support generic parameters independent of the chosen solver object. To suppress output with Gurobi, for example, one would saym = Model(solver=GurobiSolver(OutputFlag=0))When a solver is not specified, i.e., the model is created with m = Model(),  there\'s no option to suppress output. A workaround is to redirect STDOUT before  and after the call to solve(m):TT = STDOUT # save original STDOUT stream\nredirect_stdout()\nsolve(m)\nredirect_stdout(TT) # restore STDOUT"
},

]}
