var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "#![JuMP-logo](assets/jump-logo-with-text.svg)-1",
    "page": "Introduction",
    "title": "(Image: JuMP logo)",
    "category": "section",
    "text": "(Image: Powered by NumFOCUS)# These comments do not display in the HTML output.\n# See https://github.com/JuliaDocs/Documenter.jl/issues/674.warning: Warning\nBetween versions 0.18 and 0.19, JuMP underwent a major transition in its underlying solver abstraction API, from MathProgBase to MathOptInterface. See NEWS.md for a comprehensive list of changes between the two versions, many of which are breaking. This documentation is for JuMP/MathOptInterface. For the documentation of JuMP 0.18, see here.JuMP is a domain-specific modeling language for mathematical optimization embedded in Julia. It currently supports a number of open-source and commercial solvers (see below) for a variety of problem classes, including linear programming, mixed-integer programming, second-order conic programming, semidefinite programming, and nonlinear programming. JuMP\'s features include:User friendliness\nSyntax that mimics natural mathematical expressions.\nComplete documentation (WIP!)\nSpeed\nBenchmarking has shown that JuMP can create problems at similar speeds   to special-purpose modeling languages such as   AMPL.\nJuMP communicates with most solvers in memory, avoiding the need to   write intermediary files.\nSolver independence\nJuMP uses a generic solver-independent interface provided by the   MathOptInterface   package, making it easy to change between a number of open-source and   commercial optimization software packages (\"solvers\").\nCurrently supported solvers include   Artelys Knitro,   Bonmin,   Cbc,   Clp,   Couenne,   CPLEX,   ECOS,   FICO Xpress,   GLPK,   Gurobi,   Ipopt,   MOSEK,   NLopt, and   SCS.\nAccess to advanced algorithmic techniques\nIncluding efficient LP re-solves which previously required using   solver-specific and/or low-level C++ libraries.\nEase of embedding\nJuMP itself is written purely in Julia. Solvers are the only binary   dependencies.\nBeing embedded in a general-purpose programming language makes it easy   to solve optimization problems as part of a larger workflow (e.g.,   inside a simulation, behind a web server, or as a subproblem in a   decomposition algorithm).\nAs a trade-off, JuMP\'s syntax is constrained by the syntax available   in Julia.\nJuMP is MPL licensed, meaning that   it can be embedded in commercial software that complies with the terms   of the license.Although JuMP has not reached version 1.0 yet, the releases are stable enough for everyday use and are being used in a number of research projects and neat applications by a growing community of users who are early adopters. JuMP remains under active development, and we welcome your feedback, suggestions, and bug reports."
},

{
    "location": "#Contents-1",
    "page": "Introduction",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\n    \"installation.md\",\n    \"quickstart.md\",\n    \"variables.md\",\n    \"expressions.md\",\n    \"objective.md\",\n    \"constraints.md\",\n    \"containers.md\",\n    \"solvers.md\",\n    \"solutions.md\",\n    \"nlp.md\",\n    \"style.md\",\n    \"extensions.md\"\n]\nDepth = 2"
},

{
    "location": "#Citing-JuMP-1",
    "page": "Introduction",
    "title": "Citing JuMP",
    "category": "section",
    "text": "If you find JuMP useful in your work, we kindly request that you cite the following paper (pdf):@article{DunningHuchetteLubin2017,\nauthor = {Iain Dunning and Joey Huchette and Miles Lubin},\ntitle = {JuMP: A Modeling Language for Mathematical Optimization},\njournal = {SIAM Review},\nvolume = {59},\nnumber = {2},\npages = {295-320},\nyear = {2017},\ndoi = {10.1137/15M1020575},\n}For an earlier work where we presented a prototype implementation of JuMP, see here:@article{LubinDunningIJOC,\nauthor = {Miles Lubin and Iain Dunning},\ntitle = {Computing in Operations Research Using Julia},\njournal = {INFORMS Journal on Computing},\nvolume = {27},\nnumber = {2},\npages = {238-248},\nyear = {2015},\ndoi = {10.1287/ijoc.2014.0623},\n}A preprint of this paper is freely available.(Image: NumFOCUS logo)JuMP is a fiscally sponsored project of NumFOCUS, a nonprofit dedicated to supporting the open source scientific computing community.You can support JuMP by donating. JuMP\'s largest expense is the annual JuMP-dev workshop. Donations will help us provide travel support for JuMP-dev attendees and take advantage of other opportunities that arise to support JuMP development."
},

{
    "location": "installation/#",
    "page": "Installation Guide",
    "title": "Installation Guide",
    "category": "page",
    "text": ""
},

{
    "location": "installation/#Installation-Guide-1",
    "page": "Installation Guide",
    "title": "Installation Guide",
    "category": "section",
    "text": "JuMP is a package for Julia. To use JuMP, first download and install Julia or open up a remote notebook at JuliaBox or similar services.This version of JuMP is compatible with Julia 0.7 and later. The following instructions are based on Julia 1.0.From Julia, JuMP is installed by using the built-in package manager:import Pkg\nPkg.add(\"JuMP\")"
},

{
    "location": "installation/#Getting-Solvers-1",
    "page": "Installation Guide",
    "title": "Getting Solvers",
    "category": "section",
    "text": "JuMP depends on solvers to solve optimization problems. Most solvers are not written in Julia, and some require commercial licenses to use, so installation is often more complex. We list below the currently available solvers.note: Note\nThis list is open for new contributions. See also Interacting with solvers and the MathOptInterface docs for more details on how JuMP interacts with solvers. Please get in touch with any questions about connecting new solvers with JuMP.Solver Julia Package License Supports\nCbc Cbc.jl EPL MILP\nClp Clp.jl EPL LP\nCPLEX CPLEX.jl Comm. LP, MILP, SOCP, MISOCP\nCSDP CSDP.jl EPL LP, SDP\nECOS ECOS.jl GPL LP, SOCP\nFICO Xpress Xpress.jl Comm. LP, MILP, SOCP, MISOCP\nGLPK GLPK.jl GPL LP, MILP\nGurobi Gurobi.jl Comm. LP, MILP, SOCP, MISOCP\nIpopt Ipopt.jl EPL LP, QP, NLP\nMOSEK MathOptInterfaceMosek.jl Comm. LP, MILP, SOCP, MISOCP, SDP\nOSQP OSQP.jl Apache LP, QP\nSCS SCS.jl MIT LP, SOCP, SDP\nSeDuMi SeDuMi.jl GPL LP, SOCP, SDPWhere:LP = Linear programming\nQP = Quadratic programming\nSOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)\nMILP = Mixed-integer linear programming\nNLP = Nonlinear programming\nMINLP = Mixed-integer nonlinear programming\nSDP = Semidefinite programming\nMISDP = Mixed-integer semidefinite programmingTo install Gurobi, for example, and use it with a JuMP model model, run:import Pkg\nPkg.add(\"Gurobi\")\nusing JuMP\nusing Gurobi\nmodel = Model(with_optimizer(Gurobi.Optimizer))Most packages follow the ModuleName.Optimizer naming convention, but exceptions may exist. See the corresponding Julia package README for more details on how to use the solver.# TODO: Discuss setting solver options.The following solvers were compatible with JuMP up to release 0.18 but are not yet compatible with the latest version because they do not implement the new MathOptInterface API:Artelys Knitro\nBARON\nBonmin and Couenne via AmplNLWriter.jl\nCDD\nNLopt\nPavito\nPajarito\nSCIPSolver-specific notes follow below."
},

{
    "location": "installation/#Artelys-Knitro-1",
    "page": "Installation Guide",
    "title": "Artelys Knitro",
    "category": "section",
    "text": "Requires a license. The KNITRO.jl interface currently supports only nonlinear problems."
},

{
    "location": "installation/#BARON-1",
    "page": "Installation Guide",
    "title": "BARON",
    "category": "section",
    "text": "Requires a license. A trial version is available for small problem instances."
},

{
    "location": "installation/#COIN-OR-Cbc-1",
    "page": "Installation Guide",
    "title": "COIN-OR Cbc",
    "category": "section",
    "text": "Cbc supports \"SOS\" constraints."
},

{
    "location": "installation/#CPLEX-1",
    "page": "Installation Guide",
    "title": "CPLEX",
    "category": "section",
    "text": "Requires a working installation of CPLEX with a license (free for faculty members and graduate teaching assistants). The interface requires using CPLEX as a shared library, which is unsupported by the CPLEX developers. Special installation steps are required on Mac OS. CPLEX supports \"SOS\" constraints."
},

{
    "location": "installation/#ECOS-1",
    "page": "Installation Guide",
    "title": "ECOS",
    "category": "section",
    "text": "ECOS can be used by JuMP to solve LPs and SOCPs. ECOS does not support general quadratic objectives or constraints, only second-order conic constraints specified by using the SecondOrderCone set."
},

{
    "location": "installation/#Gurobi-1",
    "page": "Installation Guide",
    "title": "Gurobi",
    "category": "section",
    "text": "Requires a working installation of Gurobi with an activated license (free for academic use). Gurobi supports \"SOS\" constraints."
},

{
    "location": "installation/#FICO-Xpress-1",
    "page": "Installation Guide",
    "title": "FICO Xpress",
    "category": "section",
    "text": "Requires a working installation of Xpress with an active license (it is possible to get a license for academic use, see FICO Academic Partner Program). Supports SOCP and \"SOS\" constraints."
},

{
    "location": "installation/#MOSEK-1",
    "page": "Installation Guide",
    "title": "MOSEK",
    "category": "section",
    "text": "Requires a license (free for academic use). The Mosek interface is maintained by the Mosek team. (Thanks!)"
},

{
    "location": "installation/#SCS-1",
    "page": "Installation Guide",
    "title": "SCS",
    "category": "section",
    "text": "SCS can be used by JuMP to solve LPs and SOCPs, and SDPs. SCS is a first order solver and has low accuracy (10^4) by default; see the SCS.jl documentation for more information."
},

{
    "location": "quickstart/#",
    "page": "Quick Start Guide",
    "title": "Quick Start Guide",
    "category": "page",
    "text": ""
},

{
    "location": "quickstart/#Quick-Start-Guide-1",
    "page": "Quick Start Guide",
    "title": "Quick Start Guide",
    "category": "section",
    "text": "This quick start guide will introduce the main concepts of JuMP. If you are familiar with another modeling language embedded in a high-level language such as PuLP (Python) or a solver-specific interface you will find most of this familiar. If you are coming from an AMPL or similar background, you may find some of the concepts novel but the general appearance will still be familiar.The example in this guide is deliberately kept simple. There are more complex examples in the JuMP/examples/ folder.Once JuMP is installed, to use JuMP in your programs, you just need to say:julia> using JuMPYou also need to include a Julia package which provides an appropriate solver. One such solver is GLPK.Optimizer, which is provided by the GLPK.jl package.julia> using GLPKSee Installation Guide for a list of other solvers you can use.Models are created with the Model() function. The with_optimizer syntax is used to specify the optimizer to be used:julia> model = Model(with_optimizer(GLPK.Optimizer))\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.note: Note\nThe term \"solver\" is used as a synonym for \"optimizer\". The convention in code, however, is to always use \"optimizer\", e.g., with_optimizer and GLPK.Optimizer.DocTestSetup = quote\n    # Using a mock optimizer removes the need to load a solver such as GLPK for\n    # building the documentation.\n    const MOI = JuMP.MathOptInterface\n    model = Model(with_optimizer(MOI.Utilities.MockOptimizer,\n                                 JuMP._MOIModel{Float64}(),\n                                 eval_objective_value = false,\n                                 eval_variable_constraint_dual = false))\nendnote: Note\nYour model doesn\'t have to be called model - it\'s just a name.The following commands will create two variables, x and y, with both lower and upper bounds. Note the first argument is our model model. These variables (x and y) are associated with this model and cannot be used in another model.julia> @variable(model, 0 <= x <= 2)\nx\n\njulia> @variable(model, 0 <= y <= 30)\nySee the Variables section for more information on creating variables, including the syntax for specifying different combinations of bounds, i.e., only lower bounds, only upper bounds, or no bounds.DocTestSetup = nothingNext we\'ll set our objective. Note again the model, so we know which model\'s objective we are setting! The objective sense, Max or Min, should be provided as the second argument. Note also that we don\'t have a multiplication * symbol between 5 and our variable x - Julia is smart enough to not need it! Feel free to use * if it makes you feel more comfortable, as we have done with 3 * y. (We have been intentionally inconsistent here to demonstrate different syntax; however, it is good practice to pick one way or the other consistently in your code.)julia> @objective(model, Max, 5x + 3 * y)\n5 x + 3 yAdding constraints is a lot like setting the objective. Here we create a less-than-or-equal-to constraint using <=, but we can also create equality constraints using == and greater-than-or-equal-to constraints with >=:julia> @constraint(model, con, 1x + 5y <= 3)\ncon : x + 5 y <= 3.0Note that in a similar manner to the @variable macro, we have named the constraint con. This will bind the constraint to the Julia variable con for later analysis.Models are solved with the JuMP.optimize! function:julia> optimize!(model)DocTestSetup = quote\n    # Now we load in the solution. Using a caching optimizer removes the need to\n    # load a solver such as GLPK for building the documentation.\n    mock = JuMP.backend(model).optimizer.model\n    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)\n    MOI.set(mock, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)\n    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)\n    MOI.set(mock, MOI.ResultCount(), 1)\n    MOI.set(mock, MOI.ObjectiveValue(), 10.6)\n    MOI.set(mock, MOI.VariablePrimal(), JuMP.optimizer_index(x), 2.0)\n    MOI.set(mock, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.2)\n    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(con), -0.6)\n    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(UpperBoundRef(x)), -4.4)\n    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(LowerBoundRef(y)), 0.0)\nendAfter the call to JuMP.optimize! has finished, we need to query what happened. The solve could terminate for a number of reasons. First, the solver might have found the optimal solution or proved that the problem is infeasible. However, it might also have run into numerical difficulties or terminated due to a setting such as a time limit. We can ask the solver why it stopped using the JuMP.termination_status function:julia> termination_status(model)\nOPTIMAL::TerminationStatusCode = 1In this case, GLPK returned OPTIMAL, this mean that it has found the optimal solution.DocTestSetup = nothingAs the solver found an optimal solution, we expect the solution returned to be a primal-dual pair of feasible solutions with zero duality gap. We can verify the primal and dual status as follows to confirm this:julia> primal_status(model)\nFEASIBLE_POINT::ResultStatusCode = 1\n\njulia> dual_status(model)\nFEASIBLE_POINT::ResultStatusCode = 1Note that the primal and dual status only inform that the primal and dual solutions are feasible and it is only because we verified that the termination status is OPTIMAL that we can conclude that they form an optimal solution.Finally, we can query the result of the optimization. First, we can query the objective value:julia> objective_value(model)\n10.6We can also query the primal result values of the x and y variables:julia> value(x)\n2.0\n\njulia> value(y)\n0.2We can also query the value of the dual variable associated with the constraint con (which we bound to a Julia variable when defining the constraint):julia> dual(con)\n-0.6info: Info\nSee the duality section for a discussion of the convention that JuMP uses for signs of duals.To query the dual variables associated with the variable bounds, things are a little trickier as we first need to obtain a reference to the constraint:julia> x_upper = UpperBoundRef(x)\nx <= 2.0\n\njulia> dual(x_upper)\n-4.4A similar process can be followed to obtain the dual of the lower bound constraint on y:julia> y_lower = LowerBoundRef(y)\ny >= 0.0\n\njulia> dual(y_lower)\n0.0"
},

{
    "location": "variables/#",
    "page": "Variables",
    "title": "Variables",
    "category": "page",
    "text": "CurrentModule = JuMP\nDocTestSetup = quote\n    using JuMP\nend"
},

{
    "location": "variables/#Variables-1",
    "page": "Variables",
    "title": "Variables",
    "category": "section",
    "text": ""
},

{
    "location": "variables/#What-is-a-JuMP-variable?-1",
    "page": "Variables",
    "title": "What is a JuMP variable?",
    "category": "section",
    "text": "The term variable in mathematical optimization has many meanings. Here, we distinguish between the following three types of variables:optimization variables, which are the mathematical x in the problem maxf_0(x)  f_i(x) in S_i.\nJulia variables, which are bindings between a name and a value, for example x = 1. (See here for the Julia docs.)\nJuMP variables, which are instances of the VariableRef struct defined by JuMP that contains a reference to an optimization variable in a model. (Extra for experts: the VariableRef struct is a thin wrapper around a MOI.VariableIndex, and also contains a reference to the JuMP model.)To illustrate these three types of variables, consider the following JuMP code (the full syntax is explained below):julia> model = Model()\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\n\njulia> @variable(model, x[1:2])\n2-element Array{VariableRef,1}:\n x[1]\n x[2]This code does three things:It adds two optimization variables to model.\nIt creates two JuMP variables that act as references to those optimization variables.\nIt binds those JuMP variables as a vector with two elements to the Julia variable x.To reduce confusion, we will attempt, where possible, to always refer to variables with their corresponding prefix.JuMP variables can have attributes, such as names or an initial primal start value. We illustrate the name attribute in the following example:julia> @variable(model, y, base_name=\"decision variable\")\ndecision variableThis code does four things:It adds one optimization variable to model.\nIt creates one JuMP variable that acts as a reference to that optimization variable.\nIt binds the JuMP variable to the Julia variable y.\nIt tells JuMP that the name attribute of this JuMP variable is \"decision variable\". JuMP uses the value of base_name when it has to print the variable as a string.For example, when we print y at the REPL we get:julia> y\ndecision variableBecause y is a Julia variable, we can bind it to a different value. For example, if we write:julia> y = 1\n1y is no longer a binding to a JuMP variable. This does not mean that the JuMP variable has been destroyed. It still exists and is still a reference to the same optimization variable. The binding can be reset by querying the model for the symbol as it was written in the @variable macro. For example:julia> model[:y]\ndecision variableThis act of looking up the JuMP variable by using the symbol is most useful when composing JuMP models across multiple functions, as illustrated by the following example:function add_component_to_model(model::JuMP.Model)\n    x = model[:x]\n    # ... code that uses x\nend\nfunction build_model()\n    model = Model()\n    @variable(model, x)\n    add_component_to_model(model)\nendNow that we understand the difference between optimization, JuMP, and Julia variables, we can introduce more of the functionality of the @variable macro."
},

{
    "location": "variables/#Variable-bounds-1",
    "page": "Variables",
    "title": "Variable bounds",
    "category": "section",
    "text": "We have already seen the basic usage of the @variable macro. The next extension is to add lower- and upper-bounds to each optimization variable. This can be done as follows:julia> @variable(model, x_free)\nx_free\n\njulia> @variable(model, x_lower >= 0)\nx_lower\n\njulia> @variable(model, x_upper <= 1)\nx_upper\n\njulia> @variable(model, 2 <= x_interval <= 3)\nx_interval\n\njulia> @variable(model, x_fixed == 4)\nx_fixedIn the above examples, x_free represents an unbounded optimization variable, x_lower represents an optimization variable with a lower bound and so forth.note: Note\nWhen creating a variable with only a lower-bound or an upper-bound, and the value of the bound is not a numeric literal, the name must appear on the left-hand side. Putting the name on the right-hand side will result in an error. For example:@variable(model, 1 <= x)  # works\na = 1\n@variable(model, a <= x)  # errorsWe can query whether an optimization variable has a lower- or upper-bound via the has_lower_bound and has_upper_bound functions. For example:julia> has_lower_bound(x_free)\nfalse\n\njulia> has_upper_bound(x_upper)\ntrueIf a variable has a lower or upper bound, we can query the value of it via the lower_bound and upper_bound functions. For example:julia> lower_bound(x_interval)\n2.0\n\njulia> upper_bound(x_interval)\n3.0Querying the value of a bound that does not exist will result in an error.Instead of using the <= and >= syntax, we can also use the lower_bound and upper_bound keyword arguments. For example:julia> @variable(model, x, lower_bound=1, upper_bound=2)\nx\n\njulia> lower_bound(x)\n1.0Another option is to use the set_lower_bound and set_upper_bound functions. These can also be used to modify an existing variable bound. For example:julia> @variable(model, x >= 1)\nx\n\njulia> lower_bound(x)\n1.0\n\njulia> set_lower_bound(x, 2)\n\njulia> lower_bound(x)\n2.0We can delete variable bounds using delete_lower_bound and delete_upper_bound:julia> @variable(model, 1 <= x <= 2)\nx\n\njulia> lower_bound(x)\n1.0\n\njulia> delete_lower_bound(x)\n\njulia> has_lower_bound(x)\nfalse\n\njulia> upper_bound(x)\n2.0\n\njulia> delete_upper_bound(x)\n\njulia> has_upper_bound(x)\nfalseIn addition to upper and lower bounds, JuMP variables can also be fixed to a value using fix. See also is_fixed, fix_value, and unfix.julia> @variable(model, x == 1)\nx\n\njulia> is_fixed(x)\ntrue\n\njulia> fix_value(x)\n1.0\n\njulia> unfix(x)\n\njulia> is_fixed(x)\nfalseFixing a variable with existing bounds will throw an error. To delete the bounds prior to fixing, use fix(variable, value; force = true).julia> @variable(model, x >= 1)\nx\n\njulia> fix(x, 2)\nERROR: Unable to fix x to 2 because it has existing variable bounds. Consider calling `JuMP.fix(variable, value; force=true)` which will delete existing bounds before fixing the variable.\n\njulia> fix(x, 2; force = true)\n\n\njulia> fix_value(x)\n2.0"
},

{
    "location": "variables/#JuMP.name-Tuple{VariableRef}",
    "page": "Variables",
    "title": "JuMP.name",
    "category": "method",
    "text": "name(v::VariableRef)::String\n\nGet a variable\'s name attribute.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.set_name-Tuple{VariableRef,String}",
    "page": "Variables",
    "title": "JuMP.set_name",
    "category": "method",
    "text": "set_name(v::VariableRef, s::AbstractString)\n\nSet a variable\'s name attribute.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.variable_by_name",
    "page": "Variables",
    "title": "JuMP.variable_by_name",
    "category": "function",
    "text": "variable_by_name(model::AbstractModel,\n                 name::String)::Union{AbstractVariableRef, Nothing}\n\nReturns the reference of the variable with name attribute name or Nothing if no variable has this name attribute. Throws an error if several variables have name as their name attribute.\n\njulia> model = Model()\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\n\njulia> @variable(model, x)\nx\n\njulia> variable_by_name(model, \"x\")\nx\n\njulia> @variable(model, base_name=\"x\")\nx\n\njulia> variable_by_name(model, \"x\")\nERROR: Multiple variables have the name x.\nStacktrace:\n [1] error(::String) at ./error.jl:33\n [2] get(::JuMP._MOIModel{Float64}, ::Type{MathOptInterface.VariableIndex}, ::String) at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/model.jl:222\n [3] get at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/universalfallback.jl:201 [inlined]\n [4] get(::MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer,MathOptInterface.Utilities.UniversalFallback{JuMP._MOIModel{Float64}}}, ::Type{MathOptInterface.VariableIndex}, ::String) at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/cachingoptimizer.jl:490\n [5] variable_by_name(::Model, ::String) at /home/blegat/.julia/dev/JuMP/src/variables.jl:268\n [6] top-level scope at none:0\n\njulia> var = @variable(model, base_name=\"y\")\ny\n\njulia> variable_by_name(model, \"y\")\ny\n\njulia> set_name(var, \"z\")\n\njulia> variable_by_name(model, \"y\")\n\njulia> variable_by_name(model, \"z\")\nz\n\njulia> @variable(model, u[1:2])\n2-element Array{VariableRef,1}:\n u[1]\n u[2]\n\njulia> variable_by_name(model, \"u[2]\")\nu[2]\n\n\n\n\n\n"
},

{
    "location": "variables/#Variable-names-1",
    "page": "Variables",
    "title": "Variable names",
    "category": "section",
    "text": "The name, i.e. the value of the MOI.VariableName attribute, of a variable can be obtained by JuMP.name(::JuMP.VariableRef) and set by JuMP.set_name(::JuMP.VariableRef, ::String).name(::JuMP.VariableRef)\nset_name(::JuMP.VariableRef, ::String)The variable can also be retrieved from its name using JuMP.variable_by_name.variable_by_name"
},

{
    "location": "variables/#Variable-containers-1",
    "page": "Variables",
    "title": "Variable containers",
    "category": "section",
    "text": "In the examples above, we have mostly created scalar variables. By scalar, we mean that the Julia variable is bound to exactly one JuMP variable. However, it is often useful to create collections of JuMP variables inside more complicated data structures.JuMP provides a mechanism for creating three types of these data structures, which we refer to as containers. The three types are Arrays, DenseAxisArrays, and SparseAxisArrays. We explain each of these in the following."
},

{
    "location": "variables/#Arrays-1",
    "page": "Variables",
    "title": "Arrays",
    "category": "section",
    "text": "We have already seen the creation of an array of JuMP variables with the x[1:2] syntax. This can naturally be extended to create multi-dimensional arrays of JuMP variables. For example:julia> @variable(model, x[1:2, 1:2])\n2×2 Array{VariableRef,2}:\n x[1,1]  x[1,2]\n x[2,1]  x[2,2]Arrays of JuMP variables can be indexed and sliced as follows:julia> x[1, 2]\nx[1,2]\n\njulia> x[2, :]\n2-element Array{VariableRef,1}:\n x[2,1]\n x[2,2]Variable bounds can depend upon the indices:julia> @variable(model, x[i=1:2, j=1:2] >= 2i + j)\n2×2 Array{VariableRef,2}:\n x[1,1]  x[1,2]\n x[2,1]  x[2,2]\n\njulia> lower_bound.(x)\n2×2 Array{Float64,2}:\n 3.0  4.0\n 5.0  6.0JuMP will form an Array of JuMP variables when it can determine at compile time that the indices are one-based integer ranges. Therefore x[1:b] will create an Array of JuMP variables, but x[a:b] will not. If JuMP cannot determine that the indices are one-based integer ranges (e.g., in the case of x[a:b]), JuMP will create a DenseAxisArray instead."
},

{
    "location": "variables/#variable_jump_arrays-1",
    "page": "Variables",
    "title": "DenseAxisArrays",
    "category": "section",
    "text": "We often want to create arrays where the indices are not one-based integer ranges. For example, we may want to create a variable indexed by the name of a product or a location. The syntax is the same as that above, except with an arbitrary vector as an index as opposed to a one-based range. The biggest difference is that instead of returning an Array of JuMP variables, JuMP will return a DenseAxisArray. For example:julia> @variable(model, x[1:2, [:A,:B]])\n2-dimensional DenseAxisArray{VariableRef,2,...} with index sets:\n    Dimension 1, 1:2\n    Dimension 2, Symbol[:A, :B]\nAnd data, a 2×2 Array{VariableRef,2}:\n x[1,A]  x[1,B]\n x[2,A]  x[2,B]DenseAxisArrays can be indexed and sliced as follows:julia> x[1, :A]\nx[1,A]\n\njulia> x[2, :]\n1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:\n    Dimension 1, Symbol[:A, :B]\nAnd data, a 2-element Array{VariableRef,1}:\n x[2,A]\n x[2,B]Similarly to the Array case, bounds can depend upon indices. For example:julia> @variable(model, x[i=2:3, j=1:2:3] >= 0.5i + j)\n2-dimensional DenseAxisArray{VariableRef,2,...} with index sets:\n    Dimension 1, 2:3\n    Dimension 2, 1:2:3\nAnd data, a 2×2 Array{VariableRef,2}:\n x[2,1]  x[2,3]\n x[3,1]  x[3,3]\n\njulia> lower_bound.(x)\n2-dimensional DenseAxisArray{Float64,2,...} with index sets:\n    Dimension 1, 2:3\n    Dimension 2, 1:2:3\nAnd data, a 2×2 Array{Float64,2}:\n 2.0  4.0\n 2.5  4.5"
},

{
    "location": "variables/#variable_sparseaxisarrays-1",
    "page": "Variables",
    "title": "SparseAxisArrays",
    "category": "section",
    "text": "The third container type that JuMP natively supports is SparseAxisArray. These arrays are created when the indices do not form a rectangular set. For example, this applies when indices have a dependence upon previous indices (called triangular indexing). JuMP supports this as follows:julia> @variable(model, x[i=1:2, j=i:2])\nJuMP.Containers.SparseAxisArray{VariableRef,2,Tuple{Any,Any}} with 3 entries:\n  [1, 2]  =  x[1,2]\n  [2, 2]  =  x[2,2]\n  [1, 1]  =  x[1,1]We can also conditionally create variables via a JuMP-specific syntax. This syntax appends a comparison check that depends upon the named indices and is separated from the indices by a semi-colon (;). For example:julia> @variable(model, x[i=1:4; mod(i, 2)==0])\nJuMP.Containers.SparseAxisArray{VariableRef,1,Tuple{Any}} with 2 entries:\n  [4]  =  x[4]\n  [2]  =  x[2]"
},

{
    "location": "variables/#variable_forcing-1",
    "page": "Variables",
    "title": "Forcing the container type",
    "category": "section",
    "text": "When creating a container of JuMP variables, JuMP will attempt to choose the tightest container type that can store the JuMP variables. Thus, it will prefer to create an Array before a DenseAxisArray and a DenseAxisArray before a SparseAxisArray. However, because this happens at compile time, it does not always make the best choice. To illustrate this, consider the following example:julia> A = 1:2\n1:2\n\njulia> @variable(model, x[A])\n1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:\n    Dimension 1, 1:2\nAnd data, a 2-element Array{VariableRef,1}:\n x[1]\n x[2]Since the value (and type) of A is unknown at parsing time, JuMP is unable to infer that A is a one-based integer range. Therefore, JuMP creates a DenseAxisArray, even though it could store these two variables in a standard one-dimensional Array.We can share our knowledge that it is possible to store these JuMP variables as an array by setting the container keyword:julia> @variable(model, y[A], container=Array)\n2-element Array{VariableRef,1}:\n y[1]\n y[2]JuMP now creates a vector of JuMP variables instead of a DenseAxisArray. Note that choosing an invalid container type will throw an error."
},

{
    "location": "variables/#Integrality-shortcuts-1",
    "page": "Variables",
    "title": "Integrality shortcuts",
    "category": "section",
    "text": "Adding integrality constraints to a model such as @constraint(model, x in MOI.ZeroOne()) and @constraint(model, x in MOI.Integer()) is a common operation. Therefore, JuMP supports two shortcuts for adding such constraints."
},

{
    "location": "variables/#Binary-(ZeroOne)-constraints-1",
    "page": "Variables",
    "title": "Binary (ZeroOne) constraints",
    "category": "section",
    "text": "Binary optimization variables are constrained to the set x in 0 1. (The MOI.ZeroOne set in MathOptInterface.) Binary optimization variables can be created in JuMP by passing Bin as an optional positional argument:julia> @variable(model, x, Bin)\nxWe can check if an optimization variable is binary by calling is_binary on the JuMP variable, and binary constraints can be removed with unset_binary.julia> is_binary(x)\ntrue\n\njulia> unset_binary(x)\n\njulia> is_binary(x)\nfalseBinary optimization variables can also be created by setting the binary keyword to true.julia> @variable(model, x, binary=true)\nx"
},

{
    "location": "variables/#Integer-constraints-1",
    "page": "Variables",
    "title": "Integer constraints",
    "category": "section",
    "text": "Integer optimization variables are constrained to the set x in mathbbZ. (The MOI.Integer set in MathOptInterface.) Integer optimization variables can be created in JuMP by passing Int as an optional positional argument:julia> @variable(model, x, Int)\nxInteger optimization variables can also be created by setting the integer keyword to true.julia> @variable(model, x, integer=true)\nxWe can check if an optimization variable is integer by calling is_integer on the JuMP variable, and integer constraints can be removed with unset_integer.julia> is_integer(x)\ntrue\n\njulia> unset_integer(x)\n\njulia> is_integer(x)\nfalse"
},

{
    "location": "variables/#Semidefinite-variables-1",
    "page": "Variables",
    "title": "Semidefinite variables",
    "category": "section",
    "text": "JuMP also supports modeling with semidefinite variables. A square symmetric matrix X is positive semidefinite if all eigenvalues are nonnegative. We can declare a matrix of JuMP variables to be positive semidefinite as follows:julia> @variable(model, x[1:2, 1:2], PSD)\n2×2 LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}:\n x[1,1]  x[1,2]\n x[1,2]  x[2,2]Note that x must be a square 2-dimensional Array of JuMP variables; it cannot be a DenseAxisArray or a SparseAxisArray. (See Variable containers, above, for more on this.)You can also impose a weaker constraint that the square matrix is only symmetric (instead of positive semidefinite) as follows:julia> @variable(model, x[1:2, 1:2], Symmetric)\n2×2 LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}:\n x[1,1]  x[1,2]\n x[1,2]  x[2,2]"
},

{
    "location": "variables/#Anonymous-JuMP-variables-1",
    "page": "Variables",
    "title": "Anonymous JuMP variables",
    "category": "section",
    "text": "In all of the above examples, we have created named JuMP variables. However, it is also possible to create so called anonymous JuMP variables. To create an anonymous JuMP variable, we drop the name of the variable from the macro call. This means dropping the second positional argument if the JuMP variable is a scalar, or dropping the name before the square bracket ([) if a container is being created. For example:julia> x = @variable(model)\nnonameThis shows how (model, x) is really short for:julia> x = model[:x] = @variable(model, base_name=\"x\")\nxAn Array of anonymous JuMP variables can be created as follows:julia> y = @variable(model, [i=1:2])\n2-element Array{VariableRef,1}:\n noname\n nonameIf necessary, you can store x in model as follows:julia> model[:x] = xThe <= and >= short-hand cannot be used to set bounds on anonymous JuMP variables. Instead, you should use the lower_bound and upper_bound keywords.Passing the Bin and Int variable types are also invalid. Instead, you should use the binary and integer keywords.Thus, the anonymous variant of @variable(model, x[i=1:2] >= i, Int) is:julia> x = @variable(model, [i=1:2], base_name=\"x\", lower_bound=i, integer=true)\n2-element Array{VariableRef,1}:\n x[1]\n x[2]warn: Warn\nCreating two named JuMP variables with the same name results in an error at runtime. Use anonymous variables as an alternative."
},

{
    "location": "variables/#User-defined-containers-1",
    "page": "Variables",
    "title": "User-defined containers",
    "category": "section",
    "text": "In the section Variable containers, we explained how JuMP supports the efficient creation of collections of JuMP variables in three types of containers. However, users are also free to create collections of JuMP variables in their own datastructures. For example, the following code creates a dictionary with symmetric matrices as the values:julia> variables = Dict{Symbol, Array{VariableRef,2}}()\nDict{Symbol,Array{VariableRef,2}} with 0 entries\n\njulia> for key in [:A, :B]\n           global variables[key] = @variable(model, [1:2, 1:2])\n       end\n\njulia> variables\nDict{Symbol,Array{VariableRef,2}} with 2 entries:\n  :A => VariableRef[noname noname; noname noname]\n  :B => VariableRef[noname noname; noname noname]"
},

{
    "location": "variables/#Deleting-variables-1",
    "page": "Variables",
    "title": "Deleting variables",
    "category": "section",
    "text": "JuMP supports the deletion of optimization variables.  To delete variables, we can use the delete method. We can also check whether x is a valid JuMP variable in model using the is_valid method:julia> @variable(model, x)\nx\n\njulia> is_valid(model, x)\ntrue\n\njulia> delete(model, x)\n\njulia> is_valid(model, x)\nfalse"
},

{
    "location": "variables/#Listing-all-variables-1",
    "page": "Variables",
    "title": "Listing all variables",
    "category": "section",
    "text": "Use JuMP.all_variables to obtain a list of all variables present in the model. This is useful for performing operations like:relaxing all integrality constraints in the model\nsetting the starting values for variables to the result of the last solve"
},

{
    "location": "variables/#JuMP.@variable",
    "page": "Variables",
    "title": "JuMP.@variable",
    "category": "macro",
    "text": "@variable(model, kw_args...)\n\nAdd an anonymous variable to the model model described by the keyword arguments kw_args and returns the variable.\n\n@variable(model, expr, args..., kw_args...)\n\nAdd a variable to the model model described by the expression expr, the positional arguments args and the keyword arguments kw_args. The expression expr can either be (note that in the following the symbol <= can be used instead of ≤ and the symbol >=can be used instead of ≥)\n\nof the form varexpr creating variables described by varexpr;\nof the form varexpr ≤ ub (resp. varexpr ≥ lb) creating variables described by varexpr with upper bounds given by ub (resp. lower bounds given by lb);\nof the form varexpr == value creating variables described by varexpr with fixed values given by value; or\nof the form lb ≤ varexpr ≤ ub or ub ≥ varexpr ≥ lb creating variables described by varexpr with lower bounds given by lb and upper bounds given by ub.\n\nThe expression varexpr can either be\n\nof the form varname creating a scalar real variable of name varname;\nof the form varname[...] or [...] creating a container of variables (see Containers in macros).\n\nThe recognized positional arguments in args are the following:\n\nBin: Sets the variable to be binary, i.e. either 0 or 1.\nInt: Sets the variable to be integer, i.e. one of ..., -2, -1, 0, 1, 2, ...\nSymmetric: Only available when creating a square matrix of variables, i.e. when varexpr is of the form varname[1:n,1:n] or varname[i=1:n,j=1:n]. It creates a symmetric matrix of variable, that is, it only creates a new variable for varname[i,j] with i ≤ j and sets varname[j,i] to the same variable as varname[i,j].\nPSD: The square matrix of variable is both Symmetric and constrained to be positive semidefinite.\n\nThe recognized keyword arguments in kw_args are the following:\n\nbase_name: Sets the name prefix used to generate variable names. It corresponds to the variable name for scalar variable, otherwise, the variable names are set to base_name[...] for each index ... of the axes axes.\nlower_bound: Sets the value of the variable lower bound.\nupper_bound: Sets the value of the variable upper bound.\nstart: Sets the variable starting value used as initial guess in optimization.\nbinary: Sets whether the variable is binary or not.\ninteger: Sets whether the variable is integer or not.\nvariable_type: See the \"Note for extending the variable macro\" section below.\ncontainer: Specify the container type, see Containers in macros.\n\nExamples\n\nThe following are equivalent ways of creating a variable x of name x with lower bound 0:\n\n# Specify everything in `expr`\n@variable(model, x >= 0)\n# Specify the lower bound using a keyword argument\n@variable(model, x, lower_bound=0)\n# Specify everything in `kw_args`\nx = @variable(model, base_name=\"x\", lower_bound=0)\n\nThe following are equivalent ways of creating a DenseAxisArray of index set [:a, :b] and with respective upper bounds 2 and 3 and names x[a] and x[b].\n\nub = Dict(:a => 2, :b => 3)\n# Specify everything in `expr`\n@variable(model, x[i=keys(ub)] <= ub[i])\n# Specify the upper bound using a keyword argument\n@variable(model, x[i=keys(ub)], upper_bound=ub[i])\n\nNote for extending the variable macro\n\nThe single scalar variable or each scalar variable of the container are created using add_variable(model, build_variable(_error, info, extra_args...; extra_kw_args...)) where\n\nmodel is the model passed to the @variable macro;\n_error is an error function with a single String argument showing the @variable call in addition to the error message given as argument;\ninfo is the VariableInfo struct containing the information gathered in expr, the recognized keyword arguments (except base_name and variable_type) and the recognized positional arguments (except Symmetric and PSD);\nextra_args are the unrecognized positional arguments of args plus the value of the variable_type keyword argument if present. The variable_type keyword argument allows the user to pass a position argument to build_variable without the need to give a positional argument to @variable. In particular, this allows the user to give a positional argument to the build_variable call when using the anonymous single variable syntax @variable(model, kw_args...); and\nextra_kw_args are the unrecognized keyword argument of kw_args.\n\nExamples\n\nThe following creates a variable x of name x with lower_bound 0 as with the first example above but does it without using the @variable macro\n\ninfo = VariableInfo(true, 0, false, NaN, false, NaN, false, NaN, false, false)\nJuMP.add_variable(model, JuMP.build_variable(error, info), \"x\")\n\nThe following creates a DenseAxisArray of index set [:a, :b] and with respective upper bounds 2 and 3 and names x[a] and x[b] as with the second example above but does it without using the @variable macro\n\n# Without the `@variable` macro\ndata = Vector{JuMP.variable_type(model)}(undef, length(keys(ub)))\nx = JuMP.Containers.DenseAxisArray(data, keys(ub))\nfor i in keys(ub)\n    info = VariableInfo(false, NaN, true, ub[i], false, NaN, false, NaN, false, false)\n    x[i] = JuMP.add_variable(model, JuMP.build_variable(error, info), \"x[$i]\")\nend\n\nThe following are equivalent ways of creating a Matrix of size N x N with variables custom variables created with a JuMP extension using the Poly(X) positional argument to specify its variables:\n\n# Using the `@variable` macro\n@variable(model, x[1:N,1:N], Symmetric, Poly(X))\n# Without the `@variable` macro\nx = Matrix{JuMP.variable_type(model, Poly(X))}(N, N)\ninfo = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)\nfor i in 1:N, j in i:N\n    x[i,j] = x[j,i] = JuMP.add_variable(model, build_variable(error, info, Poly(X)), \"x[$i,$j]\")\nend\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.owner_model",
    "page": "Variables",
    "title": "JuMP.owner_model",
    "category": "function",
    "text": "owner_model(s::AbstractJuMPScalar)\n\nReturn the model owning the scalar s.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.VariableRef",
    "page": "Variables",
    "title": "JuMP.VariableRef",
    "category": "type",
    "text": "VariableRef <: AbstractVariableRef\n\nHolds a reference to the model and the corresponding MOI.VariableIndex.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.all_variables",
    "page": "Variables",
    "title": "JuMP.all_variables",
    "category": "function",
    "text": "all_variables(model::Model)::Vector{VariableRef}\n\nReturns a list of all variables currently in the model. The variables are ordered by creation time.\n\nExample\n\nmodel = Model()\n@variable(model, x)\n@variable(model, y)\nall_variables(model)\n\n# output\n\n2-element Array{VariableRef,1}:\n x\n y\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.has_lower_bound",
    "page": "Variables",
    "title": "JuMP.has_lower_bound",
    "category": "function",
    "text": "has_lower_bound(v::VariableRef)\n\nReturn true if v has a lower bound. If true, the lower bound can be queried with lower_bound. See also LowerBoundRef.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.lower_bound",
    "page": "Variables",
    "title": "JuMP.lower_bound",
    "category": "function",
    "text": "lower_bound(v::VariableRef)\n\nReturn the lower bound of a variable. Error if one does not exist. See also has_lower_bound.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.set_lower_bound",
    "page": "Variables",
    "title": "JuMP.set_lower_bound",
    "category": "function",
    "text": "set_lower_bound(v::VariableRef, lower::Number)\n\nSet the lower bound of a variable. If one does not exist, create a new lower bound constraint. See also delete_lower_bound.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.delete_lower_bound",
    "page": "Variables",
    "title": "JuMP.delete_lower_bound",
    "category": "function",
    "text": "delete_lower_bound(v::VariableRef)\n\nDelete the lower bound constraint of a variable.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.has_upper_bound",
    "page": "Variables",
    "title": "JuMP.has_upper_bound",
    "category": "function",
    "text": "has_upper_bound(v::VariableRef)\n\nReturn true if v has a upper bound. If true, the upper bound can be queried with upper_bound. See also UpperBoundRef.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.upper_bound",
    "page": "Variables",
    "title": "JuMP.upper_bound",
    "category": "function",
    "text": "upper_bound(v::VariableRef)\n\nReturn the upper bound of a variable. Error if one does not exist. See also has_upper_bound.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.set_upper_bound",
    "page": "Variables",
    "title": "JuMP.set_upper_bound",
    "category": "function",
    "text": "set_upper_bound(v::VariableRef,upper::Number)\n\nSet the upper bound of a variable. If one does not exist, create an upper bound constraint. See also delete_upper_bound.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.delete_upper_bound",
    "page": "Variables",
    "title": "JuMP.delete_upper_bound",
    "category": "function",
    "text": "delete_upper_bound(v::VariableRef)\n\nDelete the upper bound constraint of a variable.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.is_fixed",
    "page": "Variables",
    "title": "JuMP.is_fixed",
    "category": "function",
    "text": "is_fixed(v::VariableRef)\n\nReturn true if v is a fixed variable. If true, the fixed value can be queried with fix_value. See also FixRef.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.fix_value",
    "page": "Variables",
    "title": "JuMP.fix_value",
    "category": "function",
    "text": "fix_value(v::VariableRef)\n\nReturn the value to which a variable is fixed. Error if one does not exist. See also is_fixed.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.fix",
    "page": "Variables",
    "title": "JuMP.fix",
    "category": "function",
    "text": "fix(v::VariableRef, value::Number; force::Bool = false)\n\nFix a variable to a value. Update the fixing constraint if one exists, otherwise create a new one. See also unfix.\n\nIf the variable already has variable bounds and force=false, calling fix will throw an error. If force=true, existing variable bounds will be deleted, and the fixing constraint will be added. Note a variable will have no bounds after a call to unfix.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.unfix",
    "page": "Variables",
    "title": "JuMP.unfix",
    "category": "function",
    "text": "unfix(v::VariableRef)\n\nDelete the fixing constraint of a variable.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.is_integer",
    "page": "Variables",
    "title": "JuMP.is_integer",
    "category": "function",
    "text": "is_integer(v::VariableRef)\n\nReturn true if v is constrained to be integer. See also IntegerRef.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.set_integer",
    "page": "Variables",
    "title": "JuMP.set_integer",
    "category": "function",
    "text": "set_integer(variable_ref::VariableRef)\n\nAdd an integrality constraint on the variable variable_ref. See also unset_integer.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.unset_integer",
    "page": "Variables",
    "title": "JuMP.unset_integer",
    "category": "function",
    "text": "unset_integer(variable_ref::VariableRef)\n\nRemove the integrality constraint on the variable variable_ref.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.IntegerRef",
    "page": "Variables",
    "title": "JuMP.IntegerRef",
    "category": "function",
    "text": "IntegerRef(v::VariableRef)\n\nReturn a constraint reference to the constraint constrainting v to be integer. Errors if one does not exist.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.is_binary",
    "page": "Variables",
    "title": "JuMP.is_binary",
    "category": "function",
    "text": "is_binary(v::VariableRef)\n\nReturn true if v is constrained to be binary. See also BinaryRef.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.set_binary",
    "page": "Variables",
    "title": "JuMP.set_binary",
    "category": "function",
    "text": "set_binary(v::VariableRef)\n\nAdd a constraint on the variable v that it must take values in the set 01. See also unset_binary.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.unset_binary",
    "page": "Variables",
    "title": "JuMP.unset_binary",
    "category": "function",
    "text": "unset_binary(variable_ref::VariableRef)\n\nRemove the binary constraint on the variable variable_ref.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.BinaryRef",
    "page": "Variables",
    "title": "JuMP.BinaryRef",
    "category": "function",
    "text": "BinaryRef(v::VariableRef)\n\nReturn a constraint reference to the constraint constrainting v to be binary. Errors if one does not exist.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.index-Tuple{VariableRef}",
    "page": "Variables",
    "title": "JuMP.index",
    "category": "method",
    "text": "index(v::VariableRef)::MOI.VariableIndex\n\nReturn the index of the variable that corresponds to v in the MOI backend.\n\n\n\n\n\n"
},

{
    "location": "variables/#JuMP.optimizer_index-Tuple{VariableRef}",
    "page": "Variables",
    "title": "JuMP.optimizer_index",
    "category": "method",
    "text": "optimizer_index(v::VariableRef)::MOI.VariableIndex\n\nReturn the index of the variable that corresponds to v in the optimizer model. It throws NoOptimizer if no optimizer is set and throws an ErrorException if the optimizer is set but is not attached.\n\n\n\n\n\n"
},

{
    "location": "variables/#Reference-1",
    "page": "Variables",
    "title": "Reference",
    "category": "section",
    "text": "@variable\nowner_model\nVariableRef\nall_variables\n\nhas_lower_bound\nlower_bound\nset_lower_bound\ndelete_lower_bound\n\nhas_upper_bound\nupper_bound\nset_upper_bound\ndelete_upper_bound\n\nis_fixed\nfix_value\nfix\nunfix\n\nis_integer\nset_integer\nunset_integer\nIntegerRef\n\nis_binary\nset_binary\nunset_binary\nBinaryRef\n\nindex(::VariableRef)\noptimizer_index(::VariableRef)"
},

{
    "location": "expressions/#",
    "page": "Expressions",
    "title": "Expressions",
    "category": "page",
    "text": "DocTestSetup = quote\n    using JuMP\nend"
},

{
    "location": "expressions/#Expressions-1",
    "page": "Expressions",
    "title": "Expressions",
    "category": "section",
    "text": "JuMP has three types of expressions: affine, quadratic, and nonlinear. These expressions can be inserted into constraints or into the objective. This is particularly useful if an expression is used in multiple places in the model."
},

{
    "location": "expressions/#Affine-expressions-1",
    "page": "Expressions",
    "title": "Affine expressions",
    "category": "section",
    "text": "There are four ways of constructing an affine expression in JuMP: with the @expression macro, with operator overloading, with the AffExpr constructor, and with add_to_expression!."
},

{
    "location": "expressions/#Macros-1",
    "page": "Expressions",
    "title": "Macros",
    "category": "section",
    "text": "The recommended way to create an affine expression is via the @expression macro.model = Model()\n@variable(model, x)\n@variable(model, y)\nex = @expression(model, 2x + y - 1)\n\n# output\n\n2 x + y - 1This expression can be used in the objective or added to a constraint. For example:@objective(model, Min, 2 * ex - 1)\nobjective_function(model)\n\n# output\n\n4 x + 2 y - 3Just like variables and constraints, named expressions can also be created. For examplemodel = Model()\n@variable(model, x[i = 1:3])\n@expression(model, expr[i = 1:3], i * sum(x[j] for j in i:3))\nexpr\n\n# output\n\n3-element Array{GenericAffExpr{Float64,VariableRef},1}:\n x[1] + x[2] + x[3]\n 2 x[2] + 2 x[3]\n 3 x[3]"
},

{
    "location": "expressions/#Operator-overloading-1",
    "page": "Expressions",
    "title": "Operator overloading",
    "category": "section",
    "text": "Expressions can also be created without macros. However, note that in some cases, this can be much slower that constructing an expression using macros.model = Model()\n@variable(model, x)\n@variable(model, y)\nex = 2x + y - 1\n\n# output\n\n2 x + y - 1"
},

{
    "location": "expressions/#Constructors-1",
    "page": "Expressions",
    "title": "Constructors",
    "category": "section",
    "text": "A third way to create an affine expression is by the AffExpr constructor. The first argument is the constant term, and the remaining arguments are variable-coefficient pairs.model = Model()\n@variable(model, x)\n@variable(model, y)\nex = AffExpr(-1.0, x => 2.0, y => 1.0)\n\n# output\n\n2 x + y - 1"
},

{
    "location": "expressions/#add_to_expression!-1",
    "page": "Expressions",
    "title": "add_to_expression!",
    "category": "section",
    "text": "The fourth way to create an affine expression is by using add_to_expression!. Compared to the operator overloading method, this approach is faster because it avoids constructing temporary objects. The @expression macro uses add_to_expression! behind-the-scenes.model = Model()\n@variable(model, x)\n@variable(model, y)\nex = AffExpr(-1.0)\nadd_to_expression!(ex, 2.0, x)\nadd_to_expression!(ex, 1.0, y)\n\n# output\n\n2 x + y - 1"
},

{
    "location": "expressions/#Removing-zero-terms-1",
    "page": "Expressions",
    "title": "Removing zero terms",
    "category": "section",
    "text": "Use drop_zeros! to remove terms from an affine expression with a 0 coefficient.julia> model = Model();\n\njulia> @variable(model, x)\nx\n\njulia> @expression(model, ex, 0 * x + 1)\n0 x + 1\n\njulia> drop_zeros!(ex)\n\njulia> ex\n1"
},

{
    "location": "expressions/#Quadratic-expressions-1",
    "page": "Expressions",
    "title": "Quadratic expressions",
    "category": "section",
    "text": "Like affine expressions, there are four ways of constructing a quadratic expression in JuMP: macros, operator overloading, constructors, and add_to_expression!."
},

{
    "location": "expressions/#Macros-2",
    "page": "Expressions",
    "title": "Macros",
    "category": "section",
    "text": "The @expression macro can be used to create quadratic expressions by including quadratic terms.model = Model()\n@variable(model, x)\n@variable(model, y)\nex = @expression(model, x^2 + 2 * x * y + y^2 + x + y - 1)\n\n# output\n\nx² + 2 x*y + y² + x + y - 1"
},

{
    "location": "expressions/#Operator-overloading-2",
    "page": "Expressions",
    "title": "Operator overloading",
    "category": "section",
    "text": "Operator overloading can also be used to create quadratic expressions. The same performance warning (discussed in the affine expression section) applies.model = Model()\n@variable(model, x)\n@variable(model, y)\nex = x^2 + 2 * x * y + y^2 + x + y - 1\n\n# output\n\nx² + 2 x*y + y² + x + y - 1"
},

{
    "location": "expressions/#Constructors-2",
    "page": "Expressions",
    "title": "Constructors",
    "category": "section",
    "text": "Quadratic expressions can also be created using the QuadExpr constructor. The first argument is an affine expression, and the remaining arguments are pairs, where the first term is a JuMP.UnorderedPair and the second term is the coefficient.model = Model()\n@variable(model, x)\n@variable(model, y)\naff_expr = AffExpr(-1.0, x => 1.0, y => 1.0)\nquad_expr = QuadExpr(aff_expr, UnorderedPair(x, x) => 1.0,\n                     UnorderedPair(x, y) => 2.0, UnorderedPair(y, y) => 1.0)\n\n# output\n\nx² + 2 x*y + y² + x + y - 1"
},

{
    "location": "expressions/#add_to_expression!-2",
    "page": "Expressions",
    "title": "add_to_expression!",
    "category": "section",
    "text": "Finally, add_to_expression! can also be used to add quadratic terms.model = Model()\n@variable(model, x)\n@variable(model, y)\nex = QuadExpr(x + y - 1.0)\nadd_to_expression!(ex, 1.0, x, x)\nadd_to_expression!(ex, 2.0, x, y)\nadd_to_expression!(ex, 1.0, y, y)\n\n# output\n\nx² + 2 x*y + y² + x + y - 1"
},

{
    "location": "expressions/#Removing-zero-terms-2",
    "page": "Expressions",
    "title": "Removing zero terms",
    "category": "section",
    "text": "Use drop_zeros! to remove terms from a quadratic expression with a 0 coefficient.julia> model = Model();\n\njulia> @variable(model, x)\nx\n\njulia> @expression(model, ex, 0 * x^2 + x + 1)\n0 x² + x + 1\n\njulia> drop_zeros!(ex)\n\njulia> ex\nx + 1"
},

{
    "location": "expressions/#Nonlinear-expressions-1",
    "page": "Expressions",
    "title": "Nonlinear expressions",
    "category": "section",
    "text": "Nonlinear expressions can be constructed only using the @NLexpression macro and can be used only in @NLobjective, @NLconstraint, and other @NLexpressions. Moreover, quadratic and affine expressions cannot be used in the nonlinear macros. For more details, see the Nonlinear Modeling section."
},

{
    "location": "expressions/#JuMP.@expression",
    "page": "Expressions",
    "title": "JuMP.@expression",
    "category": "macro",
    "text": "@expression(args...)\n\nEfficiently builds a linear or quadratic expression but does not add to model immediately. Instead, returns the expression which can then be inserted in other constraints. For example:\n\n@expression(m, shared, sum(i*x[i] for i=1:5))\n@constraint(m, shared + y >= 5)\n@constraint(m, shared + z <= 10)\n\nThe ref accepts index sets in the same way as @variable, and those indices can be used in the construction of the expressions:\n\n@expression(m, expr[i=1:3], i*sum(x[j] for j=1:3))\n\nAnonymous syntax is also supported:\n\nexpr = @expression(m, [i=1:3], i*sum(x[j] for j=1:3))\n\n\n\n\n\n"
},

{
    "location": "expressions/#JuMP.add_to_expression!",
    "page": "Expressions",
    "title": "JuMP.add_to_expression!",
    "category": "function",
    "text": "add_to_expression!(expression, terms...)\n\nUpdates expression in place to expression + (*)(terms...). This is typically much more efficient than expression += (*)(terms...). For example, add_to_expression!(expression, a, b) produces the same result as expression += a*b, and add_to_expression!(expression, a) produces the same result as expression += a.\n\nOnly a few methods are defined, mostly for internal use, and only for the cases when (1) they can be implemented efficiently and (2) expression is capable of storing the result. For example, add_to_expression!(::AffExpr, ::VariableRef, ::VariableRef) is not defined because a GenericAffExpr cannot store the product of two variables.\n\n\n\n\n\n"
},

{
    "location": "expressions/#JuMP.drop_zeros!",
    "page": "Expressions",
    "title": "JuMP.drop_zeros!",
    "category": "function",
    "text": "drop_zeros!(expr::GenericAffExpr)\n\nRemove terms in the affine expression with 0 coefficients.\n\n\n\n\n\ndrop_zeros!(expr::GenericQuadExpr)\n\nRemove terms in the quadratic expression with 0 coefficients.\n\n\n\n\n\n"
},

{
    "location": "expressions/#Reference-1",
    "page": "Expressions",
    "title": "Reference",
    "category": "section",
    "text": "@expression\nadd_to_expression!\ndrop_zeros!"
},

{
    "location": "objective/#",
    "page": "Objective",
    "title": "Objective",
    "category": "page",
    "text": ""
},

{
    "location": "objective/#Objective-1",
    "page": "Objective",
    "title": "Objective",
    "category": "section",
    "text": "This page describes macros and functions related to linear and quadratic objective functions only, unless otherwise indicated. For nonlinear objective functions, see Nonlinear Modeling.Use the @objective macro to set linear and quadratic objective functions in a JuMP model. The functions set_objective_sense and set_objective_function provide an equivalent lower-level interface.To query the objective function from a model, see objective_sense, objective_function, and objective_function_type.To query the optimal objective value or best known bound after a solve, see objective_value and objective_bound. These two functions apply to nonlinear objectives also."
},

{
    "location": "objective/#JuMP.@objective",
    "page": "Objective",
    "title": "JuMP.@objective",
    "category": "macro",
    "text": "@objective(model::Model, sense, func)\n\nSet the objective sense to sense and objective function to func. The objective sense can be either Min, Max, MathOptInterface.MIN_SENSE, MathOptInterface.MAX_SENSE or MathOptInterface.FEASIBILITY_SENSE; see MathOptInterface.ObjectiveSense. In order to set the sense programatically, i.e., when sense is a Julia variable whose value is the sense, one of the three MathOptInterface.ObjectiveSense values should be used. The function func can be a single JuMP variable, an affine expression of JuMP variables or a quadratic expression of JuMP variables.\n\nExamples\n\nTo minimize the value of the variable x, do as follows:\n\njulia> model = Model()\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\n\njulia> @variable(model, x)\nx\n\njulia> @objective(model, Min, x)\nx\n\nTo maximize the value of the affine expression 2x - 1, do as follows:\n\njulia> @objective(model, Max, 2x - 1)\n2 x - 1\n\nTo set a quadratic objective and set the objective sense programatically, do as follows:\n\njulia> sense = MOI.MIN_SENSE\nMIN_SENSE::OptimizationSense = 0\n\njulia> @objective(model, sense, x^2 - 2x + 1)\nx² - 2 x + 1\n\n\n\n\n\n"
},

{
    "location": "objective/#JuMP.set_objective_sense",
    "page": "Objective",
    "title": "JuMP.set_objective_sense",
    "category": "function",
    "text": "set_objective_sense(model::Model, sense::MathOptInterface.OptimizationSense)\n\nSets the objective sense of the model to the given sense. See set_objective_function to set the objective function. These are low-level functions; the recommended way to set the objective is with the @objective macro.\n\n\n\n\n\n"
},

{
    "location": "objective/#JuMP.set_objective_function",
    "page": "Objective",
    "title": "JuMP.set_objective_function",
    "category": "function",
    "text": "set_objective_function(\n    model::Model,\n    func::Union{AbstractJuMPScalar, MathOptInterface.AbstractScalarFunction})\n\nSets the objective function of the model to the given function. See set_objective_sense to set the objective sense. These are low-level functions; the recommended way to set the objective is with the @objective macro.\n\n\n\n\n\n"
},

{
    "location": "objective/#JuMP.objective_sense",
    "page": "Objective",
    "title": "JuMP.objective_sense",
    "category": "function",
    "text": "objective_sense(model::Model)::MathOptInterface.OptimizationSense\n\nReturn the objective sense.\n\n\n\n\n\n"
},

{
    "location": "objective/#JuMP.objective_function",
    "page": "Objective",
    "title": "JuMP.objective_function",
    "category": "function",
    "text": "objective_function(model::Model,\n               T::Type{<:AbstractJuMPScalar}=objective_function_type(model))\n\nReturn an object of type T representing the objective function. Error if the objective is not convertible to type T.\n\nExamples\n\njulia> model = Model()\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\n\njulia> @variable(model, x)\nx\n\njulia> @objective(model, Min, 2x + 1)\n2 x + 1\n\njulia> objective_function(model, AffExpr)\n2 x + 1\n\njulia> objective_function(model, QuadExpr)\n2 x + 1\n\njulia> typeof(objective_function(model, QuadExpr))\nGenericQuadExpr{Float64,VariableRef}\n\nWe see with the last two commands that even if the objective function is affine, as it is convertible to a quadratic function, it can be queried as a quadratic function and the result is quadratic.\n\nHowever, it is not convertible to a variable.\n\njulia> objective_function(model, VariableRef)\nERROR: InexactError: convert(MathOptInterface.SingleVariable, MathOptInterface.ScalarAffineFunction{Float64}(MathOptInterface.ScalarAffineTerm{Float64}[ScalarAffineTerm{Float64}(2.0, VariableIndex(1))], 1.0))\nStacktrace:\n [1] convert at /home/blegat/.julia/dev/MathOptInterface/src/functions.jl:398 [inlined]\n [2] get(::JuMP.JuMPMOIModel{Float64}, ::MathOptInterface.ObjectiveFunction{MathOptInterface.SingleVariable}) at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/model.jl:290\n [3] get at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/universalfallback.jl:114 [inlined]\n [4] get at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/cachingoptimizer.jl:439 [inlined]\n [5] get(::MathOptInterface.Bridges.LazyBridgeOptimizer{MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer,MathOptInterface.Utilities.UniversalFallback{JuMP.JuMPMOIModel{Float64}}},MathOptInterface.Bridges.AllBridgedConstraints{Float64}}, ::MathOptInterface.ObjectiveFunction{MathOptInterface.SingleVariable}) at /home/blegat/.julia/dev/MathOptInterface/src/Bridges/bridgeoptimizer.jl:172\n [6] objective_function(::Model, ::Type{VariableRef}) at /home/blegat/.julia/dev/JuMP/src/objective.jl:129\n [7] top-level scope at none:0\n\n\n\n\n\n"
},

{
    "location": "objective/#JuMP.objective_function_type",
    "page": "Objective",
    "title": "JuMP.objective_function_type",
    "category": "function",
    "text": "objective_function_type(model::Model)::AbstractJuMPScalar\n\nReturn the type of the objective function.\n\n\n\n\n\n"
},

{
    "location": "objective/#JuMP.objective_bound",
    "page": "Objective",
    "title": "JuMP.objective_bound",
    "category": "function",
    "text": "objective_bound(model::Model)\n\nReturn the best known bound on the optimal objective value after a call to optimize!(model).\n\n\n\n\n\n"
},

{
    "location": "objective/#JuMP.objective_value",
    "page": "Objective",
    "title": "JuMP.objective_value",
    "category": "function",
    "text": "objective_value(model::Model)\n\nReturn the objective value after a call to optimize!(model).\n\n\n\n\n\n"
},

{
    "location": "objective/#Reference-1",
    "page": "Objective",
    "title": "Reference",
    "category": "section",
    "text": "@objective\nJuMP.set_objective_sense\nJuMP.set_objective_function\n\nJuMP.objective_sense\nJuMP.objective_function\nJuMP.objective_function_type\n\nJuMP.objective_bound\nJuMP.objective_value"
},

{
    "location": "constraints/#",
    "page": "Constraints",
    "title": "Constraints",
    "category": "page",
    "text": "CurrentModule = JuMP\nDocTestSetup = quote\n    using JuMP\nend\nDocTestFilters = [r\"≤|<=\", r\"≥|>=\", r\" == | = \", r\" ∈ | in \"]"
},

{
    "location": "constraints/#Constraints-1",
    "page": "Constraints",
    "title": "Constraints",
    "category": "section",
    "text": "This page explains how to write various types of constraints in JuMP. Before reading further, please make sure you are familiar with JuMP models, and JuMP Variables. For nonlinear constraints, see Nonlinear Modeling instead.JuMP is based on the MathOptInterface (MOI) API. Because of this, JuMP thinks of a constraint as the restriction that the output of a function belongs to a set. For example, instead of representing a constraint a^top x le b as a less-than-or-equal-to constraint, JuMP models this as the scalar affine function a^top x belonging to the less-than set (-infty b. Thus, instead of a less-than-or-equal-to constraint, we consider this constraint to be a scalar affine -in- less than constraint. More generally, we use the shorthand function-in-set to refer to constraints composed of different types of functions and sets. In the rest of this page, we will introduce the different types of functions and sets that JuMP knows about as needed. You can read more details about this function-in-set concept in the MOI documentation.note: Note\nThe examples use MOI as an alias for the MathOptInterface module. This alias is defined by using JuMP. You may also define it in your code byimport MathOptInterface\nconst MOI = MathOptInterface"
},

{
    "location": "constraints/#The-@constraint-macro-1",
    "page": "Constraints",
    "title": "The @constraint macro",
    "category": "section",
    "text": "Constraints are added to a JuMP model using the @constraint macro. Here is an example of how to add the constraint 2x le 1 to a JuMP model:julia> @constraint(model, con, 2x <= 1)\ncon : 2 x <= 1.0Wasn\'t that easy! Let\'s unpack what happened, because just like @variable there are a few subtle things going on.The mathematical constraint 2x le 1 was added to the model.\nA Julia variable called con was created that is a reference to the constraint.\nThis Julia variable was stored in model and can be accessed by model[:con].\nJuMP set the name attribute (the one that is shown when printing) of the constraint to \"con\".Just like the Julia variables created in @variable, con can be bound to a different value. For example:julia> con\ncon : 2 x <= 1.0\n\njulia> con = 1\n1\n\njulia> con\n1However, the reference can be retrieved by querying the model using the symbolic name:julia> con = model[:con]\ncon : 2 x <= 1.0\n\njulia> con\ncon : 2 x <= 1.0Because the named variables and constraints are stored in the same namespace, creating a constraint with the same name as a variable or an existing constraint will result in an error. To overcome this limitation, it is possible to create anonymous constraints, just like it is possible to create Anonymous JuMP variables. This is done by dropping the second argument to @constraint:julia> con = @constraint(model, 2x <= 1)\n2 x <= 1.0It is also possible use different comparison operators (e.g., >= and ==) to create the following types of constraints:julia> @constraint(model, 2x >= 1)\n2 x >= 1.0\n\njulia> @constraint(model, 2x == 1)\n2 x = 1.0\n\njulia> @constraint(model, 1 <= 2x <= 3)\n2 x ∈ [1.0, 3.0]Note that JuMP normalizes the constraints by moving all of the terms containing variables to the left-hand side, and all of the constant terms to the right-hand side. Thus, we get:julia> @constraint(model, 2x + 1 <= 4x + 4)\n-2 x <= 3.0"
},

{
    "location": "constraints/#constraint_duality-1",
    "page": "Constraints",
    "title": "Duality",
    "category": "section",
    "text": "JuMP adopts the notion of conic duality from MOI. For linear programs, a feasible dual on a >= constraint is nonnegative and a feasible dual on a <= constraint is nonpositive. If the constraint is an equality constraint, it depends on which direction is binding.note: Note\nJuMP\'s definition of duality is independent of the objective sense. That is, the sign of feasible duals associated with a constraint depends on the direction of the constraint and not whether the problem is maximization or minimization. This is a different convention from linear programming duality in some common textbooks.The dual value associated with a constraint in the most recent solution can be accessed using the dual function. You can use the has_duals function to check whether the model has a dual solution available to query. For example:julia> model = Model();\n\njulia> @variable(model, x)\nx\n\njulia> @constraint(model, con, x <= 1)\ncon : x <= 1.0\n\njulia> has_duals(model)\nfalseDocTestSetup = quote\n    using JuMP\n    model = Model(\n        with_optimizer(\n            MOI.Utilities.MockOptimizer,\n            JuMP._MOIModel{Float64}(),\n            eval_objective_value = false,\n            eval_variable_constraint_dual = false));\n    @variable(model, x);\n    @constraint(model, con, x <= 1);\n    @objective(model, Max, -2x);\n    optimize!(model);\n    mock = backend(model).optimizer.model;\n    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)\n    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)\n    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(con), -2.0)\nendjulia> @objective(model, Min, -2x)\n-2 x\n\njulia> optimize!(model)\n\njulia> has_duals(model)\ntrue\n\njulia> dual(con)\n-2.0\n\njulia> @objective(model, Max, 2x)\n2 x\n\njulia> optimize!(model)\n\njulia> dual(con)\n-2.0To help users who may be less familiar with conic duality, JuMP provides the shadow_price function which returns a value that can be interpreted as the improvement in the objective in response to an infinitesimal relaxation (on the scale of one unit) in the right-hand side of the constraint. shadow_price can be used only on linear constraints with a <=, >=, or == comparison operator.In the example above, dual(con) returned -2.0 regardless of the optimization sense. However, in the second case when the optimization sense is Max, shadow_price returns:julia> shadow_price(con)\n2.0To query the dual variables associated with a variable bound, first obtain a constraint reference using one of UpperBoundRef, LowerBoundRef, or FixRef, and then call dual on the returned constraint reference. Note that in linear programming, the duals on variable bounds are also called the reduced costs (although the sign might differ from the one you expect).DocTestSetup = quote\n    using JuMP\nend"
},

{
    "location": "constraints/#JuMP.name-Tuple{ConstraintRef{Model,#s1,Shape} where Shape<:AbstractShape where #s1<:MathOptInterface.ConstraintIndex}",
    "page": "Constraints",
    "title": "JuMP.name",
    "category": "method",
    "text": "name(v::ConstraintRef)\n\nGet a constraint\'s name attribute.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.set_name-Tuple{ConstraintRef{Model,#s1,Shape} where Shape<:AbstractShape where #s1<:MathOptInterface.ConstraintIndex,String}",
    "page": "Constraints",
    "title": "JuMP.set_name",
    "category": "method",
    "text": "set_name(v::ConstraintRef, s::AbstractString)\n\nSet a constraint\'s name attribute.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.constraint_by_name",
    "page": "Constraints",
    "title": "JuMP.constraint_by_name",
    "category": "function",
    "text": "constraint_by_name(model::AbstractModel,\n                   name::String)::Union{ConstraintRef, Nothing}\n\nReturns the reference of the constraint with name attribute name or Nothing if no constraint has this name attribute. Throws an error if several constraints have name as their name attribute.\n\nconstraint_by_name(model::AbstractModel,\n                   name::String,\n                   F::Type{<:Union{AbstractJuMPScalar,\n                                   Vector{<:AbstractJuMPScalar},\n                                   MOI.AbstactFunction}},\n                   S::Type{<:MOI.AbstractSet})::Union{ConstraintRef, Nothing}\n\nSimilar to the method above, except that it throws an error if the constraint is not an F-in-S contraint where F is either the JuMP or MOI type of the function, and S is the MOI type of the set. This method is recommended if you know the type of the function and set since its returned type can be inferred while for the method above (i.e. without F and S), the exact return type of the constraint index cannot be inferred.\n\njulia> using JuMP\n\njulia> model = Model()\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\n\njulia> @variable(model, x)\nx\n\njulia> @constraint(model, con, x^2 == 1)\ncon : x² = 1.0\n\njulia> constraint_by_name(model, \"kon\")\n\njulia> constraint_by_name(model, \"con\")\ncon : x² = 1.0\n\njulia> constraint_by_name(model, \"con\", AffExpr, MOI.EqualTo{Float64})\n\njulia> constraint_by_name(model, \"con\", QuadExpr, MOI.EqualTo{Float64})\ncon : x² = 1.0\n\n\n\n\n\n"
},

{
    "location": "constraints/#Constraint-names-1",
    "page": "Constraints",
    "title": "Constraint names",
    "category": "section",
    "text": "The name, i.e. the value of the MOI.ConstraintName attribute, of a constraint can be obtained by name(::JuMP.ConstraintRef) and set by set_name(::JuMP.ConstraintRef, ::String).name(::JuMP.ConstraintRef{Model, <:JuMP.MOI.ConstraintIndex})\nset_name(::JuMP.ConstraintRef{Model, <:JuMP.MOI.ConstraintIndex}, ::String)The constraint can also be retrieved from its name using constraint_by_name.constraint_by_name"
},

{
    "location": "constraints/#Constraint-containers-1",
    "page": "Constraints",
    "title": "Constraint containers",
    "category": "section",
    "text": "So far, we\'ve added constraints one-by-one. However, just like Variable containers, JuMP provides a mechanism for building groups of constraints compactly. References to these groups of constraints are returned in containers. Three types of constraint containers are supported: Arrays, DenseAxisArrays, and SparseAxisArrays. We explain each of these in the following."
},

{
    "location": "constraints/#constraint_arrays-1",
    "page": "Constraints",
    "title": "Arrays",
    "category": "section",
    "text": "One way of adding a group of constraints compactly is the following:julia> @constraint(model, con[i = 1:3], i * x <= i + 1)\n3-element Array{ConstraintRef{Model,C,Shape} where Shape<:AbstractShape where C,1}:\n con[1] : x <= 2.0\n con[2] : 2 x <= 3.0\n con[3] : 3 x <= 4.0JuMP returns references to the three constraints in an Array that is bound to the Julia variable con. This array can be accessed and sliced as you would with any Julia array:julia> con[1]\ncon[1] : x <= 2.0\n\njulia> con[2:3]\n2-element Array{ConstraintRef{Model,C,Shape} where Shape<:AbstractShape where C,1}:\n con[2] : 2 x <= 3.0\n con[3] : 3 x <= 4.0Anonymous containers can also be constructed by dropping the name (e.g. con) before the square brackets:julia> @constraint(model, [i = 1:2], i * x <= i + 1)\n2-element Array{ConstraintRef{Model,C,Shape} where Shape<:AbstractShape where C,1}:\n x <= 2.0\n 2 x <= 3.0Just like @variable, JuMP will form an Array of constraints when it can determine at parse time that the indices are one-based integer ranges. Therefore con[1:b] will create an Array, but con[a:b] will not. A special case is con[Base.OneTo(n)] which will produce an Array. If JuMP cannot determine that the indices are one-based integer ranges (e.g., in the case of con[a:b]), JuMP will create a DenseAxisArray instead."
},

{
    "location": "constraints/#DenseAxisArrays-1",
    "page": "Constraints",
    "title": "DenseAxisArrays",
    "category": "section",
    "text": "The syntax for constructing a DenseAxisArray of constraints is very similar to the syntax for constructing a DenseAxisArray of variables.julia> @constraint(model, con[i = 1:2, j = 2:3], i * x <= j + 1)\n2-dimensional DenseAxisArray{ConstraintRef{Model,C,Shape} where Shape<:AbstractShape where C,2,...} with index sets:\n    Dimension 1, 1:2\n    Dimension 2, 2:3\nAnd data, a 2×2 Array{ConstraintRef{Model,C,Shape} where Shape<:AbstractShape where C,2}:\n con[1,2] : x <= 3.0    con[1,3] : x <= 4.0\n con[2,2] : 2 x <= 3.0  con[2,3] : 2 x <= 4.0"
},

{
    "location": "constraints/#SparseAxisArrays-1",
    "page": "Constraints",
    "title": "SparseAxisArrays",
    "category": "section",
    "text": "The syntax for constructing a SparseAxisArray of constraints is very similar to the syntax for constructing a SparseAxisArray of variables.julia> @constraint(model, con[i = 1:2, j = 1:2; i != j], i * x <= j + 1)\nJuMP.Containers.SparseAxisArray{ConstraintRef{Model,C,Shape} where Shape<:AbstractShape where C,2,Tuple{Any,Any}} with 2 entries:\n  [1, 2]  =  con[1,2] : x <= 3.0\n  [2, 1]  =  con[2,1] : 2 x <= 2.0"
},

{
    "location": "constraints/#Forcing-the-container-type-1",
    "page": "Constraints",
    "title": "Forcing the container type",
    "category": "section",
    "text": "When creating a container of constraints, JuMP will attempt to choose the tightest container type that can store the constraints. However, because this happens at parse time, it does not always make the best choice. Just like in @variable, we can force the type of container using the container keyword. For syntax and the reason behind this, take a look at the variable docs."
},

{
    "location": "constraints/#Vectorized-constraints-1",
    "page": "Constraints",
    "title": "Vectorized constraints",
    "category": "section",
    "text": "We can also add constraints to JuMP using vectorized linear algebra. For example:julia> @variable(model, x[i=1:2])\n2-element Array{VariableRef,1}:\n x[1]\n x[2]\n\njulia> A = [1 2; 3 4]\n2×2 Array{Int64,2}:\n 1  2\n 3  4\n\njulia> b = [5, 6]\n2-element Array{Int64,1}:\n 5\n 6\n\njulia> @constraint(model, con, A * x .== b)\n2-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.EqualTo{Float64}},ScalarShape},1}:\n x[1] + 2 x[2] == 5.0\n 3 x[1] + 4 x[2] == 6.0note: Note\nMake sure to use Julia\'s dot syntax in front of the comparison operators (e.g. .==, .>=, and .<=). If you use a comparison without the dot, an error will be thrown.Instead of adding an array of ScalarAffineFunction-in-EqualTo constraints, we can instead construct a VectorAffineFunction-in-Nonnegatives constraint as follows:julia> @constraint(model, A * x - b in MOI.Nonnegatives(2))\n[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] in MathOptInterface.Nonnegatives(2)In addition to the Nonnegatives set, MOI defines a number of other vector-valued sets such as Nonpositives. See the MOI documentation for more information.Note also that for the first time we have used an explicit function-in-set description of the constraint. Read more about this representation for constraints in the MOI documentation."
},

{
    "location": "constraints/#Constraints-on-a-single-variable-1",
    "page": "Constraints",
    "title": "Constraints on a single variable",
    "category": "section",
    "text": "In Variables, we saw how to modify the variable bounds, as well as add binary and integer restrictions to the domain of each variable. This can also be achieved using the @constraint macro. For example, MOI.ZeroOne() restricts the domain to 0 1:julia> @constraint(model, x in MOI.ZeroOne())\nx binaryand MOI.Integer() restricts to the domain to the integers mathbbZ:julia> @constraint(model, x in MOI.Integer())\nx integerJuMP also supports modeling semi-continuous variables, whose domain is 0  l u, using the MOI.Semicontinuous set:julia> @constraint(model, x in MOI.Semicontinuous(1.5, 3.5))\nx in MathOptInterface.Semicontinuous{Float64}(1.5, 3.5)as well as semi-integer variables, whose domain is 0  l l+1 dots u, using the MOI.Semiinteger set:julia> @constraint(model, x in MOI.Semiinteger(1.0, 3.0))\nx in MathOptInterface.Semiinteger{Float64}(1.0, 3.0)"
},

{
    "location": "constraints/#Quadratic-constraints-1",
    "page": "Constraints",
    "title": "Quadratic constraints",
    "category": "section",
    "text": "In addition to affine functions, JuMP also supports constraints with quadratic terms. (For more general nonlinear functions, see Nonlinear Modeling.) For example:julia> @variable(model, x[i=1:2])\n2-element Array{VariableRef,1}:\n x[1]\n x[2]\n\njulia> @variable(model, t >= 0)\nt\n\njulia> @constraint(model, x[1]^2 + x[2]^2 <= t^2)\nx[1]² + x[2]² - t² <= 0.0Note that this quadratic constraint (including the lower bound on t) is equivalent to a second order cone constraint where ||x[1]^2 + x[2]^2||\\_2 ≤ t and t ≥ 0. Instead of writing out the quadratic expansion, we can pass JuMP the constraint in function-in-set form. To do so, we need to define the function and the set.The function is a vector of variables:julia> [t, x[1], x[2]]\n3-element Array{VariableRef,1}:\n t\n x[1]\n x[2]Note that the variable t comes first, followed by the x arguments. The set is an instance of SecondOrderCone: SecondOrderCone(). Thus, we can add the second order cone constraint as follows:julia> @constraint(model, [t, x[1], x[2]] in SecondOrderCone())\n[t, x[1], x[2]] in MathOptInterface.SecondOrderCone(3)JuMP also supports the RotatedSecondOrderCone which requires the addition of a perspective variable u. The rotated second order cone constraints the variables t, u, and x such that: ||x[1]^2 + x[2]^2||\\_2 ≤ t × u and t, u ≥ 0. It can be added as follows:julia> @variable(model, u)\nu\n\njulia> @constraint(model, [t, u, x[1], x[2]] in RotatedSecondOrderCone())\n[t, u, x[1], x[2]] in MathOptInterface.RotatedSecondOrderCone(4)In addition to the second order cone and rotated second order cone, MOI defines a number of other conic sets such as the exponential and power cones. See the MathOptInterface documentation for more information."
},

{
    "location": "constraints/#Constraints-on-a-collection-of-variables-1",
    "page": "Constraints",
    "title": "Constraints on a collection of variables",
    "category": "section",
    "text": "In addition to constraining the domain of a single variable, JuMP supports placing constraints of a subset of the variables. We already saw an example of this in the Quadratic constraints section when we constrained a vector of variables to belong to the second order cone.In a special ordered set of type I (often denoted SOS-I), at most one variable can take a non-zero value. We can construct SOS-I constraints using the MOI.SOS1 set:julia> @variable(model, x[1:3])\n3-element Array{VariableRef,1}:\n x[1]\n x[2]\n x[3]\n\njulia> @constraint(model, x in MOI.SOS1([1.0, 2.0, 3.0]))\n[x[1], x[2], x[3]] in MathOptInterface.SOS1{Float64}([1.0, 2.0, 3.0])Note that we have to pass MOI.SOS1 a weight vector. This vector implies an ordering on the variables. If the decision variables are related and have a physical ordering (e.g., they correspond to the size of a factory to be built, and the SOS-I constraint enforces that only one factory can be built), then the weight vector, although not used directly in the constraint, can help the solver make a better decision in the solution process.This ordering is more important in a special ordered set of type II (SOS-II), in which at most two values can be non-zero, and if there are two non-zeros, they must be consecutive according to the ordering. For example, in the following constraint, the possible non-zero pairs are (x[1] and x[3]) and (x[2] and x[3]):julia> @constraint(model, x in MOI.SOS2([3.0, 1.0, 2.0]))\n[x[1], x[2], x[3]] in MathOptInterface.SOS2{Float64}([3.0, 1.0, 2.0])"
},

{
    "location": "constraints/#Semidefinite-constraints-1",
    "page": "Constraints",
    "title": "Semidefinite constraints",
    "category": "section",
    "text": "JuMP provides a special syntax for constraining a matrix to be symmetric positive semidefinite (PSD) with the @SDconstraint macro. In the context of this macro, the inequality A >= B between two square matrices A and B is understood as constraining A - B to be symmetric positive semidefinite.julia> @variable(model, x)\nx\n\njulia> @SDconstraint(model, [x 2x; 3x 4x] >= ones(2, 2))\n[x - 1    2 x - 1;\n 3 x - 1  4 x - 1] ∈ PSDCone()Solvers supporting such constraints usually expect to be given a matrix that is symbolically symmetric, that is, for which the expression in corresponding off-diagonal entries are the same. In our example, the expressions of entries (1, 2) and (2, 1) are respectively 2x - 1 and 3x - 1 which are different. To bridge the gap between the constraint modeled and what the solver expects, JuMP creates an equality constraint 3x - 1 == 2x - 1 and constrains the symmetric matrix [x - 1, 2 x - 1, 2 x - 1, 4 x - 1] to be positive semidefinite.note: Note\nIf the matrix provided is already symbolically symmetric, the equality constrains are equivalent to 0 = 0 and are not added. In practice, if all coefficients are smaller than 1e-10, the constraint is ignored, if all coefficients are smaller than 1e-8 but some are larger than 1e-10, it is ignored but a warning is displayed, otherwise if at least one coefficient is larger than 1e-8, the constraint is added.If the matrix is known to be symmetric, the PSD constraint can be added as follows:julia> using LinearAlgebra\n\njulia> @constraint(model, Symmetric([x 2x; 2x 4x] - ones(2, 2)) in PSDCone())\n[x - 1    2 x - 1;\n 2 x - 1  4 x - 1] ∈ PSDCone()Note that the lower triangular entries are silently ignored even if they are different so use it with caution:julia> cref = @constraint(model, Symmetric([x 2x; 3x 4x]) in PSDCone())\n[x    2 x;\n 2 x  4 x] ∈ PSDCone()\n\njulia> jump_function(constraint_object(cref))\n3-element Array{GenericAffExpr{Float64,VariableRef},1}:\n x\n 2 x\n 4 x\n\njulia> moi_set(constraint_object(cref))\nMathOptInterface.PositiveSemidefiniteConeTriangle(2)"
},

{
    "location": "constraints/#Constraint-modifications-1",
    "page": "Constraints",
    "title": "Constraint modifications",
    "category": "section",
    "text": "A common paradigm, especially in linear programming, is to repeatedly solve a model with different coefficients."
},

{
    "location": "constraints/#Modifying-a-constant-term-1",
    "page": "Constraints",
    "title": "Modifying a constant term",
    "category": "section",
    "text": "Most often, modifications involve changing the \"right-hand side\" of a linear constraint. This presents a challenge for JuMP because it leads to ambiguities. For example, what is the right-hand side term of @constraint(model, 2x + 1 <= x - 3)? This applies more generally to any constant term in a function appearing in the objective or a constraint.To avoid these ambiguities, JuMP includes the ability to fix variables to a value using the fix function. Fixing a variable sets its lower and upper bound to the same value. Thus, changes in a constant term can be simulated by adding a dummy variable and fixing it to different values. Here is an example:julia> @variable(model, const_term)\nconst_term\n\njulia> @constraint(model, con, 2x <= const_term)\ncon : 2 x - const_term <= 0.0\n\njulia> fix(const_term, 1.0)note: Note\nEven though const_term is fixed, it is still a decision variable. Thus, const_term * x is bilinear. Fixed variables are not replaced with constants when communicating the problem to a solver."
},

{
    "location": "constraints/#Modifying-a-variable-coefficient-1",
    "page": "Constraints",
    "title": "Modifying a variable coefficient",
    "category": "section",
    "text": "It is also possible to modify the scalar coefficients (but notably not yet the quadratic coefficients) using the set_coefficient function. Here is an example:julia> @constraint(model, con, 2x <= 1)\ncon : 2 x <= 1.0\n\njulia> set_coefficient(con, x, 3)\n\njulia> con\ncon : 3 x <= 1.0"
},

{
    "location": "constraints/#Constraint-deletion-1",
    "page": "Constraints",
    "title": "Constraint deletion",
    "category": "section",
    "text": "Constraints can be deleted from a model using delete. Just like variable references, it is possible to check if a constraint reference is valid using is_valid. Here is an example of deleting a constraint:julia> @constraint(model, con, 2x <= 1)\ncon : 2 x <= 1.0\n\njulia> is_valid(model, con)\ntrue\n\njulia> delete(model, con)\n\njulia> is_valid(model, con)\nfalse"
},

{
    "location": "constraints/#Accessing-constraints-from-a-model-1",
    "page": "Constraints",
    "title": "Accessing constraints from a model",
    "category": "section",
    "text": "You can query the types of constraints currently present in the model by calling list_of_constraint_types. Then, given a function and set type, use num_constraints to access the number of constraints of this type and all_constraints to access a list of their references. Then use constraint_object to get an instance of an AbstractConstraint object, either ScalarConstraint or VectorConstraint, that stores the constraint data.julia> model = Model();\n\njulia> @variable(model, x[i=1:2] >= i, Int);\n\njulia> @constraint(model, x[1] + x[2] <= 1);\n\njulia> list_of_constraint_types(model)\n3-element Array{Tuple{DataType,DataType},1}:\n (VariableRef, MathOptInterface.Integer)\n (VariableRef, MathOptInterface.GreaterThan{Float64})\n (GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64})\n\njulia> num_constraints(model, VariableRef, MOI.Integer)\n2\n\njulia> all_constraints(model, VariableRef, MOI.Integer)\n2-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.Integer},ScalarShape},1}:\n x[1] integer\n x[2] integer\n\njulia> num_constraints(model, VariableRef, MOI.GreaterThan{Float64})\n2\n\njulia> all_constraints(model, VariableRef, MOI.GreaterThan{Float64})\n2-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.GreaterThan{Float64}},ScalarShape},1}:\n x[1] ≥ 1.0\n x[2] ≥ 2.0\n\njulia> num_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})\n1\n\njulia> less_than_constraints = all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})\n1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},1}:\n x[1] + x[2] ≤ 1.0\n\njulia> con = constraint_object(less_than_constraints[1])\nScalarConstraint{GenericAffExpr{Float64,VariableRef},MathOptInterface.LessThan{Float64}}(x[1] + x[2], MathOptInterface.LessThan{Float64}(1.0))\n\njulia> con.func\nx[1] + x[2]\n\njulia> con.set\nMathOptInterface.LessThan{Float64}(1.0)"
},

{
    "location": "constraints/#JuMP.@constraint",
    "page": "Constraints",
    "title": "JuMP.@constraint",
    "category": "macro",
    "text": "@constraint(m::Model, expr)\n\nAdd a constraint described by the expression expr.\n\n@constraint(m::Model, ref[i=..., j=..., ...], expr)\n\nAdd a group of constraints described by the expression expr parametrized by i, j, ...\n\nThe expression expr can either be\n\nof the form func in set constraining the function func to belong to the set set which is either a MathOptInterface.AbstractSet or one of the JuMP shortcuts SecondOrderCone, RotatedSecondOrderCone and PSDCone, e.g. @constraint(model, [1, x-1, y-2] in SecondOrderCone()) constrains the norm of [x-1, y-2] be less than 1;\nof the form a sign b, where sign is one of ==, ≥, >=, ≤ and <= building the single constraint enforcing the comparison to hold for the expression a and b, e.g. @constraint(m, x^2 + y^2 == 1) constrains x and y to lie on the unit circle;\nof the form a ≤ b ≤ c or a ≥ b ≥ c (where ≤ and <= (resp. ≥ and >=) can be used interchangeably) constraining the paired the expression b to lie between a and c;\nof the forms @constraint(m, a .sign b) or @constraint(m, a .sign b .sign c) which broadcast the constraint creation to each element of the vectors.\n\nNote for extending the constraint macro\n\nEach constraint will be created using add_constraint(m, build_constraint(_error, func, set)) where\n\n_error is an error function showing the constraint call in addition to the error message given as argument,\nfunc is the expression that is constrained\nand set is the set in which it is constrained to belong.\n\nFor expr of the first type (i.e. @constraint(m, func in set)), func and set are passed unchanged to build_constraint but for the other types, they are determined from the expressions and signs. For instance, @constraint(m, x^2 + y^2 == 1) is transformed into add_constraint(m, build_constraint(_error, x^2 + y^2, MOI.EqualTo(1.0))).\n\nTo extend JuMP to accept new constraints of this form, it is necessary to add the corresponding methods to build_constraint. Note that this will likely mean that either func or set will be some custom type, rather than e.g. a Symbol, since we will likely want to dispatch on the type of the function or set appearing in the constraint.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.@SDconstraint",
    "page": "Constraints",
    "title": "JuMP.@SDconstraint",
    "category": "macro",
    "text": "@SDconstraint(model::Model, expr)\n\nAdd a semidefinite constraint described by the expression expr.\n\n@SDconstraint(model::Model, ref[i=..., j=..., ...], expr)\n\nAdd a group of semidefinite constraints described by the expression expr parametrized by i, j, ...\n\nThe expression expr needs to be of the form a sign b where sign is ⪰, ≥, >=, ⪯, ≤ or <= and a and b are square matrices. It constrains the matrix x = a - b (or x = b - a if the sign is ⪯, ≤ or <=) to be symmetric and positive semidefinite.\n\nBy default, we check numerical symmetry of the matrix x, and if symmetry is violated by some arbitrary amount, we add explicit equality constraints. You can use Symmetric(x) in PSDCone() with the @constraint macro to skip these checks if you know the matrix must be symmetric; see PSDCone for more information.\n\nExamples\n\nThe following constrains the matrix [x-1 2x-2; -3 x-4] to be symmetric and positive semidefinite, that is, it constrains 2x-2 to be equal to -3 and constrains all eigenvalues of the matrix to be nonnegative.\n\njulia> model = Model();\n\njulia> @variable(model, x)\nx\n\njulia> a = [x 2x\n            0  x];\n\njulia> b = [1 2\n            3 4];\n\njulia> cref = @SDconstraint(model, a ⪰ b)\n[x - 1  2 x - 2;\n -3     x - 4  ] ∈ PSDCone()\n\njulia> jump_function(constraint_object(cref))\n4-element Array{GenericAffExpr{Float64,VariableRef},1}:\n x - 1\n -3\n 2 x - 2\n x - 4\n\njulia> moi_set(constraint_object(cref))\nMathOptInterface.PositiveSemidefiniteConeSquare(2)\n\nIn the set PositiveSemidefiniteConeSquare(2) in the last output, Square means that the matrix is passed as a square matrix as the corresponding off-diagonal entries need to be constrained to be equal. A similar set PositiveSemidefiniteConeTriangle exists which only uses the upper triangular part of the matrix assuming that it is symmetric, see PSDCone to see how to use it.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.SecondOrderCone",
    "page": "Constraints",
    "title": "JuMP.SecondOrderCone",
    "category": "type",
    "text": "SecondOrderCone\n\nSecond order cone object that can be used to constrain the euclidean norm of a vector x to be less than or equal to a nonnegative scalar t. This is a shortcut for the MathOptInterface.SecondOrderCone.\n\nExamples\n\nThe following constrains (x-1 x-2)_2 le t and t ge 0:\n\njulia> model = Model();\n\njulia> @variable(model, x)\nx\n\njulia> @variable(model, t)\nt\n\njulia> @constraint(model, [t, x-1, x-2] in SecondOrderCone())\n[t, x - 1, x - 2] ∈ MathOptInterface.SecondOrderCone(3)\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.RotatedSecondOrderCone",
    "page": "Constraints",
    "title": "JuMP.RotatedSecondOrderCone",
    "category": "type",
    "text": "RotatedSecondOrderCone\n\nRotated second order cone object that can be used to constrain the square of the euclidean norm of a vector x to be less than or equal to 2tu where t and u are nonnegative scalars. This is a shortcut for the MathOptInterface.RotatedSecondOrderCone.\n\nExamples\n\nThe following constrains (x-1 x-2)_2 le 2tx and t x ge 0:\n\njulia> model = Model();\n\njulia> @variable(model, x)\nx\n\njulia> @variable(model, t)\nt\n\njulia> @constraint(model, [t, x, x-1, x-2] in RotatedSecondOrderCone())\n[t, x, x - 1, x - 2] ∈ MathOptInterface.RotatedSecondOrderCone(4)\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.PSDCone",
    "page": "Constraints",
    "title": "JuMP.PSDCone",
    "category": "type",
    "text": "PSDCone\n\nPositive semidefinite cone object that can be used to constrain a square matrix to be positive semidefinite in the @constraint macro. If the matrix has type Symmetric then the columns vectorization (the vector obtained by concatenating the columns) of its upper triangular part is constrained to belong to the MOI.PositiveSemidefiniteConeTriangle set, otherwise its column vectorization is constrained to belong to the MOI.PositiveSemidefiniteConeSquare set.\n\nExamples\n\nConsider the following example:\n\njulia> model = Model();\n\njulia> @variable(model, x)\nx\n\njulia> a = [ x 2x\n            2x  x];\n\njulia> b = [1 2\n            2 4];\n\njulia> cref = @SDconstraint(model, a ⪰ b)\n[x - 1    2 x - 2;\n 2 x - 2  x - 4  ] ∈ PSDCone()\n\njulia> jump_function(constraint_object(cref))\n4-element Array{GenericAffExpr{Float64,VariableRef},1}:\n x - 1\n 2 x - 2\n 2 x - 2\n x - 4\n\njulia> moi_set(constraint_object(cref))\nMathOptInterface.PositiveSemidefiniteConeSquare(2)\n\nWe see in the output of the last command that the matrix the vectorization of the matrix is constrained to belong to the PositiveSemidefiniteConeSquare.\n\njulia> using LinearAlgebra # For Symmetric\n\njulia> cref = @constraint(model, Symmetric(a - b) in PSDCone())\n[x - 1    2 x - 2;\n 2 x - 2  x - 4  ] ∈ PSDCone()\n\njulia> jump_function(constraint_object(cref))\n3-element Array{GenericAffExpr{Float64,VariableRef},1}:\n x - 1\n 2 x - 2\n x - 4\n\njulia> moi_set(constraint_object(cref))\nMathOptInterface.PositiveSemidefiniteConeTriangle(2)\n\nAs we see in the output of the last command, the vectorization of only the upper triangular part of the matrix is constrained to belong to the PositiveSemidefiniteConeSquare.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.shadow_price",
    "page": "Constraints",
    "title": "JuMP.shadow_price",
    "category": "function",
    "text": "shadow_price(constraint::ConstraintRef)\n\nThe change in the objective from an infinitesimal relaxation of the constraint. This value is computed from dual and can be queried only when has_duals is true and the objective sense is MIN_SENSE or MAX_SENSE (not FEASIBILITY_SENSE). For linear constraints, the shadow prices differ at most in sign from the dual value depending on the objective sense.\n\nNotes\n\nThe function simply translates signs from dual and does not validate the conditions needed to guarantee the sensitivity interpretation of the shadow price. The caller is responsible, e.g., for checking whether the solver converged to an optimal primal-dual pair or a proof of infeasibility.\nThe computation is based on the current objective sense of the model. If this has changed since the last solve, the results will be incorrect.\nRelaxation of equality constraints (and hence the shadow price) is defined based on which sense of the equality constraint is active.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.set_coefficient",
    "page": "Constraints",
    "title": "JuMP.set_coefficient",
    "category": "function",
    "text": "set_coefficient(constraint::ConstraintRef, variable::VariableRef, value)\n\nSet the coefficient of variable in the constraint constraint to value.\n\nNote that prior to this step, JuMP will aggregate multiple terms containing the same variable. For example, given a constraint 2x + 3x <= 2, set_coefficient(c, x, 4) will create the constraint 4x <= 2.\n\nmodel = Model()\n@variable(model, x)\n@constraint(model, con, 2x + 3x <= 2)\nset_coefficient(con, x, 4)\ncon\n\n# output\n\ncon : 4 x <= 2.0\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.is_valid",
    "page": "Constraints",
    "title": "JuMP.is_valid",
    "category": "function",
    "text": "is_valid(model::Model, constraint_ref::ConstraintRef{Model})\n\nReturn true if constraint_ref refers to a valid constraint in model.\n\n\n\n\n\nis_valid(model::Model, variable_ref::VariableRef)\n\nReturn true if variable refers to a valid variable in model.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.delete",
    "page": "Constraints",
    "title": "JuMP.delete",
    "category": "function",
    "text": "delete(model::Model, constraint_ref::ConstraintRef)\n\nDelete the constraint associated with constraint_ref from the model model.\n\n\n\n\n\ndelete(model::Model, variable_ref::VariableRef)\n\nDelete the variable associated with variable_ref from the model model.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.LowerBoundRef",
    "page": "Constraints",
    "title": "JuMP.LowerBoundRef",
    "category": "function",
    "text": "LowerBoundRef(v::VariableRef)\n\nReturn a constraint reference to the lower bound constraint of v. Errors if one does not exist.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.UpperBoundRef",
    "page": "Constraints",
    "title": "JuMP.UpperBoundRef",
    "category": "function",
    "text": "UpperBoundRef(v::VariableRef)\n\nReturn a constraint reference to the upper bound constraint of v. Errors if one does not exist.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.FixRef",
    "page": "Constraints",
    "title": "JuMP.FixRef",
    "category": "function",
    "text": "FixRef(v::VariableRef)\n\nReturn a constraint reference to the constraint fixing the value of v. Errors if one does not exist.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.ConstraintRef",
    "page": "Constraints",
    "title": "JuMP.ConstraintRef",
    "category": "type",
    "text": "ConstraintRef\n\nHolds a reference to the model and the corresponding MOI.ConstraintIndex.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.list_of_constraint_types",
    "page": "Constraints",
    "title": "JuMP.list_of_constraint_types",
    "category": "function",
    "text": "list_of_constraint_types(model::Model)\n\nReturn a list of tuples of the form (F, S) where F is a JuMP function type and S is an MOI set type such that all_constraints(model, F, S) returns a nonempty list.\n\nExample\n\njulia> model = Model();\n\njulia> @variable(model, x >= 0, Bin);\n\njulia> @constraint(model, 2x <= 1);\n\njulia> list_of_constraint_types(model)\n3-element Array{Tuple{DataType,DataType},1}:\n (VariableRef, MathOptInterface.ZeroOne)\n (VariableRef, MathOptInterface.GreaterThan{Float64})\n (GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64})\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.all_constraints",
    "page": "Constraints",
    "title": "JuMP.all_constraints",
    "category": "function",
    "text": "all_constraints(model::Model, function_type, set_type)::Vector{<:ConstraintRef}\n\nReturn a list of all constraints currently in the model where the function has type function_type and the set has type set_type. The constraints are ordered by creation time.\n\nSee also list_of_constraint_types and num_constraints.\n\nExample\n\njulia> model = Model();\n\njulia> @variable(model, x >= 0, Bin);\n\njulia> @constraint(model, 2x <= 1);\n\njulia> all_constraints(model, VariableRef, MOI.GreaterThan{Float64})\n1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.GreaterThan{Float64}},ScalarShape},1}:\n x ≥ 0.0\n\njulia> all_constraints(model, VariableRef, MOI.ZeroOne)\n1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.ZeroOne},ScalarShape},1}:\n x binary\n\njulia> all_constraints(model, AffExpr, MOI.LessThan{Float64})\n1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},1}:\n 2 x ≤ 1.0\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.num_constraints",
    "page": "Constraints",
    "title": "JuMP.num_constraints",
    "category": "function",
    "text": "num_constraints(model::Model, function_type, set_type)::Int64\n\nReturn the number of constraints currently in the model where the function has type function_type and the set has type set_type.\n\nSee also list_of_constraint_types and all_constraints.\n\nExample\n\njulia> model = Model();\n\njulia> @variable(model, x >= 0, Bin);\n\njulia> @variable(model, y);\n\njulia> @constraint(model, y in MOI.GreaterThan(1.0));\n\njulia> @constraint(model, y <= 1.0);\n\njulia> @constraint(model, 2x <= 1);\n\njulia> num_constraints(model, VariableRef, MOI.GreaterThan{Float64})\n2\n\njulia> num_constraints(model, VariableRef, MOI.ZeroOne)\n1\n\njulia> num_constraints(model, AffExpr, MOI.LessThan{Float64})\n2\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.constraint_object",
    "page": "Constraints",
    "title": "JuMP.constraint_object",
    "category": "function",
    "text": "constraint_object(ref::ConstraintRef)\n\nReturn the underlying constraint data for the constraint referenced by ref.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.AbstractConstraint",
    "page": "Constraints",
    "title": "JuMP.AbstractConstraint",
    "category": "type",
    "text": "abstract type AbstractConstraint\n\nAn abstract base type for all constraint types. AbstractConstraints store the function and set directly, unlike ConstraintRefs that are merely references to constraints stored in a model. AbstractConstraints do not need to be attached to a model.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.ScalarConstraint",
    "page": "Constraints",
    "title": "JuMP.ScalarConstraint",
    "category": "type",
    "text": "struct ScalarConstraint\n\nThe data for a scalar constraint. The func field containts a JuMP object representing the function and the set field contains the MOI set. See also the documentation on JuMP\'s representation of constraints for more background.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.VectorConstraint",
    "page": "Constraints",
    "title": "JuMP.VectorConstraint",
    "category": "type",
    "text": "struct VectorConstraint\n\nThe data for a vector constraint. The func field containts a JuMP object representing the function and the set field contains the MOI set. The shape field contains an AbstractShape matching the form in which the constraint was constructed (e.g., by using matrices or flat vectors). See also the documentation on JuMP\'s representation of constraints.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.index-Tuple{ConstraintRef}",
    "page": "Constraints",
    "title": "JuMP.index",
    "category": "method",
    "text": "index(cr::ConstraintRef)::MOI.ConstraintIndex\n\nReturn the index of the constraint that corresponds to cr in the MOI backend.\n\n\n\n\n\n"
},

{
    "location": "constraints/#JuMP.optimizer_index-Tuple{ConstraintRef{Model,C,Shape} where Shape<:AbstractShape where C}",
    "page": "Constraints",
    "title": "JuMP.optimizer_index",
    "category": "method",
    "text": "optimizer_index(cr::ConstraintRef{Model})::MOI.ConstraintIndex\n\nReturn the index of the constraint that corresponds to cr in the optimizer model. It throws NoOptimizer if no optimizer is set and throws an ErrorException if the optimizer is set but is not attached or if the constraint is bridged.\n\n\n\n\n\n"
},

{
    "location": "constraints/#Reference-1",
    "page": "Constraints",
    "title": "Reference",
    "category": "section",
    "text": "@constraint\n@SDconstraint\nSecondOrderCone\nRotatedSecondOrderCone\nPSDCone\nshadow_price\nset_coefficient\nis_valid\nJuMP.delete\nLowerBoundRef\nUpperBoundRef\nFixRef\nConstraintRef\nlist_of_constraint_types\nall_constraints\nnum_constraints\nconstraint_object\nAbstractConstraint\nScalarConstraint\nVectorConstraint\nindex(::ConstraintRef)\noptimizer_index(::ConstraintRef{Model})"
},

{
    "location": "constraints/#JuMP.@build_constraint",
    "page": "Constraints",
    "title": "JuMP.@build_constraint",
    "category": "macro",
    "text": "@build_constraint(constraint_expr)\n\nConstructs a ScalarConstraint or VectorConstraint using the same machinery as @constraint but without adding the constraint to a model.\n\nConstraints using broadcast operators like x .<= 1 are also supported and will create arrays of ScalarConstraint or VectorConstraint.\n\nExamples\n\nmodel = Model();\n@variable(model, x);\n@build_constraint(2x >= 1)\n\n# output\nScalarConstraint{GenericAffExpr{Float64,VariableRef},MathOptInterface.GreaterThan{Float64}}(2 x, MathOptInterface.GreaterThan{Float64}(1.0))\n\n\n\n\n\n"
},

{
    "location": "constraints/#Constructing-constraints-without-adding-them-to-the-model-1",
    "page": "Constraints",
    "title": "Constructing constraints without adding them to the model",
    "category": "section",
    "text": "For advanced use cases.@build_constraint"
},

{
    "location": "containers/#",
    "page": "Containers",
    "title": "Containers",
    "category": "page",
    "text": "DocTestSetup = quote\n    using JuMP\nend"
},

{
    "location": "containers/#JuMP.Containers.DenseAxisArray",
    "page": "Containers",
    "title": "JuMP.Containers.DenseAxisArray",
    "category": "type",
    "text": "DenseAxisArray(data::Array{T, N}, axes...) where {T, N}\n\nConstruct a JuMP array with the underlying data specified by the data array and the given axes. Exactly N axes must be provided, and their lengths must match size(data) in the corresponding dimensions.\n\nExample\n\njulia> array = JuMP.Containers.DenseAxisArray([1 2; 3 4], [:a, :b], 2:3)\n2-dimensional DenseAxisArray{Int64,2,...} with index sets:\n    Dimension 1, Symbol[:a, :b]\n    Dimension 2, 2:3\nAnd data, a 2×2 Array{Int64,2}:\n 1  2\n 3  4\n\njulia> array[:b, 3]\n4\n\n\n\n\n\nDenseAxisArray{T}(undef, axes...) where T\n\nConstruct an uninitialized DenseAxisArray with element-type T indexed over the given axes.\n\nExample\n\njulia> array = JuMP.Containers.DenseAxisArray{Float64}(undef, [:a, :b], 1:2);\n\njulia> fill!(array, 1.0)\n2-dimensional DenseAxisArray{Float64,2,...} with index sets:\n    Dimension 1, Symbol[:a, :b]\n    Dimension 2, 1:2\nAnd data, a 2×2 Array{Float64,2}:\n 1.0  1.0\n 1.0  1.0\n\njulia> array[:a, 2] = 5.0\n5.0\n\njulia> array[:a, 2]\n5.0\n\njulia> array\n2-dimensional DenseAxisArray{Float64,2,...} with index sets:\n    Dimension 1, Symbol[:a, :b]\n    Dimension 2, 1:2\nAnd data, a 2×2 Array{Float64,2}:\n 1.0  5.0\n 1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "containers/#JuMP.Containers.SparseAxisArray",
    "page": "Containers",
    "title": "JuMP.Containers.SparseAxisArray",
    "category": "type",
    "text": "struct SparseAxisArray{T,N,K<:NTuple{N, Any}} <: AbstractArray{T,N}\n    data::Dict{K,T}\nend\n\nN-dimensional array with elements of type T where only a subset of the entries are defined. The entries with indices idx = (i1, i2, ..., iN) in keys(data) has value data[idx]. Note that as opposed to SparseArrays.AbstractSparseArray, the missing entries are not assumed to be zero(T), they are simply not part of the array. This means that the result of map(f, sa::SparseAxisArray) or f.(sa::SparseAxisArray) has the same sparsity structure than sa even if f(zero(T)) is not zero.\n\nExample\n\njulia> dict = Dict((:a, 2) => 1.0, (:a, 3) => 2.0, (:b, 3) => 3.0)\nDict{Tuple{Symbol,Int64},Float64} with 3 entries:\n  (:b, 3) => 3.0\n  (:a, 2) => 1.0\n  (:a, 3) => 2.0\n\njulia> array = JuMP.Containers.SparseAxisArray(dict)\nJuMP.Containers.SparseAxisArray{Float64,2,Tuple{Symbol,Int64}} with 3 entries:\n  [b, 3]  =  3.0\n  [a, 2]  =  1.0\n  [a, 3]  =  2.0\n\njulia> array[:b, 3]\n3.0\n\n\n\n\n\n"
},

{
    "location": "containers/#Containers-1",
    "page": "Containers",
    "title": "Containers",
    "category": "section",
    "text": "JuMP provides a specialized container similar to AxisArrays that enables indexing with non-integer indices. Normally these are created automatically by JuMP\'s macros. The following constructors can be used to create them manually.JuMP.Containers.DenseAxisArray\nJuMP.Containers.SparseAxisArray"
},

{
    "location": "containers/#JuMP.Containers.generate_container",
    "page": "Containers",
    "title": "JuMP.Containers.generate_container",
    "category": "function",
    "text": "generate_container(T, indexvars, indexsets, requestedtype)\n\nReturn a tuple, the first element of which is code that generates a container for objects of type T given the index variables, index sets, and requestedtype. requestedtype may be one of :Array, :DenseAxisArray, :SparseAxisArray, or :Auto. Return error-producing code if requested type is incompatible. For the case of :Auto, the following rules are used to determine the appropriate container:\n\nIf all index sets are either explicit 1:B objects for any B or symbols which refer to objects of type Base.OneTo, then an Array is generated of the appropriate size. Types of symbols/expressions are not known at compile time, so we defer to type-safe functions to check the Base.OneTo condition.\nIf condition (1) does not hold, and the index sets are independent (the index variable for one set does not appear in the definition of another), then an DenseAxisArray is generated of the appropriate size.\nOtherwise, generate an empty SparseAxisArray{T,N,NTuple{N,Any}}.\n\nThe second element of the return tuple is a Bool, true if the container type automatically checks for duplicate terms in the index sets and false otherwise.\n\nExamples\n\ngenerate_container(VariableRef, [:i,:j], [:(1:N), :(1:T)], :Auto)\n# Returns code equivalent to:\n# :(Array{VariableRef}(length(1:N), length(1:T))\n\ngenerate_container(VariableRef, [:i,:j], [:(1:N), :(2:T)], :Auto)\n# Returns code equivalent to:\n# :(JuMP.Containers.DenseAxisArray(undef, 1:N, 2:T))\n\ngenerate_container(VariableRef, [:i,:j], [:(1:N), :(S)], :Auto)\n# Returns code that generates an Array if S is of type Base.OneTo,\n# otherwise an DenseAxisArray.\n\ngenerate_container(VariableRef, [:i,:j], [:(1:N), :(1:j)], :Auto)\n# Returns code equivalent to:\n# :(Containers.SparseAxisArray(Dict{NTuple{N,Any},VariableRef}()))\n\n\n\n\n\n"
},

{
    "location": "containers/#Containers-in-macros-1",
    "page": "Containers",
    "title": "Containers in macros",
    "category": "section",
    "text": "The generate_container function encodes the logic for how containers are constructed in JuMP\'s macros.JuMP.Containers.generate_containerIn the @variable (resp. @constraint) macro, containers of variables (resp. constraints) can be created with the following syntax:name[index_set_1, index_set_2, ..., index_set_n] creating an n-dimensional container of name name; or\n[index_set_1, index_set_2, ..., index_set_n] creating an anonymous n-dimensional container.Each expression index_set_i can either beof the form index_set specifying that the ith index set of the container is index_set; or\nof the form index_name=index_set specifying that the ith index set of the container is index_set and allowing values used in the macro expression and keyword arguments to be expressions depending on the index_name.The macro then creates the container using the JuMP.Containers.generate_container function with the following arguments:VariableRef for the @variable macro and ConstraintRef for the @constraint macro.\nThe index variables and arbitrary symbols for dimensions for which no variable index is specified.\nThe index sets specified.\nThe value of the keyword argument if given or :Auto."
},

{
    "location": "solvers/#",
    "page": "Solvers",
    "title": "Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/#Interacting-with-solvers-1",
    "page": "Solvers",
    "title": "Interacting with solvers",
    "category": "section",
    "text": "A JuMP model keeps a MathOptInterface (MOI) backend of type MOI.ModelLike that stores the optimization problem and acts as the optimization solver. We call it an MOI backend and not optimizer as it can also be a wrapper around an optimization file format such as MPS that writes the JuMP model in a file. From JuMP, the MOI backend can be accessed using the backend function. JuMP can be viewed as a lightweight, user-friendly layer on top of the MOI backend, in the sense that:JuMP does not maintain any copy of the model outside this MOI backend.\nJuMP variable (resp. constraint) references are simple structures containing both a reference to the JuMP model and the MOI index of the variable (resp. constraint).\nJuMP gives the constraints to the MOI backend in the form provided by the user without doing any automatic reformulation.\nvariables additions, constraints additions/modifications and objective modifications are directly applied to the MOI backend thus expecting the backend to support such modifications.While this allows JuMP to be a thin wrapper on top of the solver API, as mentioned in the last point above, this seems rather demanding on the solver. Indeed, while some solvers support incremental building of the model and modifications before and after solve, other solvers only support the model being copied at once before solve. Moreover, it seems to require all solvers to implement all possible reformulations independently which seems both very ambitious and might generate a lot of duplicated code.These apparent limitations are addressed at level of MOI in a manner that is completely transparent to JuMP. While the MOI API may seem very demanding, it allows MOI models to be a succession of lightweight MOI layers that fill the gap between JuMP requirements and the solver capabilities. The remainder of this section describes how JuMP interacts with the MOI backend.JuMP models can be created in three different modes: AUTOMATIC, MANUAL and DIRECT."
},

{
    "location": "solvers/#JuMP.with_optimizer",
    "page": "Solvers",
    "title": "JuMP.with_optimizer",
    "category": "function",
    "text": "with_optimizer(constructor, args...; kwargs...)\n\nReturn an OptimizerFactory that creates optimizers using the constructor constructor with positional arguments args and keyword arguments kwargs.\n\nExamples\n\nThe following returns an optimizer factory that creates IpoptOptimizers using the constructor call IpoptOptimizer(print_level=0):\n\nwith_optimizer(IpoptOptimizer, print_level=0)\n\n\n\n\n\n"
},

{
    "location": "solvers/#JuMP.NoOptimizer",
    "page": "Solvers",
    "title": "JuMP.NoOptimizer",
    "category": "type",
    "text": "struct NoOptimizer <: Exception end\n\nNo optimizer is set. The optimizer can be provided at the Model constructor or at the optimize! call with with_optimizer.\n\n\n\n\n\n"
},

{
    "location": "solvers/#JuMP.optimize!",
    "page": "Solvers",
    "title": "JuMP.optimize!",
    "category": "function",
    "text": "optimize!(model::Model,\n          optimizer_factory::Union{Nothing, OptimizerFactory}=nothing;\n          ignore_optimize_hook=(model.optimize_hook === nothing))\n\nOptimize the model. If optimizer_factory is not nothing, it first sets the optimizer to a new one created using the optimizer factory. The factory can be created using the with_optimizer function. If optimizer_factory is nothing and no optimizer was set to model before calling this function, a NoOptimizer error is thrown.\n\nExamples\n\nThe optimizer factory can either be given in the Model constructor as follows:\n\nmodel = Model(with_optimizer(GLPK.Optimizer))\n# ...fill model with variables, constraints and objectives...\n# Solve the model with GLPK\noptimize!(model)\n\nor in the optimize! call as follows:\n\nmodel = Model()\n# ...fill model with variables, constraints and objectives...\n# Solve the model with GLPK\noptimize!(model, with_optimizer(GLPK.Optimizer))\n\n\n\n\n\n"
},

{
    "location": "solvers/#JuMP.Model-Tuple{}",
    "page": "Solvers",
    "title": "JuMP.Model",
    "category": "method",
    "text": "Model(; caching_mode::MOIU.CachingOptimizerMode=MOIU.AUTOMATIC,\n        bridge_constraints::Bool=true)\n\nReturn a new JuMP model without any optimizer; the model is stored the model in a cache. The mode of the CachingOptimizer storing this cache is caching_mode. The optimizer can be set later in the optimize! call. If bridge_constraints is true, constraints that are not supported by the optimizer are automatically bridged to equivalent supported constraints when an appropriate transformation is defined in the MathOptInterface.Bridges module or is defined in another module and is explicitely added.\n\n\n\n\n\n"
},

{
    "location": "solvers/#JuMP.Model-Tuple{OptimizerFactory}",
    "page": "Solvers",
    "title": "JuMP.Model",
    "category": "method",
    "text": "Model(optimizer_factory::OptimizerFactory;\n      caching_mode::MOIU.CachingOptimizerMode=MOIU.AUTOMATIC,\n      bridge_constraints::Bool=true)\n\nReturn a new JuMP model using the optimizer factory optimizer_factory to create the optimizer. The optimizer factory can be created by the with_optimizer function.\n\nExamples\n\nThe following creates a model using the optimizer IpoptOptimizer(print_level=0):\n\nmodel = Model(with_optimizer(Ipopt.Optimizer, print_level=0))\n\n\n\n\n\n"
},

{
    "location": "solvers/#Automatic-and-Manual-modes-1",
    "page": "Solvers",
    "title": "Automatic and Manual modes",
    "category": "section",
    "text": "In AUTOMATIC and MANUAL modes, two MOI layers are automatically applied to the optimizer:CachingOptimizer: maintains a cache of the model so that when the optimizer does not support an incremental change to the model, the optimizer\'s internal model can be discarded and restored from the cache just before optimization. The CachingOptimizer has two different modes: AUTOMATIC and MANUAL corresponding to the two JuMP modes with the same names.\nLazyBridgeOptimizer (this can be disabled using the bridge_constraints keyword argument to Model constructor): when a constraint added is not supported by the optimizer, it attempts to transform the constraint into an equivalent form, possibly adding new variables and constraints that are supported by the optimizer. The applied transformations are selected among known recipes which are called bridges. A few default bridges are defined in MOI but new ones can be defined and added to the LazyBridgeOptimizer used by JuMP.See the MOI documentation for more details on these two MOI layers.To attach an optimizer to a JuMP model, JuMP needs to create a new empty optimizer instance. New optimizer instances can be obtained using an OptimizerFactory that can be created using the with_optimizer function:with_optimizerThe factory can be provided either at model construction time or at optimize! time:NoOptimizer\nJuMP.optimize!New JuMP models are created using the Model constructor:Model()\nModel(::JuMP.OptimizerFactory)# TODO: how to control the caching optimizer states"
},

{
    "location": "solvers/#JuMP.direct_model",
    "page": "Solvers",
    "title": "JuMP.direct_model",
    "category": "function",
    "text": "direct_model(backend::MOI.ModelLike)\n\nReturn a new JuMP model using backend to store the model and solve it. As opposed to the Model constructor, no cache of the model is stored outside of backend and no bridges are automatically applied to backend. The absence of cache reduces the memory footprint but it is important to bear in mind the following implications of creating models using this direct mode:\n\nWhen backend does not support an operation, such as modifying constraints or adding variables/constraints after solving, an error is thrown. For models created using the Model constructor, such situations can be dealt with by storing the modifications in a cache and loading them into the optimizer when optimize! is called.\nNo constraint bridging is supported by default.\nThe optimizer used cannot be changed the model is constructed.\nThe model created cannot be copied.\n\n\n\n\n\n"
},

{
    "location": "solvers/#JuMP.backend",
    "page": "Solvers",
    "title": "JuMP.backend",
    "category": "function",
    "text": "backend(model::Model)\n\nReturn the lower-level MathOptInterface model that sits underneath JuMP. This model depends on which operating mode JuMP is in (manual, automatic, or direct), and whether there are any bridges in the model.\n\nIf JuMP is in direct mode (i.e., the model was created using direct_model), the backend with be the optimizer passed to direct_model. If JuMP is in manual or automatic mode, the backend is a MOI.Utilities.CachingOptimizer.\n\nThis function should only be used by advanced users looking to access low-level MathOptInterface or solver-specific functionality.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Direct-mode-1",
    "page": "Solvers",
    "title": "Direct mode",
    "category": "section",
    "text": "JuMP models can be created in DIRECT mode using the JuMP.direct_model function.JuMP.direct_modelJuMP.backend"
},

{
    "location": "solutions/#",
    "page": "Query Solutions",
    "title": "Query Solutions",
    "category": "page",
    "text": "CurrentModule = JuMP\nDocTestSetup = quote\n    using JuMP\nend"
},

{
    "location": "solutions/#Querying-Solutions-1",
    "page": "Query Solutions",
    "title": "Querying Solutions",
    "category": "section",
    "text": "So far we have seen all the elements and constructs related to writing a JuMP optimization model. In this section we reach the point of what to do with a solved problem. Suppose your model is named model. Right after the call to optimize!(model), it\'s natural to ask JuMP questions about the finished optimization step. Typical questions include:Why has the optimization process stopped? Did it hit the time limit or run into numerical issues?\nDo I have a solution to my problem?\nIs it optimal?\nDo I have a dual solution?JuMP follows closely the concepts defined in MathOptInterface (MOI) to answer user questions about a finished call to optimize!(model). There are three main steps in querying a solution:First, we can query the termination_status which will tell us why the optimization stopped. This could be due to a number of reasons. For example, the solver found an optimal solution, the problem was proven to be infeasible, or a user-provided limit such as a time limit was encountered. For more information, see the Termination statuses section below.Second, we can query the primal_status and the dual_status, which will tell us what kind of results we have for our primal and dual solutions. This might be an optimal primal-dual pair, a primal solution without  a corresponding dual solution, or a certificate of primal or dual infeasibility. For more information, see the Solution statuses section below.Third, we can query value and dual to obtain the primal and dual values of the optimization variables and constraints (if there are values to be queried)."
},

{
    "location": "solutions/#MathOptInterface.TerminationStatusCode",
    "page": "Query Solutions",
    "title": "MathOptInterface.TerminationStatusCode",
    "category": "type",
    "text": "TerminationStatusCode\n\nAn Enum of possible values for the TerminationStatus attribute. This attribute is meant to explain the reason why the optimizer stopped executing in the most recent call to optimize!.\n\nIf no call has been made to optimize!, then the TerminationStatus is:\n\nOPTIMIZE_NOT_CALLED: The algorithm has not started.\n\nOK\n\nThese are generally OK statuses, i.e., the algorithm ran to completion normally.\n\nOPTIMAL: The algorithm found a globally optimal solution.\nINFEASIBLE: The algorithm concluded that no feasible solution exists.\nDUAL_INFEASIBLE: The algorithm concluded that no dual bound exists for the problem. If, additionally, a feasible (primal) solution is known to exist, this status typically implies that the problem is unbounded, with some technical exceptions.\nLOCALLY_SOLVED: The algorithm converged to a stationary point, local optimal solution, could not find directions for improvement, or otherwise completed its search without global guarantees.\nLOCALLY_INFEASIBLE: The algorithm converged to an infeasible point or otherwise completed its search without finding a feasible solution, without guarantees that no feasible solution exists.\nINFEASIBLE_OR_UNBOUNDED: The algorithm stopped because it decided that the problem is infeasible or unbounded; this occasionally happens during MIP presolve.\n\nSolved to relaxed tolerances\n\nALMOST_OPTIMAL: The algorithm found a globally optimal solution to relaxed tolerances.\nALMOST_INFEASIBLE: The algorithm concluded that no feasible solution exists within relaxed tolerances.\nALMOST_DUAL_INFEASIBLE: The algorithm concluded that no dual bound exists for the problem within relaxed tolerances.\nALMOST_LOCALLY_SOLVED: The algorithm converged to a stationary point, local optimal solution, or could not find directions for improvement within relaxed tolerances.\n\nLimits\n\nThe optimizer stopped because of some user-defined limit.\n\nITERATION_LIMIT: An iterative algorithm stopped after conducting the maximum number of iterations.\nTIME_LIMIT: The algorithm stopped after a user-specified computation time.\nNODE_LIMIT: A branch-and-bound algorithm stopped because it explored a maximum number of nodes in the branch-and-bound tree.\nSOLUTION_LIMIT: The algorithm stopped because it found the required number of solutions. This is often used in MIPs to get the solver to return the first feasible solution it encounters.\nMEMORY_LIMIT: The algorithm stopped because it ran out of memory.\nOBJECTIVE_LIMIT: The algorthm stopped because it found a solution better than a minimum limit set by the user.\nNORM_LIMIT: The algorithm stopped because the norm of an iterate became too large.\nOTHER_LIMIT: The algorithm stopped due to a limit not covered by one of the above.\n\nProblematic\n\nThis group of statuses means that something unexpected or problematic happened.\n\nSLOW_PROGRESS: The algorithm stopped because it was unable to continue making progress towards the solution.\nNUMERICAL_ERROR: The algorithm stopped because it encountered unrecoverable numerical error.\nINVALID_MODEL: The algorithm stopped because the model is invalid.\nINVALID_OPTION: The algorithm stopped because it was provided an invalid option.\nINTERRUPTED: The algorithm stopped because of an interrupt signal.\nOTHER_ERROR: The algorithm stopped because of an error not covered by one of the statuses defined above.\n\n\n\n\n\n"
},

{
    "location": "solutions/#Termination-statuses-1",
    "page": "Query Solutions",
    "title": "Termination statuses",
    "category": "section",
    "text": "The reason why the optimization of model was finished is given bytermination_status(model)This function will return a MOI.TerminationStatusCode enum.MOI.TerminationStatusCode"
},

{
    "location": "solutions/#MathOptInterface.ResultStatusCode",
    "page": "Query Solutions",
    "title": "MathOptInterface.ResultStatusCode",
    "category": "type",
    "text": "ResultStatusCode\n\nAn Enum of possible values for the PrimalStatus and DualStatus attributes. The values indicate how to interpret the result vector.\n\nNO_SOLUTION: the result vector is empty.\nFEASIBLE_POINT: the result vector is a feasible point.\nNEARLY_FEASIBLE_POINT: the result vector is feasible if some constraint tolerances are relaxed.\nINFEASIBLE_POINT: the result vector is an infeasible point.\nINFEASIBILITY_CERTIFICATE: the result vector is an infeasibility certificate. If the PrimalStatus is INFEASIBILITY_CERTIFICATE, then the primal result vector is a certificate of dual infeasibility. If the DualStatus is INFEASIBILITY_CERTIFICATE, then the dual result vector is a proof of primal infeasibility.\nNEARLY_INFEASIBILITY_CERTIFICATE: the result satisfies a relaxed criterion for a certificate of infeasibility.\nUNKNOWN_RESULT_STATUS: the result vector contains a solution with an unknown interpretation.\nOTHER_RESULT_STATUS: the result vector contains a solution with an interpretation not covered by one of the statuses defined above.\n\n\n\n\n\n"
},

{
    "location": "solutions/#Solution-statuses-1",
    "page": "Query Solutions",
    "title": "Solution statuses",
    "category": "section",
    "text": "These statuses indicate what kind of result is available to be queried with value and dual. It\'s possible that no result is available to be queried.We can obtain these statuses by calling primal_status for the primal status, and dual_status for the dual status. Both will return a MOI.ResultStatusCode enum.MOI.ResultStatusCodeCommon status situations are described in the MOI docs."
},

{
    "location": "solutions/#Obtaining-solutions-1",
    "page": "Query Solutions",
    "title": "Obtaining solutions",
    "category": "section",
    "text": "Provided the primal status is not MOI.NO_SOLUTION, the primal solution can be obtained by calling value. For the dual solution, the function is dual. Calling has_values for the primal status and has_duals for the dual solution is an equivalent way to check whether the status is MOI.NO_SOLUTION. It is important to note that if has_values or has_duals return false, calls to value and dual might throw an error or return arbitrary values.The container type (e.g., scalar, vector, or matrix) of the returned solution (primal or dual) depends on the type of the variable or constraint. See AbstractShape and dual_shape for details.To call value or dual on containers of VariableRef or ConstraintRef, use the broadcast syntax, e.g., value.(x).The objective value of a solved problem can be obtained via objective_value. The best known bound on the optimal objective value can be obtained via objective_bound.The following is a recommended workflow for solving a model and querying the solution: using JuMP\nmodel = Model()\n@variable(model, x[1:10] >= 0)\n# ... other constraints ...\noptimize!(model)\n\nif termination_status(model) == MOI.OPTIMAL\n    optimal_solution = value.(x)\n    optimal_objective = objective_value(model)\nelseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)\n    suboptimal_solution = value.(x)\n    suboptimal_objective = objective_value(model)\nelse\n    error(\"The model was not solved correctly.\")\nend# TODO: How to accurately measure the solve time."
},

{
    "location": "solutions/#JuMP.termination_status",
    "page": "Query Solutions",
    "title": "JuMP.termination_status",
    "category": "function",
    "text": "termination_status(model::Model)\n\nReturn the reason why the solver stopped (i.e., the MathOptInterface model attribute TerminationStatus).\n\n\n\n\n\n"
},

{
    "location": "solutions/#JuMP.primal_status",
    "page": "Query Solutions",
    "title": "JuMP.primal_status",
    "category": "function",
    "text": "primal_status(model::Model)\n\nReturn the status of the most recent primal solution of the solver (i.e., the MathOptInterface model attribute PrimalStatus).\n\n\n\n\n\n"
},

{
    "location": "solutions/#JuMP.has_values",
    "page": "Query Solutions",
    "title": "JuMP.has_values",
    "category": "function",
    "text": "has_values(model::Model)\n\nReturn true if the solver has a primal solution available to query, otherwise return false. See also value.\n\n\n\n\n\n"
},

{
    "location": "solutions/#JuMP.value",
    "page": "Query Solutions",
    "title": "JuMP.value",
    "category": "function",
    "text": "value(cref::ConstraintRef)\n\nGet the primal value of this constraint in the result returned by a solver. That is, if cref is the reference of a constraint func-in-set, it returns the value of func evaluated at the value of the variables (given by value(::VariableRef)). Use has_values to check if a result exists before asking for values.\n\nNote\n\nFor scalar contraints, the constant is moved to the set so it is not taken into account in the primal value of the constraint. For instance, the constraint @constraint(model, 2x + 3y + 1 == 5) is transformed into 2x + 3y-in-MOI.EqualTo(4) so the value returned by this function is the evaluation of 2x + 3y. ```\n\n\n\n\n\nvalue(v::VariableRef)\n\nGet the value of this variable in the result returned by a solver. Use has_values to check if a result exists before asking for values.\n\n\n\n\n\nvalue(ex::GenericAffExpr, var_value::Function)\n\nEvaluate ex using var_value(v) as the value for each variable v.\n\n\n\n\n\nvalue(v::GenericAffExpr)\n\nEvaluate an GenericAffExpr given the result returned by a solver. Replaces getvalue for most use cases.\n\n\n\n\n\nvalue(p::NonlinearParameter)\n\nReturn the current value stored in the nonlinear parameter p.\n\nExample\n\nmodel = Model()\n@NLparameter(model, p == 10)\nvalue(p)\n\n# output\n10.0\n\n\n\n\n\nvalue(ex::NonlinearExpression, var_value::Function)\n\nEvaluate ex using var_value(v) as the value for each variable v.\n\n\n\n\n\nvalue(ex::NonlinearExpression)\n\nEvaluate ex using value as the value for each variable v.\n\n\n\n\n\n"
},

{
    "location": "solutions/#JuMP.dual_status",
    "page": "Query Solutions",
    "title": "JuMP.dual_status",
    "category": "function",
    "text": "dual_status(model::Model)\n\nReturn the status of the most recent dual solution of the solver (i.e., the MathOptInterface model attribute DualStatus).\n\n\n\n\n\n"
},

{
    "location": "solutions/#JuMP.has_duals",
    "page": "Query Solutions",
    "title": "JuMP.has_duals",
    "category": "function",
    "text": "has_duals(model::Model)\n\nReturn true if the solver has a dual solution available to query, otherwise return false.\n\nSee also dual and shadow_price.\n\n\n\n\n\n"
},

{
    "location": "solutions/#JuMP.dual",
    "page": "Query Solutions",
    "title": "JuMP.dual",
    "category": "function",
    "text": "dual(cref::ConstraintRef)\n\nGet the dual value of this constraint in the result returned by a solver. Use has_dual to check if a result exists before asking for values. See also shadow_price.\n\n\n\n\n\n"
},

{
    "location": "solutions/#JuMP.OptimizeNotCalled",
    "page": "Query Solutions",
    "title": "JuMP.OptimizeNotCalled",
    "category": "type",
    "text": "struct OptimizeNotCalled <: Exception end\n\nA result attribute cannot be queried before optimize! is called.\n\n\n\n\n\n"
},

{
    "location": "solutions/#MathOptInterface.optimize!",
    "page": "Query Solutions",
    "title": "MathOptInterface.optimize!",
    "category": "function",
    "text": "optimize!(optimizer::AbstractOptimizer)\n\nStart the solution procedure.\n\n\n\n\n\n"
},

{
    "location": "solutions/#Reference-1",
    "page": "Query Solutions",
    "title": "Reference",
    "category": "section",
    "text": "JuMP.termination_status\nJuMP.primal_status\nJuMP.has_values\nJuMP.value\nJuMP.dual_status\nJuMP.has_duals\nJuMP.dual\nOptimizeNotCalled\nMOI.optimize!"
},

{
    "location": "nlp/#",
    "page": "Nonlinear Modeling",
    "title": "Nonlinear Modeling",
    "category": "page",
    "text": ""
},

{
    "location": "nlp/#Nonlinear-Modeling-1",
    "page": "Nonlinear Modeling",
    "title": "Nonlinear Modeling",
    "category": "section",
    "text": "CurrentModule = JuMP\nDocTestSetup = quote\n    using JuMP\nendJuMP has support for general smooth nonlinear (convex and nonconvex) optimization problems. JuMP is able to provide exact, sparse second-order derivatives to solvers. This information can improve solver accuracy and performance.Nonlinear objectives and constraints are specified by using the @NLobjective and @NLconstraint macros. The familiar sum() syntax is supported within these macros, as well as prod() which analogously represents the product of the terms within. Note that the @objective and @constraint macros (and corresponding functions) do not currently support nonlinear expressions. However, a model can contain a mix of linear, quadratic, and nonlinear contraints or objective functions. Starting points may be provided by using the start keyword argument to @variable.For example, we can solve the classical Rosenbrock problem (with a twist) as follows:using Ipopt\nmodel = Model(with_optimizer(Ipopt.Optimizer))\n@variable(model, x, start = 0.0)\n@variable(model, y, start = 0.0)\n\n@NLobjective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)\n\noptimize!(model)\nprintln(\"x = \", value(x), \" y = \", value(y))\n\n# adding a (linear) constraint\n@constraint(model, x + y == 10)\noptimize!(model)\nprintln(\"x = \", value(x), \" y = \", value(y))See the JuMP examples directory for more examples (which include mle.jl, rosenbrock.jl, and clnlbeam.jl).The NLP solver tests contain additional examples."
},

{
    "location": "nlp/#Syntax-notes-1",
    "page": "Nonlinear Modeling",
    "title": "Syntax notes",
    "category": "section",
    "text": "The syntax accepted in nonlinear expressions is more restricted than the syntax for linear and quadratic expressions. We note some important points below.With the exception of the splatting syntax discussed below, all expressions must be simple scalar operations. You cannot use dot, matrix-vector products, vector slices, etc. Translate vector operations into explicit sum() operations or use the AffExpr plus auxiliary variable trick described below.\nThere is no operator overloading provided to build up nonlinear expressions. For example, if x is a JuMP variable, the code 3x will return an AffExpr object that can be used inside of future expressions and linear constraints. However, the code sin(x) is an error. All nonlinear expressions must be inside of macros.\nUser-defined Functions may be used within nonlinear expressions only after they are registered. For example, the follow code results in an error because register() must be called first to register my_function.model = Model()\nmy_function(a, b) = exp(a) * b\n@variable(model, x)\n@variable(model, y)\n@NLobjective(model, Min, my_function(x, y))\n\n# output\n\nERROR: Unrecognized function \"my_function\" used in nonlinear expression.AffExpr and QuadExpr objects cannot currently be used inside nonlinear expressions. Instead, introduce auxiliary variables, e.g.:    my_expr = dot(c, x) + 3y # where x and y are variables\n    @variable(model, aux)\n    @constraint(model, aux == my_expr)\n    @NLobjective(model, Min, sin(aux))You can declare embeddable nonlinear expressions with @NLexpression. For example:    @NLexpression(model, my_expr[i = 1:n], sin(x[i]))\n    @NLconstraint(model, my_constr[i = 1:n], my_expr[i] <= 0.5)Anonymous syntax is supported in @NLexpression and @NLconstraint:    my_expr = @NLexpression(model, [i = 1:n], sin(x[i]))\n    my_constr = @NLconstraint(model, [i = 1:n], my_expr[i] <= 0.5)The splatting operator ... is recognized in a very restricted setting for expanding function arguments. The expression splatted can be only a symbol. More complex expressions are not recognized.julia> model = Model();\n\njulia> @variable(model, x[1:3]);\n\njulia> @NLconstraint(model, *(x...) <= 1.0)\nx[1] * x[2] * x[3] - 1.0 ≤ 0\n\njulia> @NLconstraint(model, *((x / 2)...) <= 0.0)\nERROR: LoadError: Unexpected expression in (*)(x / 2...). JuMP supports splatting only symbols. For example, x... is ok, but (x + 1)..., [x; y]... and g(f(y)...) are not."
},

{
    "location": "nlp/#JuMP.@NLparameter",
    "page": "Nonlinear Modeling",
    "title": "JuMP.@NLparameter",
    "category": "macro",
    "text": "@NLparameter(model, param == value)\n\nCreate and return a nonlinear parameter param attached to the model model with initial value set to value. Nonlinear parameters may be used only in nonlinear expressions.\n\nExample\n\nmodel = Model()\n@NLparameter(model, x == 10)\nvalue(x)\n\n# output\n10.0\n\n@NLparameter(model, param_collection[...] == value_expr)\n\nCreate and return a collection of nonlinear parameters param_collection attached to the model model with initial value set to value_expr (may depend on index sets). Uses the same syntax for specifying index sets as @variable.\n\nExample\n\nmodel = Model()\n@NLparameter(model, y[i = 1:10] == 2 * i)\nvalue(y[9])\n\n# output\n18.0\n\n\n\n\n\n"
},

{
    "location": "nlp/#JuMP.value-Tuple{NonlinearParameter}",
    "page": "Nonlinear Modeling",
    "title": "JuMP.value",
    "category": "method",
    "text": "value(p::NonlinearParameter)\n\nReturn the current value stored in the nonlinear parameter p.\n\nExample\n\nmodel = Model()\n@NLparameter(model, p == 10)\nvalue(p)\n\n# output\n10.0\n\n\n\n\n\n"
},

{
    "location": "nlp/#JuMP.set_value-Tuple{NonlinearParameter,Number}",
    "page": "Nonlinear Modeling",
    "title": "JuMP.set_value",
    "category": "method",
    "text": "set_value(p::NonlinearParameter, v::Number)\n\nStore the value v in the nonlinear parameter p.\n\nExample\n\nmodel = Model()\n@NLparameter(model, p == 0)\nset_value(p, 5)\nvalue(p)\n\n# output\n5.0\n\n\n\n\n\n"
},

{
    "location": "nlp/#Nonlinear-Parameters-1",
    "page": "Nonlinear Modeling",
    "title": "Nonlinear Parameters",
    "category": "section",
    "text": "For nonlinear models only, JuMP offers a syntax for explicit \"parameter\" objects which can be used to modify a model in-place just by updating the value of the parameter. Nonlinear parameters are declared by using the @NLparameter macro and may be indexed by arbitrary sets analogously to JuMP variables and expressions. The initial value of the parameter must be provided on the right-hand side of the == sign. There is no anonymous syntax for creating parameters.@NLparameterYou may use value and set_value to query or update the value of a parameter.value(::JuMP.NonlinearParameter)\nset_value(::JuMP.NonlinearParameter, ::Number)Nonlinear parameters can be used within nonlinear expressions only:@NLparameter(model, x == 10)\n@variable(model, z)\n@objective(model, Max, x * z)             # Error: x is a nonlinear parameter.\n@NLobjective(model, Max, x * z)           # Ok.\n@expression(model, my_expr, x * z^2)      # Error: x is a nonlinear parameter.\n@NLexpression(model, my_nl_expr, x * z^2) # Ok.Nonlinear parameters are useful when solving nonlinear models in a sequence:using Ipopt\nmodel = Model(with_optimizer(Ipopt.Optimizer))\n@variable(model, z)\n@NLparameter(model, x == 1.0)\n@NLobjective(model, Min, (z - x)^2)\noptimize!(model)\nvalue(z) # Equals 1.0.\n\n# Now, update the value of x to solve a different problem.\nset_value(x, 5.0)\noptimize!(model)\nvalue(z) # Equals 5.0Using nonlinear parameters can be faster than creating a new model from scratch with updated data because JuMP is able to avoid repeating a number of steps in processing the model before handing it off to the solver."
},

{
    "location": "nlp/#User-defined-Functions-1",
    "page": "Nonlinear Modeling",
    "title": "User-defined Functions",
    "category": "section",
    "text": "JuMP\'s library of recognized univariate functions is derived from the Calculus.jl package. If you encounter a standard special function not currently supported by JuMP, consider contributing to the list of derivative rules there. In addition to this built-in list of functions, it is possible to register custom (user-defined) nonlinear functions to use within nonlinear expressions. JuMP does not support black-box optimization, so all user-defined functions must provide derivatives in some form. Fortunately, JuMP supports automatic differentiation of user-defined functions, a feature to our knowledge not available in any comparable modeling systems.note: Note\nAutomatic differentiation is not finite differencing. JuMP\'s automatically computed derivatives are not subject to approximation error.JuMP uses ForwardDiff.jl to perform automatic differentiation; see the ForwardDiff.jl documentation for a description of how to write a function suitable for automatic differentiation. The general guideline is to write code that is generic with respect to the number type; don\'t assume that the input to the function is Float64. To register a user-defined function with derivatives computed by automatic differentiation, use the register method as in the following example:my_square(x) = x^2\nmy_f(x,y) = (x - 1)^2 + (y - 2)^2\n\nmodel = Model()\n\nregister(model, :my_f, 2, my_f, autodiff=true)\nregister(model, :my_square, 1, my_square, autodiff=true)\n\n@variable(model, x[1:2] >= 0.5)\n@NLobjective(model, Min, my_f(x[1], my_square(x[2])))The above code creates a JuMP model with the objective function (x[1] - 1)^2 + (x[2]^2 - 2)^2. The first argument to register is the model for which the functions are registered. The second argument is a Julia symbol object which serves as the name of the user-defined function in JuMP expressions; the JuMP name need not be the same as the name of the corresponding Julia method. The third argument specifies how many arguments the function takes. The fourth argument is the name of the Julia method which computes the function, and autodiff=true instructs JuMP to compute exact gradients automatically.Forward-mode automatic differentiation as implemented by ForwardDiff.jl has a computational cost that scales linearly with the number of input dimensions. As such, it is not the most efficient way to compute gradients of user-defined functions if the number of input arguments is large. In this case, users may want to provide their own routines for evaluating gradients. The more general syntax for register which accepts user-provided derivative evaluation routines is:JuMP.register(model::Model, s::Symbol, dimension::Integer, f::Function,\n              ∇f::Function, ∇²f::Function)The input differs between functions which take a single input argument and functions which take more than one. For univariate functions, the derivative evaluation routines should return a number which represents the first and second-order derivatives respectively. For multivariate functions, the derivative evaluation routines will be passed a gradient vector which they must explicitly fill. Second-order derivatives of multivariate functions are not currently supported; this argument should be omitted. The following example sets up the same optimization problem as before, but now we explicitly provide evaluation routines for the user-defined functions:my_square(x) = x^2\nmy_square_prime(x) = 2x\nmy_square_prime_prime(x) = 2\n\nmy_f(x, y) = (x - 1)^2 + (y - 2)^2\nfunction ∇f(g, x, y)\n    g[1] = 2 * (x - 1)\n    g[2] = 2 * (y - 2)\nend\n\nmodel = Model()\n\nregister(model, :my_f, 2, my_f, ∇f)\nregister(model, :my_square, 1, my_square, my_square_prime,\n         my_square_prime_prime)\n\n@variable(model, x[1:2] >= 0.5)\n@NLobjective(model, Min, my_f(x[1], my_square(x[2])))Once registered, user-defined functions can also be used in constraints. For example:@NLconstraint(model, my_square(x[1]) <= 2.0)"
},

{
    "location": "nlp/#User-defined-functions-with-vector-inputs-1",
    "page": "Nonlinear Modeling",
    "title": "User-defined functions with vector inputs",
    "category": "section",
    "text": "User-defined functions which take vectors as input arguments (e.g. f(x::Vector)) are not supported. Instead, use Julia\'s splatting syntax to create a function with scalar arguments. For example, instead off(x::Vector) = sum(x[i]^i for i in 1:length(x))define:f(x...) = sum(x[i]^i for i in 1:length(x))This function f can be used in a JuMP model as follows:model = Model()\n@variable(model, x[1:5] >= 0)\nf(x...) = sum(x[i]^i for i in 1:length(x))\nregister(model, :f, 5, f; autodiff = true)\n@NLobjective(model, Min, f(x...))"
},

{
    "location": "nlp/#Factors-affecting-solution-time-1",
    "page": "Nonlinear Modeling",
    "title": "Factors affecting solution time",
    "category": "section",
    "text": "The execution time when solving a nonlinear programming problem can be divided into two parts, the time spent in the optimization algorithm (the solver) and the time spent evaluating the nonlinear functions and corresponding derivatives. Ipopt explicitly displays these two timings in its output, for example:Total CPU secs in IPOPT (w/o function evaluations)   =      7.412\nTotal CPU secs in NLP function evaluations           =      2.083For Ipopt in particular, one can improve the performance by installing advanced sparse linear algebra packages, see Installation Guide. For other solvers, see their respective documentation for performance tips.The function evaluation time, on the other hand, is the responsibility of the modeling language. JuMP computes derivatives by using reverse-mode automatic differentiation with graph coloring methods for exploiting sparsity of the Hessian matrix [1]. As a conservative bound, JuMP\'s performance here currently may be expected to be within a factor of 5 of AMPL\'s."
},

{
    "location": "nlp/#Querying-derivatives-from-a-JuMP-model-1",
    "page": "Nonlinear Modeling",
    "title": "Querying derivatives from a JuMP model",
    "category": "section",
    "text": "For some advanced use cases, one may want to directly query the derivatives of a JuMP model instead of handing the problem off to a solver. Internally, JuMP implements the AbstractNLPEvaluator interface from MathOptInterface. To obtain an NLP evaluator object from a JuMP model, use JuMP.NLPEvaluator. JuMP.index returns the MOI.VariableIndex corresponding to a JuMP variable. MOI.VariableIndex itself is a type-safe wrapper for Int64 (stored in the value field.)For example:raw_index(v::MOI.VariableIndex) = v.value\nmodel = Model()\n@variable(model, x)\n@variable(model, y)\n@NLobjective(model, Min, sin(x) + sin(y))\nvalues = zeros(2)\nx_index = raw_index(JuMP.index(x))\ny_index = raw_index(JuMP.index(y))\nvalues[x_index] = 2.0\nvalues[y_index] = 3.0\nd = NLPEvaluator(model)\nMOI.initialize(d, [:Grad])\nMOI.eval_objective(d, values) # == sin(2.0) + sin(3.0)\n\n# output\n1.0504174348855488∇f = zeros(2)\nMOI.eval_objective_gradient(d, ∇f, values)\n(∇f[x_index], ∇f[y_index]) # == (cos(2.0), cos(3.0))\n\n# output\n(-0.4161468365471424, -0.9899924966004454)Only nonlinear constraints (those added with @NLconstraint), and nonlinear objectives (added with @NLobjective) exist in the scope of the NLPEvaluator. The NLPEvaluator does not evaluate derivatives of linear or quadratic constraints or objectives. The index method applied to a nonlinear constraint reference object returns its index as a NonlinearConstraintIndex. The value field of NonlinearConstraintIndex stores the raw integer index. For example:julia> model = Model();\n\njulia> @variable(model, x);\n\njulia> @NLconstraint(model, cons1, sin(x) <= 1);\n\njulia> @NLconstraint(model, cons2, x + 5 == 10);\n\njulia> typeof(cons1)\nConstraintRef{Model,NonlinearConstraintIndex,ScalarShape}\n\njulia> index(cons1)\nNonlinearConstraintIndex(1)\n\njulia> index(cons2)\nNonlinearConstraintIndex(2)# TODO: Provide a link for how to access the linear and quadratic parts of the\n# model.Note that for one-sided nonlinear constraints, JuMP subtracts any values on the right-hand side when computing expressions. In other words, one-sided nonlinear constraints are always transformed to have a right-hand side of zero.This method of querying derivatives directly from a JuMP model is convenient for interacting with the model in a structured way, e.g., for accessing derivatives of specific variables. For example, in statistical maximum likelihood estimation problems, one is often interested in the Hessian matrix at the optimal solution, which can be queried using the NLPEvaluator."
},

{
    "location": "nlp/#Raw-expression-input-1",
    "page": "Nonlinear Modeling",
    "title": "Raw expression input",
    "category": "section",
    "text": "In addition to the @NLobjective and @NLconstraint macros, it is also possible to provide Julia Expr objects directly by using set_NL_objective and add_NL_constraint. This input form may be useful if the expressions are generated programmatically. JuMP variables should be spliced into the expression object. For example:@variable(model, 1 <= x[i = 1:4] <= 5)\nset_NL_objective(model, :Min, :($(x[1])*$(x[4])*($(x[1])+$(x[2])+$(x[3])) + $(x[3])))\nadd_NL_constraint(model, :($(x[1])*$(x[2])*$(x[3])*$(x[4]) >= 25))\n\n# Equivalent form using traditional JuMP macros:\n@NLobjective(model, Min, x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3])\n@NLconstraint(model, x[1] * x[2] * x[3] * x[4] >= 25)See the Julia documentation for more examples and description of Julia expressions."
},

{
    "location": "nlp/#JuMP.@NLconstraint",
    "page": "Nonlinear Modeling",
    "title": "JuMP.@NLconstraint",
    "category": "macro",
    "text": "@NLconstraint(m::Model, expr)\n\nAdd a constraint described by the nonlinear expression expr. See also @constraint. For example:\n\n@NLconstraint(model, sin(x) <= 1)\n@NLconstraint(model, [i = 1:3], sin(i * x) <= 1 / i)\n\n\n\n\n\n"
},

{
    "location": "nlp/#JuMP.@NLexpression",
    "page": "Nonlinear Modeling",
    "title": "JuMP.@NLexpression",
    "category": "macro",
    "text": "@NLexpression(args...)\n\nEfficiently build a nonlinear expression which can then be inserted in other nonlinear constraints and the objective. See also [@expression]. For example:\n\n@NLexpression(model, my_expr, sin(x)^2 + cos(x^2))\n@NLconstraint(model, my_expr + y >= 5)\n@NLobjective(model, Min, my_expr)\n\nIndexing over sets and anonymous expressions are also supported:\n\n@NLexpression(m, my_expr_1[i=1:3], sin(i * x))\nmy_expr_2 = @NLexpression(m, log(1 + sum(exp(x[i])) for i in 1:2))\n\n\n\n\n\n"
},

{
    "location": "nlp/#JuMP.@NLobjective",
    "page": "Nonlinear Modeling",
    "title": "JuMP.@NLobjective",
    "category": "macro",
    "text": "@NLobjective(model, sense, expression)\n\nAdd a nonlinear objective to model with optimization sense sense. sense must be Max or Min.\n\nExample\n\n@NLobjective(model, Max, 2x + 1 + sin(x))\n\n\n\n\n\n"
},

{
    "location": "nlp/#Reference-1",
    "page": "Nonlinear Modeling",
    "title": "Reference",
    "category": "section",
    "text": "@NLconstraint\n@NLexpression\n@NLobjective[1]: Dunning, Huchette, and Lubin, \"JuMP: A Modeling Language for Mathematical Optimization\", SIAM Review, PDF."
},

{
    "location": "style/#",
    "page": "Style Guide",
    "title": "Style Guide",
    "category": "page",
    "text": ""
},

{
    "location": "style/#Style-guide-and-design-principles-1",
    "page": "Style Guide",
    "title": "Style guide and design principles",
    "category": "section",
    "text": ""
},

{
    "location": "style/#Style-guide-1",
    "page": "Style Guide",
    "title": "Style guide",
    "category": "section",
    "text": "This section describes the coding style rules that apply to JuMP code and that we recommend for JuMP models and surrounding Julia code. The motivations for a style guide include:conveying best practices for writing readable and maintainable code\nreducing the amount of time spent on bike-shedding by establishing basic naming and formatting conventions\nlowering the barrier for new contributors by codifying the existing practices (e.g., you can be more confident your code will pass review if you follow the style guide)In some cases, the JuMP style guide diverges from the Julia style guide. All such cases will be explicitly noted and justified.The JuMP style guide adopts many recommendations from the Google style guides.info: Info\nThe style guide is always a work in progress, and not all JuMP code follows the rules. When modifying JuMP, please fix the style violations of the surrounding code (i.e., leave the code tidier than when you started). If large changes are needed, consider separating them into another PR."
},

{
    "location": "style/#Formatting-1",
    "page": "Style Guide",
    "title": "Formatting",
    "category": "section",
    "text": "Julia unfortunately does not have an autoformatting tool like gofmt. Until a reliable autoformatting tool is available, we adopt the following conventions."
},

{
    "location": "style/#Whitespace-1",
    "page": "Style Guide",
    "title": "Whitespace",
    "category": "section",
    "text": "For conciseness, never use more than one blank line within a function, and never begin a function with a blank line.Bad:function foo(x)\n    y = 2 * x\n\n\n    return y\nend\n\nfunction foo(x)\n\n    y = 2 * x\n    return y\nendJulia is mostly insensitive to whitespace characters within lines. For consistency:Use spaces between binary operators (with some exceptions, see below)\nUse a single space after commas and semicolons\nDo not use extra spaces for unary operators, parentheses, or braces\nIndent within new blocks (except module) using 4 spacesGood:f(x, y) = [3 * dot(x, y); x\']Bad:f(x,y) = [ 3*dot(x,y) ; x\' ]Good:module Foo\n\nfunction f(x)\n    return x + 1\nend\n\nend # module Foo"
},

{
    "location": "style/#Exceptions-1",
    "page": "Style Guide",
    "title": "Exceptions",
    "category": "section",
    "text": "For aesthetic reasons, we make an exception for whitespace surrounding the exponential operator ^.Good:f(x) = x^2Bad:f(x) = x ^ 2We also make an exception for the : operator when it is used to form a range.Good:x = 1:5Bad:x = 1 : 5One reason is that it can be confused with Julia\'s conditional statement: cond ? x : y which requires whitespace around the :.We also make an exception for juxtaposed multiplication (i.e. dropping the * between a numeric literal and an expression) when the right-hand side is a symbol.Good:2x  # Acceptable if there are space constraints.\n2 * x  # This preferred if space is not an issue.\n2 * (x + 1)Bad:2(x + 1)"
},

{
    "location": "style/#Return-statements-1",
    "page": "Style Guide",
    "title": "Return statements",
    "category": "section",
    "text": "To avoid situations in which it is unclear whether the author intended to return a certain value or not, always use an explicit return statement to exit from a function. If the return from a function is nothing, use return instead of return nothing.We make an exception for assignment-form one-line functions (f(x) = 2x).Good:foo(x) = 2x  # Acceptable if one line\nfunction foo(x)\n    return 2x\nend\nfunction foo(x)\n    x[1] += 1\n    return\nendBad:function foo(x)\n    2x\nend\nfunction foo(x)\n    x[1] += 1\n    return nothing\nend"
},

{
    "location": "style/#TODO:-Line-breaks-1",
    "page": "Style Guide",
    "title": "TODO: Line breaks",
    "category": "section",
    "text": ""
},

{
    "location": "style/#Syntax-1",
    "page": "Style Guide",
    "title": "Syntax",
    "category": "section",
    "text": "Julia sometimes provides equivalent syntax to express the same basic operation. We discuss these cases below."
},

{
    "location": "style/#for-loops-1",
    "page": "Style Guide",
    "title": "for loops",
    "category": "section",
    "text": "Julia allows both for x = 1:N and for x in 1:N. Always prefer to use in over =, because in generalizes better to other index sets like for x in eachindex(A)."
},

{
    "location": "style/#Empty-vectors-1",
    "page": "Style Guide",
    "title": "Empty vectors",
    "category": "section",
    "text": "For a type T, T[] and Vector{T}() are equivalent ways to create an empty vector with element type T. Prefer T[] because it is more concise."
},

{
    "location": "style/#Trailing-periods-in-floating-point-constants-1",
    "page": "Style Guide",
    "title": "Trailing periods in floating-point constants",
    "category": "section",
    "text": "Both 1.0 and 1. create a Float64 with value 1.0. Prefer 1.0 over 1. because it is more easily distinguished from the integer constant 1.Moreover, as recommended by the Julia style guide, never use 1.0 when 1 is okay."
},

{
    "location": "style/#Comments-1",
    "page": "Style Guide",
    "title": "Comments",
    "category": "section",
    "text": "For non-native speakers and for general clarity, comments in code must be proper English sentences with appropriate punctuation.Good:# This is a comment demonstrating a good comment.Bad:# a bad comment"
},

{
    "location": "style/#JuMP-macro-syntax-1",
    "page": "Style Guide",
    "title": "JuMP macro syntax",
    "category": "section",
    "text": "For consistency, always use parentheses.Good:@variable(model, x >= 0)Bad:@variable model x >= 0For consistency, always use constant * variable as opposed to variable * constant. This makes it easier to read models in ambiguous cases like a * x.Good:a = 4\n@constraint(model, 3 * x <= 1)\n@constraint(model, a * x <= 1)Bad:a = 4\n@constraint(model, x * 3 <= 1)\n@constraint(model, x * a <= 1)In order to reduce boilerplate code, prefer the plural form of macros over lots of repeated calls to singular forms.Good:@variables(model, begin\n    x >= 0\n    y >= 1\n    z <= 2\nend)Bad:@variable(model, x >= 0)\n@variable(model, y >= 1)\n@variable(model, z <= 2)An exception is made for calls with many keyword arguments, since these need to be enclosed in parentheses in order to parse properly.Acceptable:@variable(model, x >= 0, start = 0.0, base_name = \"my_x\")\n@variable(model, y >= 1, start = 2.0)\n@variable(model, z <= 2, start = -1.0)Also acceptable:@variables(model, begin\n    x >= 0, (start = 0.0, base_name = \"my_x\")\n    y >= 1, (start = 2.0)\n    z <= 2, (start = -1.0)\nend)"
},

{
    "location": "style/#Naming-1",
    "page": "Style Guide",
    "title": "Naming",
    "category": "section",
    "text": "module SomeModule end\nfunction some_function end\nconst SOME_CONSTANT = ...\nstruct SomeStruct end\n@enum SomeEnum ENUM_VALUE_A ENUM_VALUE_B\nsome_local_variable = ...\nsome_file.jl # Except for ModuleName.jl."
},

{
    "location": "style/#Exported-and-non-exported-names-1",
    "page": "Style Guide",
    "title": "Exported and non-exported names",
    "category": "section",
    "text": "Begin private module level functions and constants with an underscore. All other objects in the scope of a module should be exported. (See JuMP.jl for an example of how to do this.)Names beginning with an underscore should only be used for distinguishing between exported (public) and non-exported (private) objects. Therefore, never begin the name of a local variable with an underscore.module MyModule\n\nexport public_function, PUBLIC_CONSTANT\n\nfunction _private_function()\n    local_variable = 1\n    return\nend\n\nfunction public_function end\n\nconst _PRIVATE_CONSTANT = 3.14159\nconst PUBLIC_CONSTANT = 1.41421\n\nend"
},

{
    "location": "style/#Use-of-underscores-within-names-1",
    "page": "Style Guide",
    "title": "Use of underscores within names",
    "category": "section",
    "text": "The Julia style guide recommends avoiding underscores \"when readable\", for example, haskey, isequal, remotecall, and remotecall_fetch. This convention creates the potential for unnecessary bikeshedding and also forces the user to recall the presence/absence of an underscore, e.g., \"was that argument named basename or base_name?\". For consistency, always use underscores in variable names and function names to separate words."
},

{
    "location": "style/#Use-of-!-1",
    "page": "Style Guide",
    "title": "Use of !",
    "category": "section",
    "text": "Julia has a convention of appending ! to a function name if the function modifies its arguments. We recommend to:Omit ! when the name itself makes it clear that modification is taking place, e.g., add_constraint and set_name. We depart from the Julia style guide because ! does not provide a reader with any additional information in this case, and adherence to this convention is not uniform even in base Julia itself (consider Base.println and Base.finalize).\nUse ! in all other cases. In particular it can be used to distinguish between modifying and non-modifying variants of the same function like scale and scale!.Note that ! is not a self-documenting feature because it is still ambiguous which arguments are modified when multiple arguments are present. Be sure to document which arguments are modified in the method\'s docstring.See also the Julia style guide recommendations for ordering of function arguments."
},

{
    "location": "style/#Abbreviations-1",
    "page": "Style Guide",
    "title": "Abbreviations",
    "category": "section",
    "text": "Abbreviate names to make the code more readable, not to save typing. Don\'t arbitrarily delete letters from a word to abbreviate it (e.g., indx). Use abbreviations consistently within a body of code (e.g., do not mix con and constr, idx and indx).Common abbreviations:num for numberTODO: add more"
},

{
    "location": "style/#Miscellaneous-1",
    "page": "Style Guide",
    "title": "Miscellaneous",
    "category": "section",
    "text": "(TODO: Rethink categories.)"
},

{
    "location": "style/#User-facing-MethodError-1",
    "page": "Style Guide",
    "title": "User-facing MethodError",
    "category": "section",
    "text": "Specifying argument types for methods is mostly optional in Julia, which means that it\'s possible to find out that you are working with unexpected types deep in the call chain. Avoid this situation or handle it with a helpful error message. A user should see a MethodError only for methods that they called directly.Bad:internal_function(x::Integer) = x + 1\n# The user sees a MethodError for internal_function when calling\n# public_function(\"a string\"). This is not very helpful.\npublic_function(x) = internal_function(x)Good:internal_function(x::Integer) = x + 1\n# The user sees a MethodError for public_function when calling\n# public_function(\"a string\"). This is easy to understand.\npublic_function(x::Integer) = internal_function(x)If it is hard to provide an error message at the top of the call chain, then the following pattern is also ok:internal_function(x::Integer) = x + 1\nfunction internal_function(x)\n    error(\"Internal error. This probably means that you called \" *\n          \"public_function() with the wrong type.\")\nend\npublic_function(x) = internal_function(x)"
},

{
    "location": "style/#@enum-vs.-Symbol-1",
    "page": "Style Guide",
    "title": "@enum vs. Symbol",
    "category": "section",
    "text": "The @enum macro lets you define types with a finite number of values that are explicitly enumerated (like enum in C/C++). Symbols are lightweight strings that are used to represent identifiers in Julia (for example, :x).@enum provides type safety and can have docstrings attached to explain the possible values. Use @enums when applicable, e.g., for reporting statuses. Use strings to provide long-form additional information like error messages.Use of Symbol should typically be reserved for identifiers, e.g., for lookup in the JuMP model (model[:my_variable])."
},

{
    "location": "style/#using-vs.-import-1",
    "page": "Style Guide",
    "title": "using vs. import",
    "category": "section",
    "text": "using ModuleName brings all symbols exported by the module ModuleName into scope, while import ModuleName brings only the module itself into scope. (See the Julia manual) for examples and more details.For the same reason that from <module> import * is not recommended in python (PEP 8), avoid using ModuleName except in throw-away scripts or at the REPL. The using statement makes it harder to track where symbols come from and exposes the code to ambiguities when two modules export the same symbol.Prefer using ModuleName: x, p to import ModuleName.x, ModuleName.p and import MyModule: x, p because the import versions allow method extension without qualifying with the module name."
},

{
    "location": "style/#Documentation-1",
    "page": "Style Guide",
    "title": "Documentation",
    "category": "section",
    "text": "This section describes the writing style that should be used when writing documentation for JuMP (and supporting packages).We can recommend the documentation style guides by Divio, Google, and Write the Docs as general reading for those writing documentation. This guide delegates a thorough handling of the topic to those guides and instead elaborates on the more points more specific to Julia and documentation that uses Documenter.Be concise.\nUse lists instead of long sentences.\nUse numbered lists when describing a sequence, e.g., (1) do X, (2) then Y.\nUse bullet points when the items are not ordered.\nExample code should be covered by doctests. (But it\'s unclear what to do if the code depends on a solver.)\nWhen a word is a Julia symbol and not an English word, enclose it with backticks. In addition, if it has a docstring in this doc add a link using @ref. If it is a plural, add the \"s\" after the closing backtick. For example,\n[`VariableRef`](@ref)s\nUse @meta blocks for TODOs and other comments that shouldn\'t be visible to readers. For example,\n```@meta\n# TODO: Mention also X, Y, and Z.\n```"
},

{
    "location": "style/#Design-principles-1",
    "page": "Style Guide",
    "title": "Design principles",
    "category": "section",
    "text": "TODO: How to structure and test large JuMP models, libraries that use JuMP.For how to write a solver, see MOI."
},

{
    "location": "extensions/#",
    "page": "Extensions",
    "title": "Extensions",
    "category": "page",
    "text": "CurrentModule = JuMP"
},

{
    "location": "extensions/#Extending-JuMP-1",
    "page": "Extensions",
    "title": "Extending JuMP",
    "category": "section",
    "text": "# TODO: How to extend JuMP: discussion on different ways to build on top of JuMP.\n# How to extend JuMP\'s macros and how to avoid doing this."
},

{
    "location": "extensions/#Extending-MOI-1",
    "page": "Extensions",
    "title": "Extending MOI",
    "category": "section",
    "text": "# TODO: Create new MOI function/sets, how to use it in JuMP"
},

{
    "location": "extensions/#JuMP.add_bridge",
    "page": "Extensions",
    "title": "JuMP.add_bridge",
    "category": "function",
    "text": " add_bridge(model::Model,\n            BridgeType::Type{<:MOI.Bridges.AbstractBridge})\n\nAdd BridgeType to the list of bridges that can be used to transform unsupported constraints into an equivalent formulation using only constraints supported by the optimizer.\n\n\n\n\n\n"
},

{
    "location": "extensions/#JuMP.BridgeableConstraint",
    "page": "Extensions",
    "title": "JuMP.BridgeableConstraint",
    "category": "type",
    "text": "struct BridgeableConstraint{C, B} <: AbstractConstraint\n    constraint::C\n    bridge_type::B\nend\n\nConstraint constraint that can be bridged by the bridge of type bridge_type. Adding this constraint to a model is equivalent to\n\nadd_bridge(model, bridge_type)\nadd_constraint(model, constraint)\n\nExamples\n\nGiven a new scalar set type CustomSet with a bridge CustomBridge that can bridge F-in-CustomSet constraints, when the user does\n\nmodel = Model()\n@variable(model, x)\n@constraint(model, x + 1 in CustomSet())\noptimize!(model)\n\nwith an optimizer that does not support F-in-CustomSet constraints, the constraint will not be bridge unless he manually calls add_bridge(model, CustomBridge). In order to automatically add the CustomBridge to any model to which an F-in-CustomSet is added, simply add the following method:\n\nfunction JuMP.build_constraint(_error::Function, func::AbstractJuMPScalar,\n                               set::CustomSet)\n    constraint = ScalarConstraint(func, set)\n    return JuMP.BridgeableConstraint(constraint, CustomBridge)\nend\n\nNote\n\nJuMP extensions should extend JuMP.build_constraint only if they also defined CustomSet, for three reasons:\n\nIt is problematic if multiple extensions overload the same JuMP method.\nA missing method will not inform the users that they forgot to load the extension module defining the build_constraint method.\nDefining a method where neither the function nor any of the argument types are defined in the package is called type piracy and is discouraged in the Julia style guide.\n\n```\n\n\n\n\n\n"
},

{
    "location": "extensions/#Adding-a-bridge-1",
    "page": "Extensions",
    "title": "Adding a bridge",
    "category": "section",
    "text": "# TODO: create new bridgeSee the bridge section in the MOI manual.add_bridge\nBridgeableConstraint"
},

{
    "location": "extensions/#Extending-JuMP-macros-1",
    "page": "Extensions",
    "title": "Extending JuMP macros",
    "category": "section",
    "text": "In order to provide a convenient syntax for the user to create variables, constraints or set the objective of a JuMP extension, it might be required to use macros similar to @variable, @constraint and @objective. It is recommended to first check whether it is possible to extend one of these three macros before creating a new one so as to leverage all their features and provide a more consistent interface to the user.### Extending the `@variable` macro\n\n# TODO: parse/build/add"
},

{
    "location": "extensions/#Extending-the-@constraint-macro-1",
    "page": "Extensions",
    "title": "Extending the @constraint macro",
    "category": "section",
    "text": "The @constraint macro always calls the same three functions:parse_constraint: is called at parsing time, it parses the constraint expression and returns a build_constraint call expression;\nbuild_constraint: given the functions and sets involved in the constraints, it returns a AbstractConstraint;\nadd_constraint: given the model, the AbstractConstraint constructed in build_constraint and the constraint name, it stores them in the model and returns a ConstraintRef.Adding methods to these functions is the recommended way to extend the @constraint macro."
},

{
    "location": "extensions/#JuMP.sense_to_set",
    "page": "Extensions",
    "title": "JuMP.sense_to_set",
    "category": "function",
    "text": "sense_to_set(_error::Function, ::Val{sense_symbol})\n\nConverts a sense symbol to a set set such that @constraint(model, func sense_symbol 0) is equivalent to@constraint(model, func in set)for anyfunc::AbstractJuMPScalar`.\n\nExample\n\nOnce a custom set is defined you can directly create a JuMP constraint with it:\n\njulia> struct CustomSet{T} <: MOI.AbstractScalarSet\n           value::T\n       end\n\njulia> model = Model();\n\njulia> @variable(model, x)\nx\n\njulia> cref = @constraint(model, x in CustomSet(1.0))\nx ∈ CustomSet{Float64}(1.0)\n\nHowever, there might be an appropriate sign that could be used in order to provide a more convenient syntax:\n\njulia> JuMP.sense_to_set(::Function, ::Val{:⊰}) = CustomSet(0.0)\n\njulia> MOIU.shift_constant(set::CustomSet, value) = CustomSet(set.value + value)\n\njulia> cref = @constraint(model, x ⊰ 1)\nx ∈ CustomSet{Float64}(1.0)\n\nNote that the whole function is first moved to the right-hand side, then the sign is transformed into a set with zero constant and finally the constant is moved to the set with MOIU.shift_constant.\n\n\n\n\n\n"
},

{
    "location": "extensions/#Adding-parse_constraint-methods-1",
    "page": "Extensions",
    "title": "Adding parse_constraint methods",
    "category": "section",
    "text": "# TODO(Benoît): Detail how `parse_constraint` works and show how `sense_to_set`\n#               fits into the picture.sense_to_set"
},

{
    "location": "extensions/#JuMP.build_constraint",
    "page": "Extensions",
    "title": "JuMP.build_constraint",
    "category": "function",
    "text": "function build_constraint(_error::Function, Q::Symmetric{V, M},\n                          ::PSDCone) where {V <: AbstractJuMPScalar,\n                                            M <: AbstractMatrix{V}}\n\nReturn a VectorConstraint of shape SymmetricMatrixShape constraining the matrix Q to be positive semidefinite.\n\nThis function is used by the @variable macro to create a symmetric semidefinite matrix of variables and by the @constraint macros as follows:\n\n@constraint(model, Symmetric(Q) in PSDCone())\n\nThe form above is usually used when the entries of Q are affine or quadratic expressions but it can also be used when the entries are variables to get the reference of the semidefinite constraint, e.g.,\n\n@variable model Q[1:2,1:2] Symmetric\n# The type of `Q` is `Symmetric{VariableRef, Matrix{VariableRef}}`\nvar_psd = @constraint model Q in PSDCone()\n# The `var_psd` variable contains a reference to the constraint\n\n\n\n\n\nfunction build_constraint(_error::Function,\n                          Q::AbstractMatrix{<:AbstractJuMPScalar},\n                          ::PSDCone)\n\nReturn a VectorConstraint of shape SquareMatrixShape constraining the matrix Q to be symmetric and positive semidefinite.\n\nThis function is used by the @constraint and @SDconstraint macros as follows:\n\n@constraint(model, Q in PSDCone())\n@SDconstraint(model, P ⪰ Q)\n\nThe @constraint call above is usually used when the entries of Q are affine or quadratic expressions but it can also be used when the entries are variables to get the reference of the semidefinite constraint, e.g.,\n\n@variable model Q[1:2,1:2]\n# The type of `Q` is `Matrix{VariableRef}`\nvar_psd = @constraint model Q in PSDCone()\n# The `var_psd` variable contains a reference to the constraint\n\n\n\n\n\n"
},

{
    "location": "extensions/#Adding-build_constraint-methods-1",
    "page": "Extensions",
    "title": "Adding build_constraint methods",
    "category": "section",
    "text": "There is typically two choices when creating a build_constraint method, either return an AbstractConstraint already supported by the model, i.e. ScalarConstraint or VectorConstraint, or a custom AbstractConstraint with a corresponding add_constraint method (see Adding add_constraint methods).build_constraint"
},

{
    "location": "extensions/#JuMP.AbstractShape",
    "page": "Extensions",
    "title": "JuMP.AbstractShape",
    "category": "type",
    "text": "AbstractShape\n\nAbstract vectorizable shape. Given a flat vector form of an object of shape shape, the original object can be obtained by reshape_vector.\n\n\n\n\n\n"
},

{
    "location": "extensions/#JuMP.shape",
    "page": "Extensions",
    "title": "JuMP.shape",
    "category": "function",
    "text": "shape(c::AbstractConstraint)::AbstractShape\n\nReturn the shape of the constraint c.\n\n\n\n\n\n"
},

{
    "location": "extensions/#JuMP.reshape_vector",
    "page": "Extensions",
    "title": "JuMP.reshape_vector",
    "category": "function",
    "text": "reshape_vector(vectorized_form::Vector, shape::AbstractShape)\n\nReturn an object in its original shape shape given its vectorized form vectorized_form.\n\nExamples\n\nGiven a SymmetricMatrixShape of vectorized form [1, 2, 3], the following code returns the matrix Symmetric(Matrix[1 2; 2 3]):\n\njulia> reshape_vector([1, 2, 3], SymmetricMatrixShape(2))\n2×2 LinearAlgebra.Symmetric{Int64,Array{Int64,2}}:\n 1  2\n 2  3\n\n\n\n\n\n"
},

{
    "location": "extensions/#JuMP.reshape_set",
    "page": "Extensions",
    "title": "JuMP.reshape_set",
    "category": "function",
    "text": "reshape_set(vectorized_set::MOI.AbstractSet, shape::AbstractShape)\n\nReturn a set in its original shape shape given its vectorized form vectorized_form.\n\nExamples\n\nGiven a SymmetricMatrixShape of vectorized form [1, 2, 3] in MOI.PositiveSemidefinieConeTriangle(2), the following code returns the set of the original constraint Symmetric(Matrix[1 2; 2 3]) in PSDCone():\n\njulia> reshape_set(MOI.PositiveSemidefiniteConeTriangle(2), SymmetricMatrixShape(2))\nPSDCone()\n\n\n\n\n\n"
},

{
    "location": "extensions/#JuMP.dual_shape",
    "page": "Extensions",
    "title": "JuMP.dual_shape",
    "category": "function",
    "text": "dual_shape(shape::AbstractShape)::AbstractShape\n\nReturns the shape of the dual space of the space of objects of shape shape. By default, the dual_shape of a shape is itself. See the examples section below for an example for which this is not the case.\n\nExamples\n\nConsider polynomial constraints for which the dual is moment constraints and moment constraints for which the dual is polynomial constraints. Shapes for polynomials can be defined as follows:\n\nstruct Polynomial\n    coefficients::Vector{Float64}\n    monomials::Vector{Monomial}\nend\nstruct PolynomialShape <: AbstractShape\n    monomials::Vector{Monomial}\nend\nJuMP.reshape_vector(x::Vector, shape::PolynomialShape) = Polynomial(x, shape.monomials)\n\nand a shape for moments can be defined as follows:\n\nstruct Moments\n    coefficients::Vector{Float64}\n    monomials::Vector{Monomial}\nend\nstruct MomentsShape <: AbstractShape\n    monomials::Vector{Monomial}\nend\nJuMP.reshape_vector(x::Vector, shape::MomentsShape) = Moments(x, shape.monomials)\n\nThe dual_shape allows to define the shape of the dual of polynomial and moment constraints:\n\ndual_shape(shape::PolynomialShape) = MomentsShape(shape.monomials)\ndual_shape(shape::MomentsShape) = PolynomialShape(shape.monomials)\n\n\n\n\n\n"
},

{
    "location": "extensions/#JuMP.ScalarShape",
    "page": "Extensions",
    "title": "JuMP.ScalarShape",
    "category": "type",
    "text": "ScalarShape\n\nShape of scalar constraints.\n\n\n\n\n\n"
},

{
    "location": "extensions/#JuMP.VectorShape",
    "page": "Extensions",
    "title": "JuMP.VectorShape",
    "category": "type",
    "text": "VectorShape\n\nVector for which the vectorized form corresponds exactly to the vector given.\n\n\n\n\n\n"
},

{
    "location": "extensions/#JuMP.SquareMatrixShape",
    "page": "Extensions",
    "title": "JuMP.SquareMatrixShape",
    "category": "type",
    "text": "SquareMatrixShape\n\nShape object for a square matrix of side_dimension rows and columns. The vectorized form contains the entries of the the matrix given column by column (or equivalently, the entries of the lower-left triangular part given row by row).\n\n\n\n\n\n"
},

{
    "location": "extensions/#JuMP.SymmetricMatrixShape",
    "page": "Extensions",
    "title": "JuMP.SymmetricMatrixShape",
    "category": "type",
    "text": "SymmetricMatrixShape\n\nShape object for a symmetric square matrix of side_dimension rows and columns. The vectorized form contains the entries of the upper-right triangular part of the matrix given column by column (or equivalently, the entries of the lower-left triangular part given row by row).\n\n\n\n\n\n"
},

{
    "location": "extensions/#Shapes-1",
    "page": "Extensions",
    "title": "Shapes",
    "category": "section",
    "text": "Shapes allow vector constraints, which are represented as flat vectors in MOI, to retain a matrix shape at the JuMP level. There is a shape field in VectorConstraint that can be set in build_constraint and that is used to reshape the result computed in value and dual.AbstractShape\nshape\nreshape_vector\nreshape_set\ndual_shape\nScalarShape\nVectorShape\nSquareMatrixShape\nSymmetricMatrixShape"
},

{
    "location": "extensions/#JuMP.add_constraint",
    "page": "Extensions",
    "title": "JuMP.add_constraint",
    "category": "function",
    "text": "add_constraint(model::Model, c::AbstractConstraint, name::String=\"\")\n\nAdd a constraint c to Model model and sets its name.\n\n\n\n\n\n"
},

{
    "location": "extensions/#Adding-add_constraint-methods-1",
    "page": "Extensions",
    "title": "Adding add_constraint methods",
    "category": "section",
    "text": "# TODO: Introduce `add_constraint`add_constraint### Extending the [`@objective`](@ref) macro\n\n# TODO: Describe how to `@objective` macro by implementing new `JuMP.set_objective_function` methods\n\n## Defining new JuMP models\n\n# TODO: Describe how to create a new JuMP model (similar to `test/JuMPExtension.jl` and StructJuMP)."
},

]}
