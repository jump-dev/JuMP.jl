var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#JuMP-–-Julia-for-Mathematical-Optimization-1",
    "page": "Introduction",
    "title": "JuMP –- Julia for Mathematical Optimization",
    "category": "section",
    "text": "JuMP is a domain-specific modeling language for mathematical optimization embedded in Julia. It currently supports a number of open-source and commercial solvers (see below) for a variety of problem classes, including linear programming, mixed-integer programming, second-order conic programming, semidefinite programming, and nonlinear programming. JuMP's features include:User friendliness\nSyntax that mimics natural mathematical expressions.\nComplete documentation.\nSpeed\nBenchmarking has shown that JuMP can create problems at similar speeds to special-purpose modeling languages such as AMPL.\nJuMP communicates with solvers in memory, avoiding the need to write intermediary files.\nSolver independence\nJuMP uses a generic solver-independent interface provided by the MathProgBase package, making it easy to change between a number of open-source and commercial optimization software packages (\"solvers\").\nCurrently supported solvers include Artelys Knitro, Bonmin, Cbc, Clp, Couenne, CPLEX, ECOS, FICO Xpress, GLPK, Gurobi, Ipopt, MOSEK, NLopt, and SCS.\nAccess to advanced algorithmic techniques\nIncluding efficient LP re-solves &lt;probmod&gt; and callbacks for mixed-integer programming &lt;callbacks&gt; which previously required using solver-specific and/or low-level C++ libraries.\nEase of embedding\nJuMP itself is written purely in Julia. Solvers are the only binary dependencies.\nBeing embedded in a general-purpose programming language makes it easy to solve optimization problems as part of a larger workflow (e.g., inside a simulation, behind a web server, or as a subproblem in a decomposition algorithm).\nAs a trade-off, JuMP's syntax is constrained by the syntax available in Julia.\nJuMP is MPL licensed, meaning that it can be embedded in commercial software that complies with the terms of the license.While neither Julia nor JuMP have reached version 1.0 yet, the releases are stable enough for everyday use and are being used in a number of research projects and neat applications by a growing community of users who are early adopters. JuMP remains under active development, and we welcome your feedback, suggestions, and bug reports."
},

{
    "location": "index.html#Installing-JuMP-1",
    "page": "Introduction",
    "title": "Installing JuMP",
    "category": "section",
    "text": "If you are familiar with Julia you can get started quickly by using the package manager to install JuMP:julia> Pkg.add(\"JuMP\")And a solver, e.g.:julia> Pkg.add(\"Clp\")  # Will install Cbc as wellThen read the quick-start and/or see a simple-example. The subsequent sections detail the complete functionality of JuMP."
},

{
    "location": "index.html#Contents-1",
    "page": "Introduction",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"installation.md\",\n    \"quickstart.md\",\n    \"refmodel.md\",\n    \"refvariable.md\",\n    \"refexpr.md\",\n    \"probmod.md\",\n    \"callbacks.md\",\n    \"nlp.md\"]\nDepth = 2"
},

{
    "location": "index.html#Citing-JuMP-1",
    "page": "Introduction",
    "title": "Citing JuMP",
    "category": "section",
    "text": "If you find JuMP useful in your work, we kindly request that you cite the following paper:@article{DunningHuchetteLubin2017,\nauthor = {Iain Dunning and Joey Huchette and Miles Lubin},\ntitle = {JuMP: A Modeling Language for Mathematical Optimization},\njournal = {SIAM Review},\nvolume = {59},\nnumber = {2},\npages = {295-320},\nyear = {2017},\ndoi = {10.1137/15M1020575},\n}A preprint of this paper is freely available on arXiv.For an earlier work where we presented a prototype implementation of JuMP, see here:@article{LubinDunningIJOC,\nauthor = {Miles Lubin and Iain Dunning},\ntitle = {Computing in Operations Research Using Julia},\njournal = {INFORMS Journal on Computing},\nvolume = {27},\nnumber = {2},\npages = {238-248},\nyear = {2015},\ndoi = {10.1287/ijoc.2014.0623},\n}A preprint of this paper is also freely available."
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
    "text": "This guide will briefly guide you through installing Julia, JuMP and[a] solver[s] of your choice."
},

{
    "location": "installation.html#Getting-Julia-1",
    "page": "Installation Guide",
    "title": "Getting Julia",
    "category": "section",
    "text": "At the time of writing this documentation the latest release of Julia is version 0.5, which is the version required by JuMP. You can easily build from source on OS X and Linux, but the binaries will work well for most people.Download links and more detailed instructions are available on the Julia website."
},

{
    "location": "installation.html#Getting-JuMP-1",
    "page": "Installation Guide",
    "title": "Getting JuMP",
    "category": "section",
    "text": "Once you've installed Julia, installing JuMP is simple. Julia has a git-based package system. To use it, open Julia in interactive mode (i.e. julia at the command line) and use the package manager:julia> Pkg.add(\"JuMP\")This command checks METADATA.jl to determine what the most recent version of JuMP is and then downloads it from its repository on GitHub.To start using JuMP (after installing a solver), it should be imported into the local scope:julia> using JuMP"
},

{
    "location": "installation.html#Getting-Solvers-1",
    "page": "Installation Guide",
    "title": "Getting Solvers",
    "category": "section",
    "text": "Solver support in Julia is currently provided by writing a solver-specific package that provides a very thin wrapper around the solver's C interface and providing a standard interface that JuMP can call. If you are interested in providing an interface to your solver, please get in touch. The table below lists the currently supported solvers and their capabilities.Solver Julia Package solver= License LP SOCP MILP NLP MINLP SDP\nArtelys Knitro KNITRO.jl KnitroSolver() Comm.    X X \nBARON BARON.jl BaronSolver() Comm.    X X \nBonmin AmplNLWriter.jl BonminNLSolver() * EPL X  X X X \n'' CoinOptServices.jl OsilBonminSolver() ''      \nCbc Cbc.jl CbcSolver() EPL   X   \nClp Clp.jl ClpSolver() EPL X     \nCouenne AmplNLWriter.jl CouenneNLSolver() * EPL X  X X X \n'' CoinOptServices.jl OsilCouenneSolver() ''      \nCPLEX CPLEX.jl CplexSolver() Comm. X X X   \nECOS ECOS.jl ECOSSolver() GPL X X    \nFICO Xpress Xpress.jl XpressSolver() Comm. X X X   \nGLPK GLPKMath... GLPKSolver[LP|MIP]() GPL X  X   \nGurobi Gurobi.jl GurobiSolver() Comm. X X X   \nIpopt Ipopt.jl IpoptSolver() EPL X   X  \nMOSEK Mosek.jl MosekSolver() Comm. X X X X  X\nNLopt NLopt.jl NLoptSolver() LGPL    X  \nSCS SCS.jl SCSSolver() MIT X X    XWhere:LP = Linear programming\nSOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)\nMILP = Mixed-integer linear programming\nNLP = Nonlinear programming\nMINLP = Mixed-integer nonlinear programming\nSDP = Semidefinite programming* requires CoinOptServices installed, see below.To install Gurobi, for example, and use it with a JuMP model m, run:Pkg.add(\"Gurobi\")\nusing JuMP\nusing Gurobi\n\nm = Model(solver=GurobiSolver())Setting solver options is discussed in the Model &lt;ref-model&gt; section.Solver-specific notes follow below."
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
    "location": "installation.html#COIN-OR-Clp-and-Cbc-1",
    "page": "Installation Guide",
    "title": "COIN-OR Clp and Cbc",
    "category": "section",
    "text": "Binaries for Clp and Cbc are provided on OS X and Windows (32- and 64-bit) by default. On Linux, they will be compiled from source (be sure to have a C++ compiler installed). Cbc supports \"SOS\" constraints but does not support MIP callbacks."
},

{
    "location": "installation.html#CPLEX-1",
    "page": "Installation Guide",
    "title": "CPLEX",
    "category": "section",
    "text": "Requires a working installation of CPLEX with a license (free for faculty members and graduate teaching assistants). The interface requires using CPLEX as a shared library, which is unsupported by the CPLEX developers. Special installation steps are required on OS X. CPLEX supports MIP callbacks and \"SOS\" constraints."
},

{
    "location": "installation.html#ECOS-1",
    "page": "Installation Guide",
    "title": "ECOS",
    "category": "section",
    "text": "ECOS can be used by JuMP to solve LPs and SOCPs. ECOS does not support general quadratic objectives or constraints, only second-order conic constraints specified by using norm or the quadratic form x'x <= y^2."
},

{
    "location": "installation.html#FICO-Xpress-1",
    "page": "Installation Guide",
    "title": "FICO Xpress",
    "category": "section",
    "text": "Requires a working installation of Xpress with an active license (it is possible to get license for academic use, see FICO Academic Partner Program). Supports SOCP and \"SOS\" constraints. The interface is experimental, but it does pass all JuMP and MathProgBase tests. Callbacks are not yet supported.warning: Warning\nIf you are using 64-bit Xpress, you must use 64-bit Julia (and similarly with 32-bit Xpress)."
},

{
    "location": "installation.html#GLPK-1",
    "page": "Installation Guide",
    "title": "GLPK",
    "category": "section",
    "text": "GLPK binaries are provided on OS X and Windows (32- and 64-bit) by default. On Linux, it will be compiled from source. Note that GLPKSolverLP should be used for continuous problems and GLPKSolverMIP for problems with integer variables. GLPK supports MIP callbacks but does not support \"SOS\" constraints."
},

{
    "location": "installation.html#Gurobi-1",
    "page": "Installation Guide",
    "title": "Gurobi",
    "category": "section",
    "text": "Requires a working installation of Gurobi with an activated license (free for academic use). Gurobi supports MIP callbacks and \"SOS\" constraints.warning: Warning\nIf you are using 64-bit Gurobi, you must use 64-bit Julia (and similarly with 32-bit Gurobi)."
},

{
    "location": "installation.html#Ipopt-1",
    "page": "Installation Guide",
    "title": "Ipopt",
    "category": "section",
    "text": "Ipopt binaries are provided on OS X and Windows (32- and 64-bit) by default. On Linux, it will be compiled from source. The default installation of Ipopt uses the open-source MUMPS library for sparse linear algebra. Significant speedups can be obtained by manually compiling Ipopt to use proprietary sparse linear algebra libraries instead. Julia can be pointed to use a custom version of Ipopt; we suggest posting to the julia-opt mailing list with your platform details for guidance on how to do this."
},

{
    "location": "installation.html#MOSEK-1",
    "page": "Installation Guide",
    "title": "MOSEK",
    "category": "section",
    "text": "Requires a license (free for academic use). Mosek does not support the MIP callbacks used in JuMP. For nonlinear optimization, Mosek supports only convex problems. The Mosek interface is maintained by the Mosek team. (Thanks!)"
},

{
    "location": "installation.html#NLopt-1",
    "page": "Installation Guide",
    "title": "NLopt",
    "category": "section",
    "text": "NLopt supports only nonlinear models. An algorithm must be specified as an option when using NLoptSolver. NLopt is not recommended for large-scale models, because it does not currently exploit sparsity of derivative matrices."
},

{
    "location": "installation.html#SCS-1",
    "page": "Installation Guide",
    "title": "SCS",
    "category": "section",
    "text": "SCS can be used by JuMP to solve LPs and SOCPs, and SDPs. SCS is a first order solver and has low accuracy (10^4) by default; see the SCS.jl documentation for more information."
},

{
    "location": "installation.html#COIN-OR-Bonmin-and-Couenne-1",
    "page": "Installation Guide",
    "title": "COIN-OR Bonmin and Couenne",
    "category": "section",
    "text": "Binaries of Bonmin and Couenne are provided on OS X and Windows (32- and 64-bit) by the CoinOptServices.jl package. On Linux, they will be compiled from source. Once installed, they can be called either via .osil files using OsilBonminSolver and OsilCouenneSolver from CoinOptServices.jl, or via .nl files using BonminNLSolver and CouenneNLSolver from AmplNLWriter.jl. We recommend using the .nl format option, which is currently more stable and has better performance for derivative computations. Since both Bonmin and Couenne use Ipopt for continuous subproblems, the same MUMPS sparse linear algebra performance caveat applies."
},

{
    "location": "installation.html#Other-AMPL-compatible-solvers-1",
    "page": "Installation Guide",
    "title": "Other AMPL-compatible solvers",
    "category": "section",
    "text": "Any other solver not listed above that can be called from AMPL can be used by JuMP through the AmplNLWriter.jl package. The first argument to AmplNLSolver can be used to specify a solver executable name.For example, SCIP is a powerful noncommercial mixed-integer programming solver. To use SCIP within JuMP, you must first download and compile SCIP with support for AMPL. Then you may use AmplNLSolver(\"/path/to/scipampl\") where scipampl is the executable produced from the compilation process."
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
    "text": "This quick start guide will introduce the main concepts of JuMP. If you are familiar with another modeling language embedded in a high-level language such as PuLP (Python) or a solver-specific interface you will find most of this familiar, with the exception of macros. A deep understanding of macros is not essential, but if you would like to know more please see the Julia documentation. If you are coming from an AMPL or similar background, you may find some of the concepts novel but the general appearance will still be familiar."
},

{
    "location": "quickstart.html#Creating-a-Model-1",
    "page": "Quick Start Guide",
    "title": "Creating a Model",
    "category": "section",
    "text": "Models are Julia objects. They are created by calling the constructor:m = Model()All variables and constraints are associated with a Model object. Usually, you'll also want to provide a solver object here by using the solver= keyword argument; see the simple example below. For a list of all functions related to Model, see Models."
},

{
    "location": "quickstart.html#Defining-Variables-1",
    "page": "Quick Start Guide",
    "title": "Defining Variables",
    "category": "section",
    "text": "Variables are also Julia objects, and are defined using the @variable macro. The first argument will always be the Model to associate this variable with. In the examples below we assume m is already defined. The second argument is an expression that declares the variable name and optionally allows specification of lower and upper bounds. For example:@variable(m, x )              # No bounds\n@variable(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)\n@variable(m, x <= ub )        # Upper bound only\n@variable(m, lb <= x <= ub )  # Lower and upper boundsAll these variations introduce a new variable x in the local scope. The names of your variables must be valid Julia variable names. For information about common operations on variables, e.g. changing their bounds, see the Variables section.Integer and binary restrictions can optionally be specified with a third argument, Int or Bin.To create arrays of variables we append brackets to the variable name. For example:@variable(m, x[1:M,1:N] >= 0 )will create an M by N array of variables. Both ranges and arbitrary iterable sets are supported as index sets. Currently we only support ranges of the form a:b where a is an explicit integer, not a variable. Using ranges will generally be faster than using arbitrary symbols. You can mix both ranges and lists of symbols, as in the following example:s = [\"Green\", \"Blue\"]\n@variable(m, x[-10:10,s], Int )\n# e.g. x[-4, \"Green\"]Finally, bounds can depend on variable indices:@variable(m, x[i=1:10] >= i )"
},

{
    "location": "quickstart.html#Objective-and-Constraints-1",
    "page": "Quick Start Guide",
    "title": "Objective and Constraints",
    "category": "section",
    "text": "JuMP allows users to use a natural notation to describe linear expressions. To add constraints, use the @constraint() and @objective() macros, e.g.:@constraint(m, x[i] - s[i] <= 0)  # Other options: == and >=\n@constraint(m, sum(x[i] for i=1:numLocation) == 1)\n@objective(m, Max, 5x + 22y + (x+y)/2) # or Minnote: Note\nThe sense passed to @objective must be a symbol type: :Min or :Max, although the macro accepts :Min and :Max, as well as Min and Max (without the colon) directly.The sum() syntax directly follows Julia's own generator expression syntax. You may use conditions within sums, e.g.:sum(expression for i = I1, j = I2 if cond)which is equivalent to:a = zero(AffExpr)\nfor i = I1\n    for j = I2\n        ...\n        if cond\n            a += expression\n        end\n        ...\n    end\nendnote: Note\nJuMP previously used a special curly brace syntax for sum{}, prod{}, and norm2{}. This has been entirely replaced by sum(), prod(), and norm() since Julia 0.5. The curly brace syntax is deprecated and will be removed in a future release."
},

{
    "location": "quickstart.html#Simple-Example-1",
    "page": "Quick Start Guide",
    "title": "Simple Example",
    "category": "section",
    "text": "In this section we will construct a simple model and explain every step along the way. The are more complex examples in the JuMP/examples/ folder. Here is the code we will walk through:using JuMP\nusing Clp\n\nm = Model(solver = ClpSolver())\n@variable(m, 0 <= x <= 2 )\n@variable(m, 0 <= y <= 30 )\n\n@objective(m, Max, 5x + 3*y )\n@constraint(m, 1x + 5y <= 3.0 )\n\nprint(m)\n\nstatus = solve(m)\n\nprintln(\"Objective value: \", getobjectivevalue(m))\nprintln(\"x = \", getvalue(x))\nprintln(\"y = \", getvalue(y))Once JuMP is installed Installation Guide, to use JuMP in your programs, you just need to say:using JuMPWe also need to include a Julia package which provides an appropriate solver. In this case, we'll use Clp:using ClpModels are created with the Model() function. The solver= keyword argument is used to specify the solver to be used:m = Model(solver = ClpSolver())note: Note\nYour model doesn't have to be called m - it's just a name.There are a few options for defining a variable, depending on whether you want to have lower bounds, upper bounds, both bounds, or even no bounds. The following commands will create two variables, x and y, with both lower and upper bounds. Note the first argument is our model variable m. These variables are associated with this model and cannot be used in another model.:@variable(m, 0 <= x <= 2 )\n@variable(m, 0 <= y <= 30 )Next we'll set our objective. Note again the m, so we know which model's objective we are setting! The objective sense, Max or Min, should be provided as the second argument. Note also that we don't have a multiplication * symbol between 5 and our variable x - Julia is smart enough to not need it! Feel free to stick with * if it makes you feel more comfortable, as we have done with 3*y:@objective(m, Max, 5x + 3*y )Adding constraints is a lot like setting the objective. Here we create a less-than-or-equal-to constraint using <=, but we can also create equality constraints using == and greater-than-or-equal-to constraints with >=:@constraint(m, 1x + 5y <= 3.0 )If you want to see what your model looks like in a human-readable format, the print function is defined for models.print(m)Models are solved with the solve() function. This function will not raise an error if your model is infeasible - instead it will return a flag. In this case, the model is feasible so the value of status will be :Optimal, where : again denotes a symbol. The possible values of status are described in Solve Status.status = solve(m)Finally, we can access the results of our optimization. Getting the objective value is simple:println(\"Objective value: \", getobjectivevalue(m))To get the value from a variable, we call the getvalue() function. If x is not a single variable, but instead a range of variables, getvalue() will return a list. In this case, however, it will just return a single value.println(\"x = \", getvalue(x))\nprintln(\"y = \", getvalue(y))"
},

{
    "location": "refmodel.html#",
    "page": "Models",
    "title": "Models",
    "category": "page",
    "text": "CurrentModule = JuMP"
},

{
    "location": "refmodel.html#Models-1",
    "page": "Models",
    "title": "Models",
    "category": "section",
    "text": ""
},

{
    "location": "refmodel.html#Constructor-1",
    "page": "Models",
    "title": "Constructor",
    "category": "section",
    "text": "Model is a type defined by JuMP. All variables and constraints are associated with a Model object. It has a constructor that has no required arguments:m = Model()The constructor also accepts an optional keyword argument, solver. You may specify a solver either here or later on by calling setsolver. JuMP will throw an error if you try to solve a problem without specifying a solver.solver must be an AbstractMathProgSolver object, which is constructed as follows:solver = solvername(Option1=Value1, Option2=Value2, ...)where solvername is one of the supported solvers. See the solver table &lt;jump-solvertable&gt; for the list of available solvers and corresponding parameter names. All options are solver-dependent; see corresponding solver packages for more information.note: Note\nBe sure that the solver provided supports the problem class of the model. For example ClpSolver and GLPKSolverLP support only linear programming problems. CbcSolver and GLPKSolverMIP support only mixed-integer programming problems.As an example, we can create a Model object that will use GLPK's exact solver for LPs as follows:m = Model(solver = GLPKSolverLP(method=:Exact))"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.numvar",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.numvar",
    "category": "Function",
    "text": "MathProgBase.numvar(m::Model)\n\nreturns the number of variables associated with the Model m.\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.numlinconstr",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.numlinconstr",
    "category": "Function",
    "text": "MathProgBase.numlinconstr(m::Model)\n\nreturns the number of linear constraints associated with the Model m\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.numquadconstr",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.numquadconstr",
    "category": "Function",
    "text": "MathProgBase.numquadconstr(m::Model)\n\nreturns the number of quadratic constraints associated with the Model m\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.numsocconstr",
    "page": "Models",
    "title": "JuMP.numsocconstr",
    "category": "Function",
    "text": "numsocconstr(m::Model)\n\nreturns the number of second order cone constraints associated with the Model m\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.numsosconstr",
    "page": "Models",
    "title": "JuMP.numsosconstr",
    "category": "Function",
    "text": "numsosconstr(m::Model)\n\nreturns the number of sos constraints associated with the Model m\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.numsdconstr",
    "page": "Models",
    "title": "JuMP.numsdconstr",
    "category": "Function",
    "text": "numsdconstr(m::Model)\n\nreturns the number of semi-definite constraints associated with the Model m\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.numnlconstr",
    "page": "Models",
    "title": "JuMP.numnlconstr",
    "category": "Function",
    "text": "numnlconstr(m::Model)\n\nreturns the number of nonlinear constraints associated with the Model m\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.numconstr",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.numconstr",
    "category": "Function",
    "text": "MathProgBase.numconstr(m::Model)\n\nreturns the total number of constraints associated with the Model m\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.internalmodel",
    "page": "Models",
    "title": "JuMP.internalmodel",
    "category": "Function",
    "text": "internalmodel(m::Model)\n\nreturns the internal low-level AbstractMathProgModel object which can be used to access any functionality that is not exposed by JuMP. See the MathProgBase documentation\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.solve",
    "page": "Models",
    "title": "JuMP.solve",
    "category": "Function",
    "text": "solve(m::Model; suppress_warnings=false,\n            ignore_solve_hook=(m.solvehook===nothing),\n            relaxation=false,\n            kwargs...)\n\nsolves the model using the selected solver (or a default for the problem class), and takes two optional arguments that are disabled by default. Setting suppress_warnings to true will suppress all JuMP-specific output (e.g. warnings about infeasibility and lack of dual information) but will not suppress solver output (which should be done by passing options to the solver). Setting relaxation=true solves the standard continuous relaxation for the model: that is, integrality is dropped, special ordered set constraints are not enforced, and semi-continuous and semi-integer variables with bounds [l,u] are replaced with bounds [min(l,0),max(u,0)]\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.build",
    "page": "Models",
    "title": "JuMP.build",
    "category": "Function",
    "text": "build(m::Model; suppress_warnings=false, \n    relaxation=false, \n    traits=ProblemTraits(m,relaxation=relaxation))\n\nbuilds the model in memory at the MathProgBase level without optimizing.\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.setsolver",
    "page": "Models",
    "title": "JuMP.setsolver",
    "category": "Function",
    "text": "setsolver(m::Model, solver::MathProgBase.AbstractMathProgSolver)\n\nchanges the solver which will be used for the next call to solve(), discarding the current internal model if present.\n\n\n\n"
},

{
    "location": "refmodel.html#Base.getindex-Tuple{JuMP.Model,Symbol}",
    "page": "Models",
    "title": "Base.getindex",
    "category": "Method",
    "text": "Base.getindex(m::JuMP.Model, name::Symbol)\n\nTo allow easy accessing of JuMP Variables and Constraints via [] syntax. Returns the variable, or group of variables, or constraint, or group of constraints, of the given name which were added to the model. This errors if multiple variables or constraints share the same name.\n\n\n\n"
},

{
    "location": "refmodel.html#Base.setindex!-Tuple{JuMP.Model,Any,Symbol}",
    "page": "Models",
    "title": "Base.setindex!",
    "category": "Method",
    "text": "Base.setindex!(m::JuMP.Model, value, name::Symbol)\n\nstores the object value in the model m using so that it can be accessed via getindex.  Can be called with [] syntax.\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.getobjective",
    "page": "Models",
    "title": "JuMP.getobjective",
    "category": "Function",
    "text": "getobjective(m::Model)\n\nreturns the objective function as a QuadExpr\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.getobjectivesense",
    "page": "Models",
    "title": "JuMP.getobjectivesense",
    "category": "Function",
    "text": "getobjectivesense(m::Model)\n\nreturns objective sense, either :Min or :Max\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.setobjectivesense",
    "page": "Models",
    "title": "JuMP.setobjectivesense",
    "category": "Function",
    "text": "setobjectivesense(m::Model, newSense::Symbol)\n\nsets the objective sense (newSense is either :Min or :Max)\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.getobjectivevalue",
    "page": "Models",
    "title": "JuMP.getobjectivevalue",
    "category": "Function",
    "text": "getobjectivevalue(m::Model)\n\nreturns objective value after a call to solve\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.getobjectivebound",
    "page": "Models",
    "title": "JuMP.getobjectivebound",
    "category": "Function",
    "text": "getobjectivebound(m::Model)\n\nreturns the best known bound on the optimal objective value after a call to solve\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.getsolvetime",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.getsolvetime",
    "category": "Function",
    "text": "getsolvetime(m::Model)\n\nreturns the solve time reported by the solver if it is implemented.\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.getnodecount",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.getnodecount",
    "category": "Function",
    "text": "getnodecount(m::Model)\n\nreturns the number of explored branch-and-bound nodes, if it is implemented.\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.getobjbound",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.getobjbound",
    "category": "Function",
    "text": "getobjbound(m::Model)\n\nreturns the best known bound on the optimal objective value. This is used, for example, when a branch-and-bound method is stopped before finishing.\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.getobjgap",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.getobjgap",
    "category": "Function",
    "text": "getobjgap(m::Model)\n\nreturns the final relative optimality gap as optimization terminated. That is, it returns fracb-ff, where b is the best bound and f is the best feasible objective value.\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.getrawsolver",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.getrawsolver",
    "category": "Function",
    "text": "getrawsolver(m::Model)\n\nreturns an object that may be used to access a solver-specific API.\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.getsimplexiter",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.getsimplexiter",
    "category": "Function",
    "text": "getsimplexiter(m::Model)\n\nreturns the cumulative number of simplex iterations during the optimization process. In particular, for a MIP it returns the total simplex iterations for all nodes.\n\n\n\n"
},

{
    "location": "refmodel.html#MathProgBase.SolverInterface.getbarrieriter",
    "page": "Models",
    "title": "MathProgBase.SolverInterface.getbarrieriter",
    "category": "Function",
    "text": "getbarrieriter(m::Model)\n\nreturns the cumulative number of barrier iterations during the optimization process.\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.writeLP",
    "page": "Models",
    "title": "JuMP.writeLP",
    "category": "Function",
    "text": "writeLP(m::Model, fname::AbstractString; genericnames=true)\n\nwrite the model to fname in the LP file format. Set genericnames=false for user-defined variable names.\n\n\n\n"
},

{
    "location": "refmodel.html#JuMP.writeMPS",
    "page": "Models",
    "title": "JuMP.writeMPS",
    "category": "Function",
    "text": "writeMPS(m::Model, fname::AbstractString)\n\nwrite the model to fname in the MPS file format.\n\n\n\n"
},

{
    "location": "refmodel.html#Methods-1",
    "page": "Models",
    "title": "Methods",
    "category": "section",
    "text": "GeneralMathProgBase.numvar\nMathProgBase.numlinconstr\nMathProgBase.numquadconstr\nnumsocconstr\nnumsosconstr\nnumsdconstr\nnumnlconstr\nMathProgBase.numconstr\ninternalmodel\nsolve\nbuild\nsetsolver\nBase.getindex(m::Model, name::Symbol)\nBase.setindex!(m::Model, value, name::Symbol)Objectivegetobjective\ngetobjectivesense\nsetobjectivesense\ngetobjectivevalue\ngetobjectiveboundSolverThese functions are JuMP versions of the similarly named functions in MathProgBase.getsolvetime\ngetnodecount\ngetobjbound\ngetobjgap\ngetrawsolver\ngetsimplexiter\ngetbarrieriterOutputwriteLP\nwriteMPS"
},

{
    "location": "refmodel.html#Solve-Status-1",
    "page": "Models",
    "title": "Solve Status",
    "category": "section",
    "text": "The call status = solve(m) returns a symbol recording the status of the optimization process, as reported by the solver. Typical values are listed in the table below, although the code can take solver-dependent values. For instance, certain solvers prove infeasibility or unboundedness during presolve, but do not report which of the two cases holds. See your solver interface documentation (as linked to in the solver table &lt;jump-solvertable&gt;) for more information.Status Meaning\n:Optimal Problem solved to optimality\n:Unbounded Problem is unbounded\n:Infeasible Problem is infeasible\n:UserLimit Iteration limit or timeout\n:Error Solver exited with an error\n:NotSolved Model built in memory but not optimized"
},

{
    "location": "refmodel.html#Quadratic-Objectives-1",
    "page": "Models",
    "title": "Quadratic Objectives",
    "category": "section",
    "text": "Quadratic objectives are supported by JuMP using a solver which implements the corresponding extensions of the MathProgBase interface. Add them in the same way you would a linear objective:using Ipopt\nm = Model(solver=IpoptSolver())\n@variable(m, 0 <= x <= 2 )\n@variable(m, 0 <= y <= 30 )\n\n@objective(m, Min, x*x+ 2x*y + y*y )\n@constraint(m, x + y >= 1 )\n\nprint(m)\n\nstatus = solve(m)"
},

{
    "location": "refmodel.html#Second-order-cone-constraints-1",
    "page": "Models",
    "title": "Second-order cone constraints",
    "category": "section",
    "text": "Second-order cone constraints of the form Axb_2+a^Tx+cleq0 can be added directly using the norm function:@constraint(m, norm(A*x) <= 2w - 1)You may use generator expressions within norm() to build up normed expressions with complex indexing operations in much the same way as with sum(...):@constraint(m, norm(2x[i] - i for i=1:n if c[i] == 1) <= 1)"
},

{
    "location": "refmodel.html#Accessing-the-low-level-model-1",
    "page": "Models",
    "title": "Accessing the low-level model",
    "category": "section",
    "text": "It is possible to construct the internal low-level model before optimizing. To do this, call the JuMP.build function. It is then possible to obtain this model by using the internalmodel function. This may be useful when it is necessary to access some functionality that is not exposed by JuMP. When you are ready to optimize, simply call solve in the normal fashion."
},

{
    "location": "refvariable.html#",
    "page": "Variables",
    "title": "Variables",
    "category": "page",
    "text": "CurrentModule = JuMP"
},

{
    "location": "refvariable.html#Variables-1",
    "page": "Variables",
    "title": "Variables",
    "category": "section",
    "text": "Variables, also known as columns or decision variables, are the results of the optimization."
},

{
    "location": "refvariable.html#Constructors-1",
    "page": "Variables",
    "title": "Constructors",
    "category": "section",
    "text": "The primary way to create variables is with the @variable macro. The first argument will always be a Model. In the examples below we assume m is already defined. The second argument is an expression that declares the variable name and optionally allows specification of lower and upper bounds. Adding variables \"column-wise\", e.g., as in column generation, is supported as well; see the syntax discussed in the Problem Modification section.@variable(m, x )              # No bounds\n@variable(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)\n@variable(m, x <= ub )        # Upper bound only\n@variable(m, lb <= x <= ub )  # Lower and upper bounds\n@variable(m, x == fixedval )  # Fixed to a value (lb == ub)All these variations create a new local variable, in this case x. The names of your variables must be valid Julia variable names. Integer and binary restrictions can optionally be specified with a third argument, Int or Bin. For advanced users, SemiCont and SemiInt may be used to create semicontinuous or semi-integer variables, respectively.To create arrays of variables we append brackets to the variable name.@variable(m, x[1:M,1:N] >= 0 )will create an M by N array of variables. Both ranges and arbitrary iterable sets are supported as index sets. Currently we only support ranges of the form a:b where a is an explicit integer, not a variable. Using ranges will generally be faster than using arbitrary symbols. You can mix both ranges and lists of symbols, as in the following example:s = [\"Green\",\"Blue\"]\n@variable(m, x[-10:10,s] , Int)\nx[-4,\"Green\"]Bounds can depend on variable indices:@variable(m, x[i=1:10] >= i )And indices can have dependencies on preceding indices (e.g. \"triangular indexing\"):@variable(m, x[i=1:10,j=i:10] >= 0)Note the dependency must be on preceding indices, going from left to right. That is, @variable(m, x[i=j:10,i=1:10] >= 0) is not valid JuMP code.Conditions can be placed on the index values for which variables are created; the condition follows the statement of the index sets and is separated with a semicolon:@variable(m, x[i=1:10,j=1:10; isodd(i+j)] >= 0)Note that only one condition can be added, although expressions can be built up by using the usual && and || logical operators.An initial value of each variable may be provided with the start keyword to @variable:@variable(m, x[i=1:10], start=(i/2))Is equivalent to:@variable(m, x[i=1:10])\nfor i in 1:10\n    setvalue(x[i], i/2)\nendFor more complicated variable bounds, it may be clearer to specify them using the lowerbound and upperbound keyword arguments to @variable:@variable(m, x[i=1:3], lowerbound=my_complex_function(i))\n@variable(m, x[i=1:3], lowerbound=my_complex_function(i), upperbound=another_function(i))Variable categories may be set in a more programmatic way by providing the appropriate symbol to the category keyword argument:t = [:Bin,:Int]\n@variable(m, x[i=1:2], category=t[i])\n@variable(m, y, category=:SemiCont)The constructor Variable(m::Model,idx::Int) may be used to create a variable object corresponding to an existing variable in the model (the constructor does not add a new variable to the model). The variable indices correspond to those of the internal MathProgBase model. The inverse of this operation is linearindex(x::Variable), which returns the flattened out (linear) index of a variable as JuMP provides it to a solver. We guarantee that Variable(m,linearindex(x)) returns x itself. These methods are only useful if you intend to interact with solver properties which are not directly exposed through JuMP.note: Note\n@variable is equivalent to a simple assignment x = ... in Julia and therefore redefines variables. The following code will generate a warning and may lead to unexpected results:@variable(m, x[1:10,1:10])\n@variable(m, x[1:5])After the second line, the Julia variable x refers to a set of variables indexed by the range 1:5. The reference to the first set of variables has been lost, although they will remain in the model. See also the section on anonymous variables."
},

{
    "location": "refvariable.html#Anonymous-variables-1",
    "page": "Variables",
    "title": "Anonymous variables",
    "category": "section",
    "text": "We also provide a syntax for constructing \"anonymous\" variables. In @variable, you may omit the name of the variable and instead assign the return value as you would like:x = @variable(m) # Equivalent to @variable(m, x)\nx = @variable(m, [i=1:3], lowerbound = i, upperbound = 2i) # Equivalent to @variable(m, i <= x[i=1:3] <= 2i)The lowerbound and upperbound keywords must be used instead of comparison operators for specifying variable bounds within the anonymous syntax. For creating noncontinuous anonymous variables, the category keyword must be used to avoid ambiguity, e.g.:x = @variable(m, Bin) # error\nx = @variable(m, category = :Bin) # okBesides these syntax restrictions in the @variable macro, the only differences between anonymous and named variables are:For the purposes of printing a model, JuMP will not have a name for anonymous variables and will instead use __anon__. You may set the name of a variable for printing by using setname or the basename keyword argument described below.\nAnonymous variables cannot be retrieved by using getindex or m[name].If you would like to change the name used when printing a variable or group of variables, you may use the basename keyword argument:i = 3\n@variable(m, x[1:3], basename=\"myvariable-$i\")\n# OR:\nx = @variable(m, [1:3], basename=\"myvariable-$i\")Printing x[2] will display myvariable-3[2]."
},

{
    "location": "refvariable.html#Semidefinite-and-symmetric-variables-1",
    "page": "Variables",
    "title": "Semidefinite and symmetric variables",
    "category": "section",
    "text": "JuMP supports modeling with semidefinite variables. A square symmetric matrix X is positive semidefinite if all eigenvalues are nonnegative; this is typically denoted by Xsucceq0. You can declare a matrix of variables to be positive semidefinite as follows:@variable(m, X[1:3,1:3], SDP)Note in particular the indexing: 1) exactly two index sets must be specified, 2) they must both be unit ranges starting at 1, 3) no bounds can be provided alongside the SDP tag. If you wish to impose more complex semidefinite constraints on the variables, e.g. XIsucceq0, you may instead use the Symmetric tag, along with a semidefinite constraint:@variable(m, X[1:n,1:n], Symmetric)\n@SDconstraint(m, X >= eye(n))Bounds can be provided as normal when using the Symmetric tag, with the stipulation that the bounds are symmetric themselves."
},

{
    "location": "refvariable.html#@variables-blocks-1",
    "page": "Variables",
    "title": "@variables blocks",
    "category": "section",
    "text": "JuMP provides a convenient syntax for defining multiple variables in a single block:@variables m begin\n    x\n    y >= 0\n    Z[1:10], Bin\n    X[1:3,1:3], SDP\n    q[i=1:2], (lowerbound = i, start = 2i, upperbound = 3i)\n    t[j=1:3], (Int, start = j)\nend\n\n# Equivalent to:\n@variable(m, x)\n@variable(m, y >= 0)\n@variable(m, Z[1:10], Bin)\n@variable(m, X[1:3,1:3], SDP)\n@variable(m, q[i=1:2], lowerbound = i, start = 2i, upperbound = 3i)\n@variable(m, t[j=1:3], Int, start = j)The syntax follows that of @variable with each declaration separated by a new line. Note that unlike in @variable, keyword arguments must be specified within parentheses."
},

{
    "location": "refvariable.html#JuMP.setlowerbound",
    "page": "Variables",
    "title": "JuMP.setlowerbound",
    "category": "Function",
    "text": "setlowerbound(v::Variable,lower::Number)\n\nset the lower bound of a variable.\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.getlowerbound",
    "page": "Variables",
    "title": "JuMP.getlowerbound",
    "category": "Function",
    "text": "getlowerbound(v::Variable)\n\nget the lower bound of a variable.\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.setupperbound",
    "page": "Variables",
    "title": "JuMP.setupperbound",
    "category": "Function",
    "text": "setupperbound(v::Variable,upper::Number)\n\nset the upper bound of a variable.\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.getupperbound",
    "page": "Variables",
    "title": "JuMP.getupperbound",
    "category": "Function",
    "text": "getupperbound(v::Variable)\n\nget the upper bound of a variable.\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.setcategory",
    "page": "Variables",
    "title": "JuMP.setcategory",
    "category": "Function",
    "text": "setcategory(v::Variable, cat::Symbol)\n\nSet the variable category for v after construction. Possible categories are :Cont, :Int, :Bin, :SemiCont, :SemiInt.\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.getcategory",
    "page": "Variables",
    "title": "JuMP.getcategory",
    "category": "Function",
    "text": "getcategory(v::Variable)\n\nGet the variable category for v\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.getvalue-Tuple{JuMP.Variable}",
    "page": "Variables",
    "title": "JuMP.getvalue",
    "category": "Method",
    "text": "getvalue(v::Variable)\n\nGet the value of this variable in the solution.\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.getvalue-Tuple{Array{JuMP.Variable,N} where N}",
    "page": "Variables",
    "title": "JuMP.getvalue",
    "category": "Method",
    "text": "getvalue(arr::Array{Variable})\n\nReturns an indexable dictionary of values. When the model is unbounded, returns the corresponding components of an unbounded ray, if available from the solver.\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.setvalue",
    "page": "Variables",
    "title": "JuMP.setvalue",
    "category": "Function",
    "text": "setvalue(v::Variable, val::Number)\n\nProvide an initial value v for this variable that can be used by supporting MILP solvers. If v is NaN, the solver may attempt to fill in this value to construct a feasible solution. setvalue cannot be used with fixed variables; instead their value may be set with JuMP.fix(x,v)\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.getdual",
    "page": "Variables",
    "title": "JuMP.getdual",
    "category": "Function",
    "text": "getdual(v::Variable)\n\nGet the reduced cost of this variable in the solution. Similar behavior to getvalue for indexable variables.\n\n\n\ngetdual(c::ConstraintRef{Model,SDConstraint})\n\n\n\ngetdual(c::LinConstrRef)\n\n\n\ngetdual(c::ConstraintRef{Model,SOCConstraint})\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.setname",
    "page": "Variables",
    "title": "JuMP.setname",
    "category": "Function",
    "text": "setname(v::Variable,n::AbstractString)\n\nSet the variable's internal name.\n\n\n\n"
},

{
    "location": "refvariable.html#JuMP.getname",
    "page": "Variables",
    "title": "JuMP.getname",
    "category": "Function",
    "text": "getname(m::Model, col)\n\nGet the variable's internal name.\n\n\n\ngetname(v::Variable)\n\nGet the variable's internal name.\n\n\n\n"
},

{
    "location": "refvariable.html#Methods-1",
    "page": "Variables",
    "title": "Methods",
    "category": "section",
    "text": "Boundssetlowerbound\ngetlowerbound\nsetupperbound\ngetupperboundVariable Categorysetcategory\ngetcategoryValuesgetvalue(v::Variable)\ngetvalue(arr::Array{Variable})\nsetvalue\ngetdualnote: Note\nThe getvalue function always returns a floating-point value, even when a variable is constrained to take integer values, as most solvers only guarantee integrality up to a particular numerical tolerance. The built-in round function should be used to obtain integer values, e.g., by calling round(Integer, getvalue(x)).NamesVariables (in the sense of columns) can have internal names (different from the Julia variable name) that can be used for writing models to file. This feature is disabled for performance reasons, but will be added if there is demand or a special use case.setname\ngetname"
},

{
    "location": "refvariable.html#Helper-functions-1",
    "page": "Variables",
    "title": "Helper functions",
    "category": "section",
    "text": "The following built-in functions are overloaded to provide easy construction of expressions from variables,sum(x) - Operates on arrays of variables, efficiently produces an affine expression. Available in macros.\ndot(x, coeffs) - Performs a generalized \"dot product\" for arrays of variables and coefficients up to three dimensions, or equivalently the sum of the elements of the Hadamard product. Available in macros, and also as dot(coeffs, x)."
},

{
    "location": "refvariable.html#Fixed-variables-1",
    "page": "Variables",
    "title": "Fixed variables",
    "category": "section",
    "text": "Fixed variables, created with the x == fixedval syntax, have slightly special semantics. First, it is important to note that fixed variables are considered optimization variables, not constants, for the purpose of determining the problem class. For example, in:@variable(m, x == 5)\n@variable(m, y)\n@constraint(m, x*y <= 10)the constraint added is a nonconvex quadratic constraint. For efficiency reasons, JuMP will not substitute the constant 5 for x and then provide the resulting linear constraint to the solver. Two possible uses for fixed variables are:For computing sensitivities. When available from the solver, the sensitivity of the objective with respect to the fixed value may be queried with getdual(x).\nFor solving a sequence of problems with varying parameters. One may call JuMP.fix(x, val) to change the value of a fixed variable or to fix a previously unfixed variable. For LPs in particular, most solvers are able to efficiently hot-start when solving the resulting modified problem."
},

{
    "location": "refexpr.html#",
    "page": "Expressions and Constraints",
    "title": "Expressions and Constraints",
    "category": "page",
    "text": "CurrentModule = JuMP"
},

{
    "location": "refexpr.html#Expressions-and-Constraints-1",
    "page": "Expressions and Constraints",
    "title": "Expressions and Constraints",
    "category": "section",
    "text": ""
},

{
    "location": "refexpr.html#Constructor-1",
    "page": "Expressions and Constraints",
    "title": "Constructor",
    "category": "section",
    "text": "AffExpr is an affine expression type defined by JuMP. It has three fields: a vector of coefficients, a vector of variables, and a constant. Apart from a default constructor that takes no arguments, it also has a full constructor that can be useful if you want to manually build an affine expression:aff = AffExpr([x, z], [3.0, 4.0], 2.0)  # 3x + 4z + 2Note that the coefficients must be floating point numbers. The matching constraint for AffExpr is LinearConstraint which is defined by an AffExpr and a lower and upper bound. If a solver interface does not support range constraints, this will automatically translated into two constraints at solve time. Constructing constraints manually is not an expected behavior and won't add the constraint to a model automatically. See below for the correct methods.There is also QuadExpr for quadratic expressions type that also provides a default constructor that takes no arguments and a full constructor. There are four fields: two vectors of variables, a vector of coefficients, and the affine part of the expression. This is best explained by example:aff = AffExpr([x, z], [3.0, 4.0], 2.0)  # 3x + 4z + 2\nquad = QuadExpr([x,y],[x,z],[3.0,4.0],aff)  # 3x^2 + 4yz + 3x + 4z + 2The corresponding constraint is QuadConstraint, which is expected to be a convex quadratic constraint."
},

{
    "location": "refexpr.html#JuMP.addconstraint",
    "page": "Expressions and Constraints",
    "title": "JuMP.addconstraint",
    "category": "Function",
    "text": "addconstraint(m::Model, c::LinearConstraint)\n\nAdd a linear constraint to Model m.\n\n\n\naddconstraint(m::Model, c::QuadConstraint)\n\nAdd a quadratic constraint to Model m.\n\n\n\naddconstraint(m::Model, c::SOCConstraint)\n\nAdd a SOC constraint to Model m.\n\n\n\naddconstraint(m::Model, c::SDConstraint)\n\nAdd a SD constraint to Model m.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.@constraint",
    "page": "Expressions and Constraints",
    "title": "JuMP.@constraint",
    "category": "Macro",
    "text": "@constraint(m::Model, con)\n\nadd linear or quadratic constraints.\n\n@constraint(m::Model, ref, con)\n\nadd groups of linear or quadratic constraints.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.@constraints",
    "page": "Expressions and Constraints",
    "title": "JuMP.@constraints",
    "category": "Macro",
    "text": "@constraints(m, args...)\n\nadds groups of constraints at once, in the same fashion as @constraint. The model must be the first argument, and multiple constraints can be added on multiple lines wrapped in a begin ... end block. For example:\n\n@constraints(m, begin\n  x >= 1\n  y - w <= 2\n  sum_to_one[i=1:3], z[i] + y == 1\nend)\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.@expression",
    "page": "Expressions and Constraints",
    "title": "JuMP.@expression",
    "category": "Macro",
    "text": "@expression(args...)\n\nefficiently builds a linear, quadratic, or second-order cone expression but does not add to model immediately. Instead, returns the expression which can then be inserted in other constraints. For example:\n\n@expression(m, shared, sum(i*x[i] for i=1:5))\n@constraint(m, shared + y >= 5)\n@constraint(m, shared + z <= 10)\n\nThe ref accepts index sets in the same way as @variable, and those indices can be used in the construction of the expressions:\n\n@expression(m, expr[i=1:3], i*sum(x[j] for j=1:3))\n\nAnonymous syntax is also supported:\n\nexpr = @expression(m, [i=1:3], i*sum(x[j] for j=1:3))\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.@SDconstraint",
    "page": "Expressions and Constraints",
    "title": "JuMP.@SDconstraint",
    "category": "Macro",
    "text": "@SDconstraint(m, x)\n\nAdds a semidefinite constraint to the Model m. The expression x must be a square, two-dimensional array.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.addSOS1",
    "page": "Expressions and Constraints",
    "title": "JuMP.addSOS1",
    "category": "Function",
    "text": "addSOS1(m::Model, coll::Vector{AffExpr})\n\nAdds special ordered set constraint of type 1 (SOS1). Specify the set as a vector of weighted variables, e.g. coll = [3x, y, 2z]. Note that solvers expect the weights to be unique. See here for more details. If there is no inherent weighting in your model, an SOS constraint is probably unnecessary.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.addSOS2",
    "page": "Expressions and Constraints",
    "title": "JuMP.addSOS2",
    "category": "Function",
    "text": "addSOS2(m::Model, coll::Vector{AffExpr})\n\nAdds special ordered set constraint of type 2 (SOS2). Specify the set as a vector of weighted variables, e.g. coll = [3x, y, 2z]. Note that solvers expect the weights to be unique. See here for more details.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.@LinearConstraint",
    "page": "Expressions and Constraints",
    "title": "JuMP.@LinearConstraint",
    "category": "Macro",
    "text": "@LinearConstraint(x)\n\nConstructs a LinearConstraint instance efficiently by parsing the x. The same as @constraint, except it does not attach the constraint to any model.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.@LinearConstraints",
    "page": "Expressions and Constraints",
    "title": "JuMP.@LinearConstraints",
    "category": "Macro",
    "text": "@LinearConstraints(m, args...)\n\nConstructs a vector of LinearConstraint objects. Similar to @LinearConstraint, except it accepts multiple constraints as input as long as they are separated by newlines.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.@QuadConstraint",
    "page": "Expressions and Constraints",
    "title": "JuMP.@QuadConstraint",
    "category": "Macro",
    "text": "@QuadConstraint(x)\n\nConstructs a QuadConstraint instance efficiently by parsing the x. The same as @constraint, except it does not attach the constraint to any model.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.@QuadConstraints",
    "page": "Expressions and Constraints",
    "title": "JuMP.@QuadConstraints",
    "category": "Macro",
    "text": "@QuadConstraints(m, args...)\n\nConstructs a vector of QuadConstraint objects. Similar to @QuadConstraint, except it accepts multiple constraints as input as long as they are separated by newlines.\n\n\n\n"
},

{
    "location": "refexpr.html#Base.push!-Union{Tuple{C}, Tuple{JuMP.GenericAffExpr{C,V},C,V}, Tuple{V}} where V where C",
    "page": "Expressions and Constraints",
    "title": "Base.push!",
    "category": "Method",
    "text": "Base.push!{C,V}(aff::GenericAffExpr{C,V}, new_coeff::C, new_var::V)\n\nAn efficient way to grow an affine expression by one term. For example, to add 5x to an existing expression aff, use push!(aff, 5.0, x). This is significantly more efficient than aff += 5.0*x.\n\n\n\n"
},

{
    "location": "refexpr.html#Base.append!-Union{Tuple{C}, Tuple{JuMP.GenericAffExpr{C,V},JuMP.GenericAffExpr{C,V}}, Tuple{V}} where V where C",
    "page": "Expressions and Constraints",
    "title": "Base.append!",
    "category": "Method",
    "text": "Base.append!{C,V}(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V})\n\nEfficiently append the terms of an affine expression to an existing affine expression. For example, given aff = 5.0*x and other = 7.0*y + 3.0*z, we can grow aff using append!(aff, other) which results in aff equaling 5x + 7y + 3z. This is significantly more efficient than using aff += other.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.linearterms",
    "page": "Expressions and Constraints",
    "title": "JuMP.linearterms",
    "category": "Function",
    "text": "linearterms(aff::GenericAffExpr)\n\nProvides an iterator over the (a_i::C,x_i::V) terms in affine expression sum_i a_i x_i+b.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.getvalue-Tuple{JuMP.GenericAffExpr{Float64,JuMP.Variable}}",
    "page": "Expressions and Constraints",
    "title": "JuMP.getvalue",
    "category": "Method",
    "text": "getvalue(a::AffExpr)\n\nEvaluate an AffExpr given the current solution values.\n\n\n\n"
},

{
    "location": "refexpr.html#JuMP.getvalue-Tuple{JuMP.GenericQuadExpr{Float64,JuMP.Variable}}",
    "page": "Expressions and Constraints",
    "title": "JuMP.getvalue",
    "category": "Method",
    "text": "getvalue(a::QuadExpr)\n\nEvaluate a QuadExpr given the current solution values.\n\n\n\n"
},

{
    "location": "refexpr.html#Methods-1",
    "page": "Expressions and Constraints",
    "title": "Methods",
    "category": "section",
    "text": "addconstraint\n@constraint\n@constraints\n@expression\n@SDconstraint\naddSOS1\naddSOS2\n@LinearConstraint\n@LinearConstraints\n@QuadConstraint\n@QuadConstraints\nBase.push!{C,V}(aff::GenericAffExpr{C,V}, new_coeff::C, new_var::V)\nBase.append!{C,V}(aff::GenericAffExpr{C,V}, other::GenericAffExpr{C,V})\nlinearterms\ngetvalue(a::AffExpr)\ngetvalue(a::QuadExpr)"
},

{
    "location": "refexpr.html#Constraint-References-1",
    "page": "Expressions and Constraints",
    "title": "Constraint References",
    "category": "section",
    "text": "In order to manipulate constraints after creation, it is necessary to maintain a reference. The simplest way to do this is to use the special three-argument named constraint syntax for @constraint, which additionally allows you to create groups of constraints indexed by sets analogously to @variable. For example:@variable(m, x[1:3])\n@variable(m, y[2:2:6])\n@constraint(m, xyconstr[i=1:3,j=6:-2:2], x[i] - y[j] == 1)adds 9 constraints to the model m of the expected form. The variable xyconstr is a collection of ConstraintRef{Model,LinearConstraint} instances indexed by the ranges 1:3 and 6:-2:2 (the ordered tuple (6,4,2)), so, for example xyconstr[2,4] is a reference to the constraint x[2] - y[4] == 1. Indices can have dependencies on preceding indices, e.g. triangular indexing is allowed:@constraint(m, triconstr[i=1:3,j=2i:2:6], x[i] - y[j] == 1)A condition can be added following the indices; a semicolon is used to separate index sets from the condition:@constraint(m, constr[i=1:5,j=1:5; i+j >= 3], x[i] - y[j] == 1)Note that only one condition can be added, although expressions can be built up by using the usual && and || logical operators.Anonymous syntax is supported:constr = @constraint(m, [i=1:5,j=1:5; i+j >= 3], x[i] - y[j] == 1)To obtain the dual of a constraint, call getdual on the constraint reference:println(getdual(xyconstr[1,6]))When an LP model is infeasible, getdual will return the corresponding component of the infeasibility ray (Farkas proof), if available from the solver.Dual information is also accessible for second-order cone problems as described below. Duals are unavailable for MIPs.One may retrieve the corresponding internal LinearConstraint object from a ConstraintRef{Model,LinearConstraint} object constr by calling LinearConstraint(constr). This functionality is not yet implemented for other classes of constraints.For users who prefer to generate constraints in an explicit loop, we also provide the @constraintref convenience macro, e.g.:@constraintref constraintName[1:3]You can then iterate over constraints and store references in this structure, e.g.:@variable(m, x[1:5] >= 0)\n@constraintref myCons[1:5]\nfor i = 1:5\n  myCons[i] = @constraint(m, x[i] >= i)\nend"
},

{
    "location": "refexpr.html#Conic-constraint-duals-1",
    "page": "Expressions and Constraints",
    "title": "Conic constraint duals",
    "category": "section",
    "text": "JuMP supports accessing the dual solutions to second-order cone problems. Dual multipliers on variable bounds, linear constraints, and second-order cone constraints are accessible through getdual() given the corresponding variable or constraint reference object. For second-order cone constraints, getdual(c::ConstraintRef{Model,SOCConstraint}) returns a vector of dual variables in the dimension of the corresponding cone. Duals are defined such that they are consistent in sign with linear programming duals in the case that the second-order cone constraints are inactive.For example:m = Model()\n@variable(m, x[1:2] >= 1)\n@variable(m, t)\n@objective(m, Min, t)\n@constraint(m, soc, norm( x[i] for i=1:2 ) <= t)\nstatus = solve(m)\n\n@show getvalue(x) # [1.000000000323643,1.0000000003235763]\n@show getvalue(t) # 1.4142135583106126\n@show getdual(x)  # [0.7071067807797846,0.7071067802906756]\n@show getdual(soc)# [-1.0000000004665652,0.707106779497123,0.707106779008014]Note that the negative of the dual vector getdual(soc) belongs to the second-order cone. See the MathProgBase documentation for more on the definition of the dual problem. The dual solutions returned by JuMP agree with the definitions from MathProgBase up to a possible change in sign."
},

{
    "location": "refexpr.html#Vectorized-operations-1",
    "page": "Expressions and Constraints",
    "title": "Vectorized operations",
    "category": "section",
    "text": "JuMP supports vectorized expressions and constraints for linear and quadratic models. Although this syntax may be familiar for users coming from MATLAB-based modeling languages, we caution that this syntax may be slower than the scalar versions using loops–-especially for large operations. Nevertheless, the syntax often proves useful, for example in constraints involving small, dense matrix-vector products.Linear algebraic operators are available to give meaning to expressions like A*x where A is a matrix of numbers and x is a vector of Variable objects. You may also use objects of type Array{Variable} in these kinds of expressions; for example, any object you construct with @variable where each of the index sets are of the form 1:n. For example:@variable(m, x[1:3,1:4])\nexpr = rand(3,3)*xis allowed, while:@variable(m, x[2:4])\nexpr = rand(3,3)*xis not. Addition and subtraction are also defined in similar ways, following the usual Julia rules for linear algebra over arrays.Vectorized constraints can be added to the model, using the elementwise comparison operators .==, .>=, and .<=. For instance, you can write constraints of the form:@variable(m, x[1:10])\nA = rand(5,10)\nb = rand(5)\n@constraint(m, A*x + b .<= 1)Note that scalar literals (such as 1 or 0) are allowed in expressions.Concatenation is also possible for these arrays of variables or expressions. For instance, the following will create a matrix of QuadExpr that you can use elsewhere in your model:@variable(m, x[1:3])\nA = [1 x'\n     x x*x']Finally, note that this feature is not currently supported directly in nonlinear expressions; for example, a matrix–vector product will not work inside a call to the @NLconstraint macro."
},

{
    "location": "probmod.html#",
    "page": "Problem Modification",
    "title": "Problem Modification",
    "category": "page",
    "text": ""
},

{
    "location": "probmod.html#Problem-Modification-1",
    "page": "Problem Modification",
    "title": "Problem Modification",
    "category": "section",
    "text": "It can be useful to modify models after they have been created and solved, for example when we are solving many similar models in succession or generating the model dynamically (e.g. column generation). Additionally it is sometimes desirable for the solver to re-start from the last solution to reduce running times for successive solves (\"hot-start\"). Where available, JuMP exposes this functionality."
},

{
    "location": "probmod.html#Differences-in-Solvers-1",
    "page": "Problem Modification",
    "title": "Differences in Solvers",
    "category": "section",
    "text": "Some solvers do not expose the ability to modify a model after creation - the model must be constructed from scratch each time. JuMP will use the ability to modify problems exposed by the solver if possible, and will still work even if the solver does not support this functionality by passing the complete problem to the solver every time."
},

{
    "location": "probmod.html#Modifying-variables-1",
    "page": "Problem Modification",
    "title": "Modifying variables",
    "category": "section",
    "text": "As before, variables can be added using the @variable macro. To remove a variable, one can set the bounds on that variable to zero, e.g.:setlowerbound(x, 0.0)\nsetupperbound(x, 0.0)While bound updates are applied immediately in JuMP, variable bound changes are not transmitted to the solver until solve is called again.To add variables that appear in existing constraints, e.g. in column generation, there is an alternative form of the @variable macro:@variable(m, x, objective = objcoef, inconstraints = constrrefs, coefficients = values)\n@variable(m, x >= lb, objective = objcoef, inconstraints = constrrefs, coefficients = values)\n@variable(m, x <= ub, objective = objcoef, inconstraints = constrrefs, coefficients = values)\n@variable(m, lb <= x <= ub, objective = objcoef, inconstraints = constrrefs, coefficients = values)\n@variable(m, lb <= x <= ub, Int, objective = objcoef, inconstraints = constrrefs, coefficients = values)  # Types are supportedwhere objcoef is the coefficient of the variable in the new problem, constrrefs is a vector of ConstraintRef, and values is a vector of numbers. To give an example, consider the following code snippet:m = Model()\n@variable(m, 0 <= x <= 1)\n@variable(m, 0 <= y <= 1)\n@objective(m, Max, 5x + 1y)\n@constraint(m, con, x + y <= 1)\nsolve(m)  # x = 1, y = 0\n@variable(m, 0 <= z <= 1, objective = 10.0, inconstraints = [con], coefficients = [1.0])\n# The constraint is now x + y + z <= 1\n# The objective is now 5x + 1y + 10z\nsolve(m)  # z = 1In some situations you may be adding all variables in this way. To do so, first define a set of empty constraints, e.g. :m = Model()\n@constraint(m, con, 0 <= 1)\n@objective(m, Max, 0)\n@variable(m, 0 <= x <= 1, objective = 5, inconstraints = [con], coefficients = [1.0])\n@variable(m, 0 <= y <= 1, objective = 1, inconstraints = [con], coefficients = [1.0])\n@variable(m, 0 <= z <= 1, objective = 10, inconstraints = [con], coefficients = [1.0])\nsolve(m)"
},

{
    "location": "probmod.html#Modifying-constraints-1",
    "page": "Problem Modification",
    "title": "Modifying constraints",
    "category": "section",
    "text": "JuMP does not currently support changing constraint coefficients. For less-than and greater-than constraints, the right-hand-side can be changed, e.g.:@constraint(m, mycon, x + y <= 4)\nsolve(m)\nJuMP.setRHS(mycon, 3)  # Now x + y <= 3\nsolve(m)  # Hot-start for LPs"
},

{
    "location": "probmod.html#Modifying-the-objective-1",
    "page": "Problem Modification",
    "title": "Modifying the objective",
    "category": "section",
    "text": "To change the objective, simply call @objective again - the previous objective function and sense will be replaced."
},

{
    "location": "probmod.html#Modifying-nonlinear-models-1",
    "page": "Problem Modification",
    "title": "Modifying nonlinear models",
    "category": "section",
    "text": "See nonlinear parameters Nonlinear Parameters."
},

{
    "location": "callbacks.html#",
    "page": "Solver Callbacks",
    "title": "Solver Callbacks",
    "category": "page",
    "text": ""
},

{
    "location": "callbacks.html#Solver-Callbacks-1",
    "page": "Solver Callbacks",
    "title": "Solver Callbacks",
    "category": "section",
    "text": "Many mixed-integer (linear, conic, and nonlinear) programming solvers offer the ability to modify the solve process. Examples include changing branching decisions in branch-and-bound, adding custom cutting planes, providing custom heuristics to find feasible solutions, or implementing on-demand separators to add new constraints only when they are violated by the current solution (also known as lazy constraints).While historically this functionality has been limited to solver-specific interfaces, JuMP provides solver-independent support for a number of commonly used solver callbacks. Currently, we support lazy constraints, user-provided cuts, and user-provided heuristics for the Gurobi, CPLEX, GLPK, and SCIP solvers. We do not yet support any other class of callbacks, but they may be accessible by using the solver's low-level interface."
},

{
    "location": "callbacks.html#Lazy-Constraints-1",
    "page": "Solver Callbacks",
    "title": "Lazy Constraints",
    "category": "section",
    "text": "Lazy constraints are useful when the full set of constraints is too large to explicitly include in the initial formulation. When a MIP solver reaches a new solution, for example with a heuristic or by solving a problem at a node in the branch-and-bound tree, it will give the user the chance to provide constraint(s) that would make the current solution infeasible. For some more information about lazy constraints, see this blog post by Paul Rubin.There are three important steps to providing a lazy constraint callback. First we must write a function that will analyze the current solution that takes a single argument, e.g. function myLazyConGenerator(cb), where cb is a reference to the callback management code inside JuMP. Next you will do whatever analysis of the solution you need to inside your function to generate the new constraint before adding it to the model with @lazyconstraint(cb, myconstraint). There is an optional keyword option localcut to @lazyconstraint, which indicates if the lazy constraint that will be added will only apply at the current node and the tree rooted at that node. For example, @lazyconstraint(cb, myconstraint, localcut=true). By default, localcut is set to false. Finally we notify JuMP that this function should be used for lazy constraint generation using the addlazycallback(m, myLazyConGenerator) function before we call solve(m).The following is a simple example to make this more clear. In this two-dimensional problem we have a set of box constraints explicitly provided and a set of two lazy constraints we can add on the fly. The solution without the lazy constraints will be either (0,2) or (2,2), and the final solution will be (1,2):using JuMP\nusing Gurobi\n\n# We will use Gurobi\nm = Model(solver=GurobiSolver())\n\n# Define our variables to be inside a box, and integer\n@variable(m, 0 <= x <= 2, Int)\n@variable(m, 0 <= y <= 2, Int)\n\n@objective(m, Max, y)\n\n# We now define our callback function that takes one argument,\n# the callback handle. Note that we can access m, x, and y because\n# this function is defined inside the same scope\nfunction corners(cb)\n    x_val = getvalue(x)\n    y_val = getvalue(y)\n    println(\"In callback function, x=$x_val, y=$y_val\")\n\n    # We have two constraints, one cutting off the top\n    # left corner and one cutting off the top right corner, e.g.\n    # (0,2) +---+---+ (2,2)\n    #       |xx/ \\xx|\n    #       |x/   \\x|\n    #       |/     \\|\n    #       +       +\n    #       |       |\n    #       |       |\n    #       |       |\n    # (0,0) +---+---+ (2,0)\n\n    # Allow for some impreciseness in the solution\n    TOL = 1e-6\n\n    # Check top left, allowing some tolerance\n    if y_val - x_val > 1 + TOL\n        # Cut off this solution\n        println(\"Solution was in top left, cut it off\")\n        # Use the original variables, but not m - cb instead\n        @lazyconstraint(cb, y - x <= 1)\n    # Check top right\n    elseif y_val + x_val > 3 + TOL\n        # Cut off this solution\n        println(\"Solution was in top right, cut it off\")\n        # Use the original variables, but not m - cb instead\n        @lazyconstraint(cb, y + x <= 3)\n    end\nend  # End of callback function\n\n# Tell JuMP/Gurobi to use our callback function\naddlazycallback(m, corners)\n\n# Solve the problem\nsolve(m)\n\n# Print our final solution\nprintln(\"Final solution: [ $(getvalue(x)), $(getvalue(y)) ]\")The code should print something like (amongst the output from Gurobi):In callback function, x=2.0, y=2.0\nSolution was in top right, cut it off\nIn callback function, x=0.0, y=2.0\nSolution was in top left, cut it off\nIn callback function, x=1.0, y=2.0\nFinal solution: [ 1.0, 2.0 ]This code can also be found in /JuMP/examples/simplelazy.jl.There is an optional fractional keyword option to addlazycallback which indicates that the callback may be called at solutions that do not satisfy integrality constraints. For example, addlazycallback(m, myLazyConGenerator, fractional=true). Depending on the solver, this may invoke the callback after solving each LP relaxation in the Branch and Bound tree. By default, fractional is set to false."
},

{
    "location": "callbacks.html#User-Cuts-1",
    "page": "Solver Callbacks",
    "title": "User Cuts",
    "category": "section",
    "text": "User cuts, or simply cuts, provide a way for the user to tighten the LP relaxation using problem-specific knowledge that the solver cannot or is unable to infer from the model. Just like with lazy constraints, when a MIP solver reaches a new node in the branch-and-bound tree, it will give the user the chance to provide cuts to make the current relaxed (fractional) solution infeasible in the hopes of obtaining an integer solution. For more details about the difference between user cuts and lazy constraints see the aforementioned blog post.Your user cuts should not change the set of integer feasible solutions. Equivalently, your cuts can only remove fractional solutions - that is, \"tighten\" the LP relaxation of the MILP. If you add a cut that removes an integer solution, the solver may return an incorrect solution.Adding a user cut callback is similar to adding a lazy constraint callback. First we must write a function that will analyze the current solution that takes a single argument, e.g. function myUserCutGenerator(cb), where cb is a reference to the callback management code inside JuMP. Next you will do whatever analysis of the solution you need to inside your function to generate the new constraint before adding it to the model with the JuMP macro @usercut(cb, myconstraint) (same limitations as addConstraint). There is an optional keyword option localcut to @usercut, which indicates if the user cut that will be added will only apply at the current node and the tree rooted at that node. For example, @usercut(cb, myconstraint, localcut=true). By default, localcut is set to false. Finally we notify JuMP that this function should be used for user cut generation using the addcutcallback(m, myUserCutGenerator) function before we call solve(m).Consider the following example which is related to the lazy constraint example. The problem is two-dimensional, and the objective sense prefers solution in the top-right of a 2-by-2 square. There is a single constraint that cuts off the top-right corner to make the LP relaxation solution fractional. We will exploit our knowledge of the problem structure to add a user cut that will make the LP relaxation integer, and thus solve the problem at the root node:using JuMP\nusing Gurobi\n\n# We will use Gurobi, which requires that we manually set the attribute\n# PreCrush to 1 if we have user cuts. We will also disable PreSolve, Cuts,\n# and Heuristics so only our cut will be used\nm = Model(solver=GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0))\n\n# Define our variables to be inside a box, and integer\n@variable(m, 0 <= x <= 2, Int)\n@variable(m, 0 <= y <= 2, Int)\n\n# Optimal solution is trying to go towards top-right corner (2.0, 2.0)\n@objective(m, Max, x + 2y)\n\n# We have one constraint that cuts off the top right corner\n@constraint(m, y + x <= 3.5)\n\n# Optimal solution of relaxed problem will be (1.5, 2.0)\n# We can add a user cut that will cut of this fractional solution.\n\n# We now define our callback function that takes one argument,\n# the callback handle. Note that we can access m, x, and y because\n# this function is defined inside the same scope\nfunction mycutgenerator(cb)\n    x_val = getvalue(x)\n    y_val = getvalue(y)\n    println(\"In callback function, x=$x_val, y=$y_val\")\n\n    # Allow for some impreciseness in the solution\n    TOL = 1e-6\n\n    # Check top right\n    if y_val + x_val > 3 + TOL\n        # Cut off this solution\n        println(\"Fractional solution was in top right, cut it off\")\n        # Use the original variables\n        @usercut(cb, y + x <= 3)\n    end\nend  # End of callback function\n\n# Tell JuMP/Gurobi to use our callback function\naddcutcallback(m, mycutgenerator)\n\n# Solve the problem\nsolve(m)\n\n# Print our final solution\nprintln(\"Final solution: [ $(getvalue(x)), $(getvalue(y)) ]\")The code should print something like (amongst the output from Gurobi):In callback function, x=1.5, y=2.0\nFractional solution was in top right, cut it off\nIn callback function, x=1.0, y=2.0\nFinal solution: [ 1.0, 2.0 ]This code can also be found in /JuMP/examples/simpleusercut.jl."
},

{
    "location": "callbacks.html#User-Heuristics-1",
    "page": "Solver Callbacks",
    "title": "User Heuristics",
    "category": "section",
    "text": "Integer programming solvers frequently include heuristics that run at the nodes of the branch-and-bound tree. They aim to find integer solutions quicker than plain branch-and-bound would to tighten the bound, allowing us to fathom nodes quicker and to tighten the integrality gap. Some heuristics take integer solutions and explore their \"local neighborhood\" (e.g. flipping binary variables, fix some variables and solve a smaller MILP, ...) and others take fractional solutions and attempt to round them in an intelligent way. You may want to add a heuristic of your own if you have some special insight into the problem structure that the solver is not aware of, e.g. you can consistently take fractional solutions and intelligently guess integer solutions from them.The user heuristic callback is somewhat different from the previous two heuristics. The general concept is that we can create multiple partial solutions and submit them back to the solver - each solution must be submitted before a new solution is constructed. As before we provide a function that analyzes the current solution and takes a single argument, e.g. function myHeuristic(cb), where cb is a reference to the callback management code inside JuMP. You can build your solutions using setsolutionvalue(cb, x, value) and submit them with addsolution(cb). Note that addsolution will \"wipe\" the previous (partial) solution. Notify JuMP that this function should be used as a heuristic using the addheuristiccallback(m, myHeuristic) function before calling solve(m).There is some unavoidable (for performance reasons) solver-dependent behavior - you should check your solver documentation for details. For example: GLPK will not check the feasibility of your heuristic solution. If you need to submit many heuristic solutions in one callback, there may be performance impacts from the \"wiping\" behavior of addsolution - please file an issue and we can address this issue.Consider the following example, which is the same problem as seen in the user cuts section. The heuristic simply rounds the fractional variable to generate integer solutions.:using JuMP\nusing Gurobi\n\n# We will use Gurobi and disable PreSolve, Cuts, and (in-built) Heuristics so\n# only our heuristic will be used\nm = Model(solver=GurobiSolver(Cuts=0, Presolve=0, Heuristics=0.0))\n\n# Define our variables to be inside a box, and integer\n@variable(m, 0 <= x <= 2, Int)\n@variable(m, 0 <= y <= 2, Int)\n\n# Optimal solution is trying to go towards top-right corner (2.0, 2.0)\n@objective(m, Max, x + 2y)\n\n# We have one constraint that cuts off the top right corner\n@constraint(m, y + x <= 3.5)\n\n# Optimal solution of relaxed problem will be (1.5, 2.0)\n\n# We now define our callback function that takes one argument,\n# the callback handle. Note that we can access m, x, and y because\n# this function is defined inside the same scope\nfunction myheuristic(cb)\n    x_val = getvalue(x)\n    y_val = getvalue(y)\n    println(\"In callback function, x=$x_val, y=$y_val\")\n\n    setsolutionvalue(cb, x, floor(x_val))\n    # Leave y undefined - solver should handle as it sees fit. In the case\n    # of Gurobi it will try to figure out what it should be.\n    addsolution(cb)\n\n    # Submit a second solution\n    setsolutionvalue(cb, x, ceil(x_val))\n    addsolution(cb)\nend  # End of callback function\n\n# Tell JuMP/Gurobi to use our callback function\naddheuristiccallback(m, myheuristic)\n\n# Solve the problem\nsolve(m)\n\n# Print our final solution\nprintln(\"Final solution: [ $(getvalue(x)), $(getvalue(y)) ]\")The code should print something like:In callback function, x=1.5, y=2.0\n     0     0    5.50000    0    1          -    5.50000     -      -    0s\nH    1     0                       5.0000000    5.50000  10.0%   0.0    0swhere the H denotes a solution found with a heuristic - our heuristic in this case. This code can also be found in /JuMP/examples/simpleheur.jl."
},

{
    "location": "callbacks.html#Querying-Solver-Progress-1",
    "page": "Solver Callbacks",
    "title": "Querying Solver Progress",
    "category": "section",
    "text": "All JuMP callback methods must take a single argument, called cb by convention. cb is a handle to the internal callback system used by the underlying solver, and allows the user to query solver state. There are a variety of methods available which are listed in the MathProgBase documentation including:cbgetobj(cb)\ncbgetbestbound(cb)\ncbgetexplorednodes(cb)\ncbgetstate(cb)"
},

{
    "location": "callbacks.html#Informational-Callbacks-1",
    "page": "Solver Callbacks",
    "title": "Informational Callbacks",
    "category": "section",
    "text": "Sometimes it can be useful to track solver progress without actually changing the algorithm by adding cuts or heuristic solutions. In these cases, informational callbacks can be added, wherein statistics can be tracked via the cbget functions discussed in the previous section. Informational callbacks are added to a JuMP model with the addinfocallback(m::Model, f::Function; when::Symbol) function, where the when argument should be one of :MIPNode, :MIPSol or :Intermediate (listed under cbgetstate() in the MathProgBase documentation)For a simple example, we can add a function that tracks the best bound and incumbent objective value as the solver progresses through the branch-and-bound tree:type NodeData\n    time::Float64  # in seconds since the epoch\n    node::Int\n    obj::Float64\n    bestbound::Float64\nend\n\n# build model ``m`` up here\n\nbbdata = NodeData[]\n\nfunction infocallback(cb)\n    node      = MathProgBase.cbgetexplorednodes(cb)\n    obj       = MathProgBase.cbgetobj(cb)\n    bestbound = MathProgBase.cbgetbestbound(cb)\n    push!(bbdata, NodeData(time(),node,obj,bestbound))\nend\naddinfocallback(m, infocallback, when = :Intermediate)\n\nsolve(m)\n\n# Save results to file for analysis later\nopen(\"bbtrack.csv\",\"w\") do fp\n    println(fp, \"time,node,obj,bestbound\")\n    for bb in bbdata\n        println(fp, bb.time, \",\", bb.node, \",\",\n                    bb.obj, \",\", bb.bestbound)\n    end\nendFor a second example, we can add a function that tracks the intermediate solutions at each integer-feasible solution in the Branch-and-Bound tree:solutionvalues = Vector{Float64}[]\n\n# build model ``m`` up here\n\nfunction infocallback(cb)\n    push!(solutionvalues, JuMP.getvalue(x))\nend\naddinfocallback(m, infocallback, when = :MIPSol)\n\nsolve(m)\n\n# all the intermediate solutions are now stored in `solutionvalues`"
},

{
    "location": "callbacks.html#Code-Design-Considerations-1",
    "page": "Solver Callbacks",
    "title": "Code Design Considerations",
    "category": "section",
    "text": "In the above examples the callback function is defined in the same scope as the model and variable definitions, allowing us to access them. If we defined the function in some other scope, or even file, we would not be able to access them directly. The proposed solution to this design problem is to separate the logic of analyzing the current solution values from the callback itself. This has many benefits, including writing unit tests for the callback function to check its correctness. The callback function passed to JuMP is then simply a stub that extracts the current solution and any other relevant information and passes that to the constraint generation logic. To apply this to our previous lazy constraint example, consider the following code:using JuMP\nusing Gurobi\nusing Base.Test\n\nfunction cornerChecker(x_val, y_val)\n    # This function does not depend on the model, and could\n    # be written anywhere. Instead, it returns a tuple of\n    # values (newcut, x_coeff, y_coeff, rhs) where newcut is a\n    # boolean if a cut is needed, x_coeff is the coefficient\n    # on the x variable, y_coeff is the coefficient on\n    # the y variable, and rhs is the right hand side\n    TOL = 1e-6\n    if y_val - x_val > 1 + TOL\n        return (true, -1.0, 1.0, 1.0)  # Top left\n    elseif y_val + x_val > 3 + TOL\n        return (true,  1.0, 1.0, 3.0)  # Top right\n    else\n        return (false, 0.0, 0.0, 0.0)  # No cut\n    end\nend\n\n# A unit test for the cornerChecker function\nfunction test_cornerChecker()\n    # Test the four corners - only two should produce cuts\n\n    newcut, x_coeff, y_coeff, rhs = cornerChecker(0, 0)\n    @test !newcut\n\n    newcut, x_coeff, y_coeff, rhs = cornerChecker(2, 0)\n    @test !newcut\n\n    newcut, x_coeff, y_coeff, rhs = cornerChecker(0, 2)\n    @test newcut\n    @test x_coeff == -1.0\n    @test y_coeff ==  1.0\n    @test rhs == 1.0\n\n    newcut, x_coeff, y_coeff, rhs = cornerChecker(2, 2)\n    @test newcut\n    @test x_coeff ==  1.0\n    @test y_coeff ==  1.0\n    @test rhs == 3.0\nend\n\nfunction solveProblem()\n    m = Model(solver=GurobiSolver())\n\n    @variable(m, 0 <= x <= 2, Int)\n    @variable(m, 0 <= y <= 2, Int)\n    @objective(m, Max, y)\n\n    # Note that the callback is now a stub that passes off\n    # the work to the \"algorithm\"\n    function corners(cb)\n        x_val = getvalue(x)\n        y_val = getvalue(y)\n        println(\"In callback function, x=$x_val, y=$y_val\")\n\n        newcut, x_coeff, y_coeff, rhs = cornerChecker(x_val, y_val)\n\n        if newcut\n            @lazyconstraint(cb, x_coeff*x + y_coeff*y <= rhs)\n        end\n    end  # End of callback function\n\n    addlazycallback(m, corners)\n    solve(m)\n    println(\"Final solution: [ $(getvalue(x)), $(getvalue(y)) ]\")\nend\n\n# Run tests\ntest_cornerChecker()\n\n# Solve it\nsolveProblem()This code can also be found in /JuMP/examples/simplelazy2.jl."
},

{
    "location": "callbacks.html#Exiting-a-callback-early-1",
    "page": "Solver Callbacks",
    "title": "Exiting a callback early",
    "category": "section",
    "text": "If you need to exit the optimization process earlier than a solver otherwise would, it is possible to return JuMP.StopTheSolver from the callback code:return JuMP.StopTheSolverThis will trigger the solver to exit immediately and return a :UserLimit status."
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
    "text": "JuMP has support for general smooth nonlinear (convex and nonconvex) optimization problems. JuMP is able to provide exact, sparse second-order derivatives to solvers. This information can improve solver accuracy and performance.Nonlinear objectives and constraints are specified by using the @NLobjective and @NLconstraint macros. The familiar sum() syntax is supported within these macros, as well as prod() which analogously represents the product of the terms within. Note that the @objective and @constraint macros (and corresponding functions) do not currently support nonlinear expressions. However, a model can contain a mix of linear, quadratic, and nonlinear constraints or objective functions. Starting points may be provided by using the start keyword argument to @variable. If a starting value is not provided for a variable, it will be set to the projection of zero onto the interval defined by the variable bounds. For nonconvex problems, the returned solution is only guaranteed to be locally optimal. Convexity detection is not currently provided.For example, we can solve the classical Rosenbrock problem (with a twist) as follows:using JuMP\nm = Model()\n@variable(m, x, start = 0.0)\n@variable(m, y, start = 0.0)\n\n@NLobjective(m, Min, (1-x)^2 + 100(y-x^2)^2)\n\nsolve(m)\nprintln(\"x = \", getvalue(x), \" y = \", getvalue(y))\n\n# adding a (linear) constraint\n@constraint(m, x + y == 10)\nsolve(m)\nprintln(\"x = \", getvalue(x), \" y = \", getvalue(y))Examples: optimal control, maximum likelihood estimation, and Hock-Schittkowski tests."
},

{
    "location": "nlp.html#Syntax-notes-1",
    "page": "Nonlinear Modeling",
    "title": "Syntax notes",
    "category": "section",
    "text": "The syntax accepted in nonlinear expressions is more restricted than the syntax for linear and quadratic expressions. We note some important points below.All expressions must be simple scalar operations. You cannot use dot, matrix-vector products, vector slices, etc. Translate vector operations into explicit sum() operations or use the AffExpr plus auxiliary variable trick described below.\nThere is no operator overloading provided to build up nonlinear expressions. For example, if x is a JuMP variable, the code 3x will return an AffExpr object that can be used inside of future expressions and linear constraints. However, the code sin(x) is an error. All nonlinear expressions must be inside of macros.\nUser-defined Functions may be used within nonlinear expressions only after they are registered. For example:    myfunction(a,b) = exp(a)*b\n    @variable(m, x); @variable(m, y)\n    @NLobjective(m, Min, myfunction(x,y)) # ERROR. Needs JuMP.register() first.\n    @NLobjective(m, Min, exp(x)*y) # OkayAffExpr and QuadExpr objects cannot currently be used inside nonlinear expressions. Instead, introduce auxiliary variables, e.g.:    myexpr = dot(c,x) + 3y # where x and y are variables\n    @variable(m, aux)\n    @constraint(m, aux == myexpr)\n    @NLobjective(m, Min, sin(aux))You can declare embeddable nonlinear expressions with @NLexpression. For example:    @NLexpression(m, myexpr[i=1:n], sin(x[i]))\n    @NLconstraint(m, myconstr[i=1:n], myexpr[i] <= 0.5)Anonymous syntax is supported in @NLexpression and @NLconstraint:    myexpr = @NLexpression(m, [i=1:n], sin(x[i]))\n    myconstr = @NLconstraint(m, [i=1:n], myexpr[i] <= 0.5)"
},

{
    "location": "nlp.html#Nonlinear-Parameters-1",
    "page": "Nonlinear Modeling",
    "title": "Nonlinear Parameters",
    "category": "section",
    "text": "For nonlinear models only, JuMP offers a syntax for explicit \"parameter\" objects which can be used to modify a model in-place just by updating the value of the parameter. Nonlinear parameters are declared by using the @NLparameter macro and may be indexed by arbitrary sets analogously to JuMP variables and expressions. The initial value of the parameter must be provided on the right-hand side of the == sign as seen below:@NLparameter(m, x == 10)\n@NLparameter(m, y[i=1:10] == my_data[i]) # set of parameters indexed from 1 to 10You may use getvalue and setvalue to query or update the value of a parameter:getvalue(x) # 10, from above\nsetvalue(y[4], 54.3) # y[4] now holds the value 54.3Nonlinear parameters can be used within nonlinear expressions only:@variable(m, z)\n@objective(m, Max, x*z)       # error: x is a nonlinear parameter\n@NLobjective(m, Max, x*z)     # ok\n@expression(m, my_expr, x*z^2)      # error: x is a nonlinear parameter\n@NLexpression(m, my_nl_expr, x*z^2) # okNonlinear parameters are useful when solving nonlinear models in a sequence:m = Model()\n@variable(m, z)\n@NLparameter(m, x == 1.0)\n@NLobjective(m, Min, (z-x)^2)\nsolve(m)\ngetvalue(z) # equals 1.0\n\n# Now, update the value of x to solve a different problem\nsetvalue(x, 5.0)\nsolve(m)\ngetvalue(z) # equals 5.0Using nonlinear parameters can be faster than creating a new model from scratch with updated data because JuMP is able to avoid repeating a number of steps in processing the model before handing it off to the solver."
},

{
    "location": "nlp.html#User-defined-Functions-1",
    "page": "Nonlinear Modeling",
    "title": "User-defined Functions",
    "category": "section",
    "text": "JuMP's library of recognized univariate functions is derived from the Calculus.jl package. If you encounter a standard special function not currently supported by JuMP, consider contributing to the list of derivative rules there. In addition to this built-in list of functions, it is possible to register custom (user-defined) nonlinear functions to use within nonlinear expressions. JuMP does not support black-box optimization, so all user-defined functions must provide derivatives in some form. Fortunately, JuMP supports automatic differentiation of user-defined functions, a feature to our knowledge not available in any comparable modeling systems.Automatic differentiation is not finite differencing. JuMP's automatically computed derivatives are not subject to approximation error.JuMP uses ForwardDiff.jl to perform automatic differentiation; see the ForwardDiff.jl documentation for a description of how to write a function suitable for automatic differentiation. The general guideline is to write code that is generic with respect to the number type; don't assume that the input to the function is Float64. To register a user-defined function with derivatives computed by automatic differentiation, use the JuMP.register method as in the following example:mysquare(x) = x^2\nmyf(x,y) = (x-1)^2+(y-2)^2\n\nm = Model()\n\nJuMP.register(m, :myf, 2, myf, autodiff=true)\nJuMP.register(m, :mysquare, 1, mysquare, autodiff=true)\n\n@variable(m, x[1:2] >= 0.5)\n@NLobjective(m, Min, myf(x[1],mysquare(x[2])))The above code creates a JuMP model with the objective function (x[1]-1)^2 + (x[2]^2-2)^2. The first argument to JuMP.register the model for which the functions are registered. The second argument is a Julia symbol object which serves as the name of the user-defined function in JuMP expressions; the JuMP name need not be the same as the name of the corresponding Julia method. The third argument specifies how many arguments the function takes. The fourth argument is the name of the Julia method which computes the function, and autodiff=true instructs JuMP to compute exact gradients automatically.note: Note\nAll arguments to user-defined functions are scalars, not vectors. To define a function which takes a large number of arguments, you may use the splatting syntax f(x...) = ....Forward-mode automatic differentiation as implemented by ForwardDiff.jl has a computational cost that scales linearly with the number of input dimensions. As such, it is not the most efficient way to compute gradients of user-defined functions if the number of input arguments is large. In this case, users may want to provide their own routines for evaluating gradients. The more general syntax for JuMP.register which accepts user-provided derivative evaluation routines is:JuMP.register(m::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)The input differs for functions which take a single input argument and functions which take more than one. For univariate functions, the derivative evaluation routines should return a number which represents the first and second-order derivatives respectively. For multivariate functions, the derivative evaluation routines will be passed a gradient vector which they must explicitly fill. Second-order derivatives of multivariate functions are not currently supported; this argument should be omitted. The following example sets up the same optimization problem as before, but now we explicitly provide evaluation routines for the user-defined functions:mysquare(x) = x^2\nmysquare_prime(x) = 2x\nmysquare_primeprime(x) = 2\n\nmyf(x,y) = (x-1)^2+(y-2)^2\nfunction ∇f(g,x,y)\n    g[1] = 2*(x-1)\n    g[2] = 2*(y-2)\nend\n\nm = Model()\n\nJuMP.register(m, :myf, 2, myf, ∇f)\nJuMP.register(m, :mysquare, 1, mysquare, mysquare_prime, mysquare_primeprime)\n\n@variable(m, x[1:2] >= 0.5)\n@NLobjective(m, Min, myf(x[1],mysquare(x[2])))"
},

{
    "location": "nlp.html#Factors-affecting-solution-time-1",
    "page": "Nonlinear Modeling",
    "title": "Factors affecting solution time",
    "category": "section",
    "text": "The execution time when solving a nonlinear programming problem can be divided into two parts, the time spent in the optimization algorithm (the solver) and the time spent evaluating the nonlinear functions and corresponding derivatives. Ipopt explicitly displays these two timings in its output, for example:Total CPU secs in IPOPT (w/o function evaluations)   =      7.412\nTotal CPU secs in NLP function evaluations           =      2.083For Ipopt in particular, one can improve the performance by installing advanced sparse linear algebra packages, see Installation Guide. For other solvers, see their respective documentation for performance tips.The function evaluation time, on the other hand, is the responsibility of the modeling language. JuMP computes derivatives by using the ReverseDiffSparse package, which implements, in pure Julia, reverse-mode automatic differentiation with graph coloring methods for exploiting sparsity of the Hessian matrix [1]. As a conservative bound, JuMP's performance here currently may be expected to be within a factor of 5 of AMPL's."
},

{
    "location": "nlp.html#Querying-derivatives-from-a-JuMP-model-1",
    "page": "Nonlinear Modeling",
    "title": "Querying derivatives from a JuMP model",
    "category": "section",
    "text": "For some advanced use cases, one may want to directly query the derivatives of a JuMP model instead of handing the problem off to a solver. Internally, JuMP implements the AbstractNLPEvaluator interface from MathProgBase. To obtain an NLP evaluator object from a JuMP model, use JuMP.NLPEvaluator. The linearindex method maps from JuMP variables to the variable indices at the MathProgBase level.For example:m = Model()\n@variable(m, x)\n@variable(m, y)\n\n@NLobjective(m, Min, sin(x) + sin(y))\nvalues = zeros(2)\nvalues[linearindex(x)] = 2.0\nvalues[linearindex(y)] = 3.0\n\nd = JuMP.NLPEvaluator(m)\nMathProgBase.initialize(d, [:Grad])\nobjval = MathProgBase.eval_f(d, values) # == sin(2.0) + sin(3.0)\n\n∇f = zeros(2)\nMathProgBase.eval_grad_f(d, ∇f, values)\n# ∇f[linearindex(x)] == cos(2.0)\n# ∇f[linearindex(y)] == cos(3.0)The ordering of constraints in a JuMP model corresponds to the following ordering at the MathProgBase nonlinear abstraction layer. There are three groups of constraints: linear, quadratic, and nonlinear. Linear and quadratic constraints, to be recognized as such, must be added with the @constraint macros. All constraints added with the @NLconstraint macros are treated as nonlinear constraints. Linear constraints are ordered first, then quadratic, then nonlinear. The linearindex method applied to a constraint reference object returns the index of the constraint within its corresponding constraint class. For example:m = Model()\n@variable(m, x)\n@constraint(m, cons1, x^2 <= 1)\n@constraint(m, cons2, x + 1 == 3)\n@NLconstraint(m, cons3, x + 5 == 10)\n\ntypeof(cons1) # JuMP.ConstraintRef{JuMP.Model,JuMP.GenericQuadConstraint{JuMP.GenericQuadExpr{Float64,JuMP.Variable}}} indicates a quadratic constraint\ntypeof(cons2) # JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}} indicates a linear constraint\ntypeof(cons3) # JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.NonlinearExprData}} indicates a nonlinear constraint\nlinearindex(cons1) == linearindex(cons2) == linearindex(cons3) == 1When querying derivatives, cons2 will appear first, because it is the first linear constraint, then cons1, because it is the first quadratic constraint, then cons3, because it is the first nonlinear constraint. Note that for one-sided nonlinear constraints, JuMP subtracts any values on the right-hand side when computing expressions. In other words, one-sided nonlinear constraints are always transformed to have a right-hand side of zero.The JuMP.constraintbounds(m::Model) method returns the lower and upper bounds of all the constraints in the model, concatenated in the order discussed above.This method of querying derivatives directly from a JuMP model is convenient for interacting with the model in a structured way, e.g., for accessing derivatives of specific variables. For example, in statistical maximum likelihood estimation problems, one is often interested in the Hessian matrix at the optimal solution, which can be queried using the JuMP.NLPEvaluator.If you are writing a \"solver\", we highly encourage use of the MathProgBase nonlinear interface over querying derivatives using the above methods. These methods are provided for convenience but do not fully integrate with JuMP's solver infrastructure. In particular, they do not allow users to specify your solver to the Model() constructor nor to call it using solve() nor to populate the solution back into the model. Use of the MathProgBase interface also has the advantage of being independent of JuMP itself; users of MathProgBase solvers are free to implement their own evaluation routines instead of expressing their model in JuMP. You may use the JuMP.build method to ask JuMP to populate the \"solver\" without calling optimize!."
},

{
    "location": "nlp.html#Raw-expression-input-1",
    "page": "Nonlinear Modeling",
    "title": "Raw expression input",
    "category": "section",
    "text": "In addition to the @NLobjective and @NLconstraint macros, it is also possible to provide Julia Expr objects directly by using JuMP.setNLobjective and JuMP.addNLconstraint. This input form may be useful if the expressions are generated programmatically. JuMP variables should be spliced into the expression object. For example:@variable(m, 1 <= x[i=1:4] <= 5)\nJuMP.setNLobjective(m, :Min, :($(x[1])*$(x[4])*($(x[1])+$(x[2])+$(x[3])) + $(x[3])))\nJuMP.addNLconstraint(m, :($(x[1])*$(x[2])*$(x[3])*$(x[4]) >= 25))\n\n# Equivalent form using traditional JuMP macros:\n@NLobjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])\n@NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)See the Julia documentation for more examples and description of Julia expressions.[1]: Dunning, Huchette, and Lubin, \"JuMP: A Modeling Language for Mathematical Optimization\", arXiv."
},

]}
