var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#JuMP-1",
    "page": "Introduction",
    "title": "JuMP",
    "category": "section",
    "text": "# These comments do not display in the HTML output.\n# See https://github.com/JuliaDocs/Documenter.jl/issues/674.\n\n# Style conventions for the JuMP documentation:\n# - Respect the 80-character line limit whenever possible.\n# - Be concise.\n# - Use lists instead of long sentences.\n#   - Use numbered lists when describing a sequence, e.g., (1) do X, (2) then Y.\n#   - Use bullet points when the items are not ordered.\n# - Example code should be covered by doctests.\n#   - But it\'s unclear what to do if the code depends on a solver, see\n#     https://github.com/JuliaOpt/JuMP.jl/issues/1175.warning: Warning\nThis documentation is for the development version of JuMP. JuMP is undergoing a major transition to MathOptInterface, and the documentation is in the process of being rewritten. We recommend using the development version only for (1) developers of packages upstream or downstream of JuMP or (2) early adopters willing to provide feedback and file issues.JuMP is a domain-specific modeling language for mathematical optimization embedded in Julia. It currently supports a number of open-source and commercial solvers (see below) for a variety of problem classes, including linear programming, mixed-integer programming, second-order conic programming, semidefinite programming, and nonlinear programming. JuMP\'s features include:User friendliness\nSyntax that mimics natural mathematical expressions.\nComplete documentation (WIP!)\nSpeed\nBenchmarking has shown that JuMP can create problems at similar speeds   to special-purpose modeling languages such as   AMPL.\nJuMP communicates with most solvers in memory, avoiding the need to   write intermediary files.\nSolver independence\nJuMP uses a generic solver-independent interface provided by the   MathOptInterface   package, making it easy to change between a number of open-source and   commercial optimization software packages (\"solvers\").\nCurrently supported solvers include   Artelys Knitro,   Bonmin,   Cbc,   Clp,   Couenne,   CPLEX,   ECOS,   FICO Xpress,   GLPK,   Gurobi,   Ipopt,   MOSEK,   NLopt, and   SCS.\nAccess to advanced algorithmic techniques\nIncluding efficient LP re-solves which previously required using   solver-specific and/or low-level C++ libraries.\nEase of embedding\nJuMP itself is written purely in Julia. Solvers are the only binary   dependencies.\nBeing embedded in a general-purpose programming language makes it easy   to solve optimization problems as part of a larger workflow (e.g.,   inside a simulation, behind a web server, or as a subproblem in a   decomposition algorithm).\nAs a trade-off, JuMP\'s syntax is constrained by the syntax available   in Julia.\nJuMP is MPL licensed, meaning that   it can be embedded in commercial software that complies with the terms   of the license.While neither Julia nor JuMP have reached version 1.0 yet, the releases are stable enough for everyday use and are being used in a number of research projects and neat applications by a growing community of users who are early adopters. JuMP remains under active development, and we welcome your feedback, suggestions, and bug reports."
},

{
    "location": "index.html#Contents-1",
    "page": "Introduction",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"installation.md\",\n    \"quickstart.md\",\n    \"concepts.md\",\n    \"variables.md\",\n    \"expressions.md\",\n    \"constraints.md\",\n    \"solvers.md\",\n    \"nlp.md\",\n    \"style.md\",\n    \"extensions.md\",\n    \"updating.md\",\n    \"howdoi.md\"]\nDepth = 2"
},

{
    "location": "index.html#Citing-JuMP-1",
    "page": "Introduction",
    "title": "Citing JuMP",
    "category": "section",
    "text": "If you find JuMP useful in your work, we kindly request that you cite the following paper (pdf):@article{DunningHuchetteLubin2017,\nauthor = {Iain Dunning and Joey Huchette and Miles Lubin},\ntitle = {JuMP: A Modeling Language for Mathematical Optimization},\njournal = {SIAM Review},\nvolume = {59},\nnumber = {2},\npages = {295-320},\nyear = {2017},\ndoi = {10.1137/15M1020575},\n}For an earlier work where we presented a prototype implementation of JuMP, see here:@article{LubinDunningIJOC,\nauthor = {Miles Lubin and Iain Dunning},\ntitle = {Computing in Operations Research Using Julia},\njournal = {INFORMS Journal on Computing},\nvolume = {27},\nnumber = {2},\npages = {238-248},\nyear = {2015},\ndoi = {10.1287/ijoc.2014.0623},\n}A preprint of this paper is freely available."
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
    "text": "DRAFT: Install Julia. Install JuMP."
},

{
    "location": "installation.html#Getting-Solvers-1",
    "page": "Installation Guide",
    "title": "Getting Solvers",
    "category": "section",
    "text": "TODO: Not yet updated for MOI.Solver support in Julia is currently provided by writing a solver-specific package that provides a very thin wrapper around the solver\'s C interface and providing a standard interface that JuMP can call. If you are interested in providing an interface to your solver, please get in touch. The table below lists the currently supported solvers and their capabilities.Solver Julia Package solver= License LP SOCP MILP NLP MINLP SDP\nArtelys Knitro KNITRO.jl KnitroSolver() Comm.    X X \nBARON BARON.jl BaronSolver() Comm.    X X \nBonmin AmplNLWriter.jl AmplNLWriter(CoinOptServices.bonmin) * EPL X  X X X \n\'\' CoinOptServices.jl OsilBonminSolver() \'\'      \nCbc Cbc.jl CbcSolver() EPL   X   \nClp Clp.jl ClpSolver() EPL X     \nCouenne AmplNLWriter.jl AmplNLWriter(CoinOptServices.couenne) * EPL X  X X X \n\'\' CoinOptServices.jl OsilCouenneSolver() \'\'      \nCPLEX CPLEX.jl CplexSolver() Comm. X X X   \nECOS ECOS.jl ECOSSolver() GPL X X    \nFICO Xpress Xpress.jl XpressSolver() Comm. X X X   \nGLPK GLPKMath... GLPKSolver[LP|MIP]() GPL X  X   \nGurobi Gurobi.jl GurobiSolver() Comm. X X X   \nIpopt Ipopt.jl IpoptSolver() EPL X   X  \nMOSEK Mosek.jl MosekSolver() Comm. X X X X  X\nNLopt NLopt.jl NLoptSolver() LGPL    X  \nSCS SCS.jl SCSSolver() MIT X X    XWhere:LP = Linear programming\nSOCP = Second-order conic programming (including problems with convex quadratic constraints and/or objective)\nMILP = Mixed-integer linear programming\nNLP = Nonlinear programming\nMINLP = Mixed-integer nonlinear programming\nSDP = Semidefinite programming* requires CoinOptServices installed, see below.To install Gurobi, for example, and use it with a JuMP model m, run:Pkg.add(\"Gurobi\")\nusing JuMP\nusing Gurobi\n\nm = Model(solver=GurobiSolver())Setting solver options is discussed in the Model &lt;ref-model&gt; section.Solver-specific notes follow below."
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
    "text": "ECOS can be used by JuMP to solve LPs and SOCPs. ECOS does not support general quadratic objectives or constraints, only second-order conic constraints specified by using norm or the quadratic form x\'x <= y^2."
},

{
    "location": "installation.html#FICO-Xpress-1",
    "page": "Installation Guide",
    "title": "FICO Xpress",
    "category": "section",
    "text": "Requires a working installation of Xpress with an active license (it is possible to get license for academic use, see FICO Academic Partner Program). Supports SOCP and \"SOS\" constraints. The interface is experimental, but it does pass all JuMP and MathProgBase tests. Callbacks are not yet supported.warning: Warning\nIf you are using 64-bit Xpress, you must use 64-bit Julia (and similarly with 32-bit Xpress)."
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
    "text": "TODO: Quick example of solving an LP and getting the solution back."
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
    "text": ""
},

{
    "location": "variables.html#Variables-1",
    "page": "Variables",
    "title": "Variables",
    "category": "section",
    "text": ""
},

{
    "location": "variables.html#What-is-a-JuMP-VariableRef?-1",
    "page": "Variables",
    "title": "What is a JuMP VariableRef?",
    "category": "section",
    "text": "DRAFT: A JuMP variable is a reference to an index in a model. It\'s a thin wrapper around MOI.VariableIndex. Variables have (1) names, and (2) attributes. Describe the different scopes of a variable (e.g., as Julia variables, lookup by string name, and lookup by symbol)."
},

{
    "location": "variables.html#The-@variable-macro-1",
    "page": "Variables",
    "title": "The @variable macro",
    "category": "section",
    "text": "DRAFT: Describe the complete syntax of the @variable macro. Anonymous versus named variables. Describe the three possible container types returned and how to use them (Array, JuMPArray, and Dict).How to delete variables."
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
    "location": "constraints.html#Constraints-1",
    "page": "Constraints",
    "title": "Constraints",
    "category": "section",
    "text": "DRAFT: Describe how constraints are represented (link to MOI docs). Constraints are very similar to variables in (1) how names work (2) how attributes work, and (3) the macro syntax for constructing them. They\'re a bit different because they\'re parameterized by function-set type. Describe constraints vs. ConstraintRefs. Describe JuMP.constraintobject. How to delete constraints. How to modify constraints by setting attributes and MOI.modifyconstraint!. Describe semidefinite constraints and symmetry handling. Refer to NLP docs for nonlinear constraints."
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
    "text": "TODO: Describe the connection between JuMP and solvers. Automatic vs. Manual mode. CachingOptimizer. How to set/change solvers. How to set parameters (solver specific and generic). Status codes. Accessing the result. How to accurately measure the solve time."
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
    "text": "TODO: This is out of date and has not yet been updated for MOI.JuMP has support for general smooth nonlinear (convex and nonconvex) optimization problems. JuMP is able to provide exact, sparse second-order derivatives to solvers. This information can improve solver accuracy and performance.Nonlinear objectives and constraints are specified by using the @NLobjective and @NLconstraint macros. The familiar sum() syntax is supported within these macros, as well as prod() which analogously represents the product of the terms within. Note that the @objective and @constraint macros (and corresponding functions) do not currently support nonlinear expressions. However, a model can contain a mix of linear, quadratic, and nonlinear constraints or objective functions. Starting points may be provided by using the start keyword argument to @variable. If a starting value is not provided for a variable, it will be set to the projection of zero onto the interval defined by the variable bounds. For nonconvex problems, the returned solution is only guaranteed to be locally optimal. Convexity detection is not currently provided.For example, we can solve the classical Rosenbrock problem (with a twist) as follows:using JuMP\nm = Model()\n@variable(m, x, start = 0.0)\n@variable(m, y, start = 0.0)\n\n@NLobjective(m, Min, (1-x)^2 + 100(y-x^2)^2)\n\nsolve(m)\nprintln(\"x = \", getvalue(x), \" y = \", getvalue(y))\n\n# adding a (linear) constraint\n@constraint(m, x + y == 10)\nsolve(m)\nprintln(\"x = \", getvalue(x), \" y = \", getvalue(y))Examples: optimal control, maximum likelihood estimation, and Hock-Schittkowski tests."
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
    "text": "JuMP\'s library of recognized univariate functions is derived from the Calculus.jl package. If you encounter a standard special function not currently supported by JuMP, consider contributing to the list of derivative rules there. In addition to this built-in list of functions, it is possible to register custom (user-defined) nonlinear functions to use within nonlinear expressions. JuMP does not support black-box optimization, so all user-defined functions must provide derivatives in some form. Fortunately, JuMP supports automatic differentiation of user-defined functions, a feature to our knowledge not available in any comparable modeling systems.Automatic differentiation is not finite differencing. JuMP\'s automatically computed derivatives are not subject to approximation error.JuMP uses ForwardDiff.jl to perform automatic differentiation; see the ForwardDiff.jl documentation for a description of how to write a function suitable for automatic differentiation. The general guideline is to write code that is generic with respect to the number type; don\'t assume that the input to the function is Float64. To register a user-defined function with derivatives computed by automatic differentiation, use the JuMP.register method as in the following example:mysquare(x) = x^2\nmyf(x,y) = (x-1)^2+(y-2)^2\n\nm = Model()\n\nJuMP.register(m, :myf, 2, myf, autodiff=true)\nJuMP.register(m, :mysquare, 1, mysquare, autodiff=true)\n\n@variable(m, x[1:2] >= 0.5)\n@NLobjective(m, Min, myf(x[1],mysquare(x[2])))The above code creates a JuMP model with the objective function (x[1]-1)^2 + (x[2]^2-2)^2. The first argument to JuMP.register the model for which the functions are registered. The second argument is a Julia symbol object which serves as the name of the user-defined function in JuMP expressions; the JuMP name need not be the same as the name of the corresponding Julia method. The third argument specifies how many arguments the function takes. The fourth argument is the name of the Julia method which computes the function, and autodiff=true instructs JuMP to compute exact gradients automatically.note: Note\nAll arguments to user-defined functions are scalars, not vectors. To define a function which takes a large number of arguments, you may use the splatting syntax f(x...) = ....Forward-mode automatic differentiation as implemented by ForwardDiff.jl has a computational cost that scales linearly with the number of input dimensions. As such, it is not the most efficient way to compute gradients of user-defined functions if the number of input arguments is large. In this case, users may want to provide their own routines for evaluating gradients. The more general syntax for JuMP.register which accepts user-provided derivative evaluation routines is:JuMP.register(m::Model, s::Symbol, dimension::Integer, f::Function, ∇f::Function, ∇²f::Function)The input differs for functions which take a single input argument and functions which take more than one. For univariate functions, the derivative evaluation routines should return a number which represents the first and second-order derivatives respectively. For multivariate functions, the derivative evaluation routines will be passed a gradient vector which they must explicitly fill. Second-order derivatives of multivariate functions are not currently supported; this argument should be omitted. The following example sets up the same optimization problem as before, but now we explicitly provide evaluation routines for the user-defined functions:mysquare(x) = x^2\nmysquare_prime(x) = 2x\nmysquare_primeprime(x) = 2\n\nmyf(x,y) = (x-1)^2+(y-2)^2\nfunction ∇f(g,x,y)\n    g[1] = 2*(x-1)\n    g[2] = 2*(y-2)\nend\n\nm = Model()\n\nJuMP.register(m, :myf, 2, myf, ∇f)\nJuMP.register(m, :mysquare, 1, mysquare, mysquare_prime, mysquare_primeprime)\n\n@variable(m, x[1:2] >= 0.5)\n@NLobjective(m, Min, myf(x[1],mysquare(x[2])))"
},

{
    "location": "nlp.html#Factors-affecting-solution-time-1",
    "page": "Nonlinear Modeling",
    "title": "Factors affecting solution time",
    "category": "section",
    "text": "The execution time when solving a nonlinear programming problem can be divided into two parts, the time spent in the optimization algorithm (the solver) and the time spent evaluating the nonlinear functions and corresponding derivatives. Ipopt explicitly displays these two timings in its output, for example:Total CPU secs in IPOPT (w/o function evaluations)   =      7.412\nTotal CPU secs in NLP function evaluations           =      2.083For Ipopt in particular, one can improve the performance by installing advanced sparse linear algebra packages, see Installation Guide. For other solvers, see their respective documentation for performance tips.The function evaluation time, on the other hand, is the responsibility of the modeling language. JuMP computes derivatives by using the ReverseDiffSparse package, which implements, in pure Julia, reverse-mode automatic differentiation with graph coloring methods for exploiting sparsity of the Hessian matrix [1]. As a conservative bound, JuMP\'s performance here currently may be expected to be within a factor of 5 of AMPL\'s."
},

{
    "location": "nlp.html#Querying-derivatives-from-a-JuMP-model-1",
    "page": "Nonlinear Modeling",
    "title": "Querying derivatives from a JuMP model",
    "category": "section",
    "text": "For some advanced use cases, one may want to directly query the derivatives of a JuMP model instead of handing the problem off to a solver. Internally, JuMP implements the AbstractNLPEvaluator interface from MathProgBase. To obtain an NLP evaluator object from a JuMP model, use JuMP.NLPEvaluator. The linearindex method maps from JuMP variables to the variable indices at the MathProgBase level.For example:m = Model()\n@variable(m, x)\n@variable(m, y)\n\n@NLobjective(m, Min, sin(x) + sin(y))\nvalues = zeros(2)\nvalues[linearindex(x)] = 2.0\nvalues[linearindex(y)] = 3.0\n\nd = JuMP.NLPEvaluator(m)\nMathProgBase.initialize(d, [:Grad])\nobjval = MathProgBase.eval_f(d, values) # == sin(2.0) + sin(3.0)\n\n∇f = zeros(2)\nMathProgBase.eval_grad_f(d, ∇f, values)\n# ∇f[linearindex(x)] == cos(2.0)\n# ∇f[linearindex(y)] == cos(3.0)The ordering of constraints in a JuMP model corresponds to the following ordering at the MathProgBase nonlinear abstraction layer. There are three groups of constraints: linear, quadratic, and nonlinear. Linear and quadratic constraints, to be recognized as such, must be added with the @constraint macros. All constraints added with the @NLconstraint macros are treated as nonlinear constraints. Linear constraints are ordered first, then quadratic, then nonlinear. The linearindex method applied to a constraint reference object returns the index of the constraint within its corresponding constraint class. For example:m = Model()\n@variable(m, x)\n@constraint(m, cons1, x^2 <= 1)\n@constraint(m, cons2, x + 1 == 3)\n@NLconstraint(m, cons3, x + 5 == 10)\n\ntypeof(cons1) # JuMP.ConstraintRef{JuMP.Model,JuMP.GenericQuadConstraint{JuMP.GenericQuadExpr{Float64,JuMP.VariableRef}}} indicates a quadratic constraint\ntypeof(cons2) # JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}}} indicates a linear constraint\ntypeof(cons3) # JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.NonlinearExprData}} indicates a nonlinear constraint\nlinearindex(cons1) == linearindex(cons2) == linearindex(cons3) == 1When querying derivatives, cons2 will appear first, because it is the first linear constraint, then cons1, because it is the first quadratic constraint, then cons3, because it is the first nonlinear constraint. Note that for one-sided nonlinear constraints, JuMP subtracts any values on the right-hand side when computing expressions. In other words, one-sided nonlinear constraints are always transformed to have a right-hand side of zero.The JuMP.constraintbounds(m::Model) method returns the lower and upper bounds of all the constraints in the model, concatenated in the order discussed above.This method of querying derivatives directly from a JuMP model is convenient for interacting with the model in a structured way, e.g., for accessing derivatives of specific variables. For example, in statistical maximum likelihood estimation problems, one is often interested in the Hessian matrix at the optimal solution, which can be queried using the JuMP.NLPEvaluator.If you are writing a \"solver\", we highly encourage use of the MathProgBase nonlinear interface over querying derivatives using the above methods. These methods are provided for convenience but do not fully integrate with JuMP\'s solver infrastructure. In particular, they do not allow users to specify your solver to the Model() constructor nor to call it using solve() nor to populate the solution back into the model. Use of the MathProgBase interface also has the advantage of being independent of JuMP itself; users of MathProgBase solvers are free to implement their own evaluation routines instead of expressing their model in JuMP. You may use the JuMP.build method to ask JuMP to populate the \"solver\" without calling optimize!."
},

{
    "location": "nlp.html#Raw-expression-input-1",
    "page": "Nonlinear Modeling",
    "title": "Raw expression input",
    "category": "section",
    "text": "In addition to the @NLobjective and @NLconstraint macros, it is also possible to provide Julia Expr objects directly by using JuMP.setNLobjective and JuMP.addNLconstraint. This input form may be useful if the expressions are generated programmatically. JuMP variables should be spliced into the expression object. For example:@variable(m, 1 <= x[i=1:4] <= 5)\nJuMP.setNLobjective(m, :Min, :($(x[1])*$(x[4])*($(x[1])+$(x[2])+$(x[3])) + $(x[3])))\nJuMP.addNLconstraint(m, :($(x[1])*$(x[2])*$(x[3])*$(x[4]) >= 25))\n\n# Equivalent form using traditional JuMP macros:\n@NLobjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])\n@NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)See the Julia documentation for more examples and description of Julia expressions.[1]: Dunning, Huchette, and Lubin, \"JuMP: A Modeling Language for Mathematical Optimization\", arXiv."
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
    "text": "TODO: A style guide for JuMP, JuMP models and surrounding Julia code. Formatting, naming, use of macros, comments, TODOs, docstrings, etc."
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
