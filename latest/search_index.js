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
    "text": "(Image: Powered by NumFOCUS)# These comments do not display in the HTML output.\n# See https://github.com/JuliaDocs/Documenter.jl/issues/674.\n\n# Style conventions for the JuMP documentation:\n# - Respect the 80-character line limit whenever possible.\n# - Be concise.\n# - Use lists instead of long sentences.\n#   - Use numbered lists when describing a sequence, e.g., (1) do X, (2) then Y.\n#   - Use bullet points when the items are not ordered.\n# - Example code should be covered by doctests.\n#   - But it\'s unclear what to do if the code depends on a solver, see\n#     https://github.com/JuliaOpt/JuMP.jl/issues/1175.warning: Warning\nThis documentation is for the development version of JuMP. JuMP is undergoing a major transition to MathOptInterface, and the documentation is in the process of being rewritten. We recommend using the development version only for (1) developers of packages upstream or downstream of JuMP or (2) early adopters willing to provide feedback and file issues.JuMP is a domain-specific modeling language for mathematical optimization embedded in Julia. It currently supports a number of open-source and commercial solvers (see below) for a variety of problem classes, including linear programming, mixed-integer programming, second-order conic programming, semidefinite programming, and nonlinear programming. JuMP\'s features include:User friendliness\nSyntax that mimics natural mathematical expressions.\nComplete documentation (WIP!)\nSpeed\nBenchmarking has shown that JuMP can create problems at similar speeds   to special-purpose modeling languages such as   AMPL.\nJuMP communicates with most solvers in memory, avoiding the need to   write intermediary files.\nSolver independence\nJuMP uses a generic solver-independent interface provided by the   MathOptInterface   package, making it easy to change between a number of open-source and   commercial optimization software packages (\"solvers\").\nCurrently supported solvers include   Artelys Knitro,   Bonmin,   Cbc,   Clp,   Couenne,   CPLEX,   ECOS,   FICO Xpress,   GLPK,   Gurobi,   Ipopt,   MOSEK,   NLopt, and   SCS.\nAccess to advanced algorithmic techniques\nIncluding efficient LP re-solves which previously required using   solver-specific and/or low-level C++ libraries.\nEase of embedding\nJuMP itself is written purely in Julia. Solvers are the only binary   dependencies.\nBeing embedded in a general-purpose programming language makes it easy   to solve optimization problems as part of a larger workflow (e.g.,   inside a simulation, behind a web server, or as a subproblem in a   decomposition algorithm).\nAs a trade-off, JuMP\'s syntax is constrained by the syntax available   in Julia.\nJuMP is MPL licensed, meaning that   it can be embedded in commercial software that complies with the terms   of the license.While neither Julia nor JuMP have reached version 1.0 yet, the releases are stable enough for everyday use and are being used in a number of research projects and neat applications by a growing community of users who are early adopters. JuMP remains under active development, and we welcome your feedback, suggestions, and bug reports."
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
    "location": "variables.html#JuMP.@variable",
    "page": "Variables",
    "title": "JuMP.@variable",
    "category": "macro",
    "text": "@variable(model, kwargs...)\n\nAdd an anonymous (see Names) variable to the model model described by the keyword arguments kwargs and returns the variable.\n\n@variable(model, expr, args..., kwargs...)\n\nAdd a variable to the model model described by the expression expr, the positional arguments args and the keyword arguments kwargs. The expression expr can either be (note that in the following the symbol <= can be used instead of ≤ and the symbol >=can be used instead of ≥)\n\nof the form varexpr creating variables described by varexpr;\nof the form varexpr ≤ ub (resp. varexpr ≥ lb) creating variables described by varexpr with upper bounds given by ub (resp. lower bounds given by lb);\nof the form varexpr == value creating variables described by varexpr with fixed values given by value; or\nof the form lb ≤ varexpr ≤ ub or ub ≥ varexpr ≥ lb creating variables described by varexpr with lower bounds given by lb and upper bounds given by ub.\n\nThe expression varexpr can either be\n\nof the form varname creating a scalar real variable of name varname;\nof the form varname[...] or [...] creating a container of variables (see Containers in macro.\n\nThe recognized positional arguments in args are the following:\n\nBin: Sets the variable to be binary, i.e. either 0 or 1.\nInt: Sets the variable to be integer, i.e. one of ..., -2, -1, 0, 1, 2, ...\nSymmetric: Only available when creating a square matrix of variables, i.e. when varexpr is of the form varname[1:n,1:n] or varname[i=1:n,j=1:n]. It creates a symmetric matrix of variable, that is, it only creates a new variable for varname[i,j] with i ≤ j and sets varname[j,i] to the same variable as varname[i,j].\nPSD: The square matrix of variable is both Symmetric and constrained to be positive semidefinite.\n\nThe recognized keyword arguments in kwargs are the following:\n\nbasename: Sets the base name used to generate variable names. It corresponds to the variable name for scalar variable, otherwise, the variable names are basename[...] for each indices ... of the axes axes.\nlowerbound: Sets the value of the variable lower bound.\nupperbound: Sets the value of the variable upper bound.\nstart: Sets the variable starting value used as initial guess in optimization.\nbinary: Sets whether the variable is binary or not.\ninteger: Sets whether the variable is integer or not.\nvariabletype: See the \"Note for extending the variable macro\" section below.\ncontainer: Specify the container type, see Containers in macro.\n\nExamples\n\nThe following are equivalent ways of creating a variable x of name x with lowerbound 0:\n\n# Specify everything in `expr`\n@variable(model, x >= 0)\n# Specify the lower bound using a keyword argument\n@variable(model, x, lowerbound=0)\n# Specify everything in `kwargs`\nx = @variable(model, basename=\"x\", lowerbound=0)\n\nThe following are equivalent ways of creating a JuMPArray of index set [:a, :b] and with respective upper bounds 2 and 3 and names x[a] and `x[b].\n\nub = Dict(:a => 2, :b => 3)\n# Specify everything in `expr`\n@variable(model, x[i=keys(ub)] <= ub[i])\n# Specify the upper bound using a keyword argument\n@variable(model, x[i=keys(ub)], upperbound=ub[i])\n\nNote for extending the variable macro\n\nThe single scalar variable or each scalar variable of the container are created using addvariable(model, buildvariable(_error, info, extra_args...; extra_kwargs...)) where\n\nmodel is the model passed to the @variable macro;\n_error is an error function with a single String argument showing the @variable call in addition to the error message given as argument;\ninfo is the VariableInfo struct containing the information gathered in expr, the recognized keyword arguments (except basename and variabletype) and the recognized positional arguments (except Symmetric and PSD);\nextra_args are the unrecognized positional arguments of args plus the value of the variabletype keyword argument if present. The variabletype keyword argument allows the user to pass a position argument to buildvariable without the need to give a positional argument to @variable. In particular, this allows the user to give a positional argument to the buildvariable call when using the anonymous single variable syntax @variable(model, kwargs...); and\nextra_kwargs are the unrecognized keyword argument of kwargs.\n\nExamples\n\nThe following creates a variable x of name x with lowerbound 0 as with the first example above but does it without using the @variable macro\n\ninfo = VariableInfo(true, 0, false, NaN, false, NaN, false, NaN, false, false)\nJuMP.addvariable(model, JuMP.buildvariable(error, info), \"x\")\n\nThe following creates a JuMPArray of index set [:a, :b] and with respective upper bounds 2 and 3 and names x[a] and x[b] as with the second example above but does it without using the @variable macro\n\n# Without the `@variable` macro\ndata = Vector{JuMP.variabletype(model)}(undef, length(keys(ub)))\nx = JuMPArray(data, keys(ub))\nfor i in keys(ub)\n    info = VariableInfo(false, NaN, true, ub[i], false, NaN, false, NaN, false, false)\n    x[i] = JuMP.addvariable(model, JuMP.buildvariable(error, info), \"x[$i]\")\nend\n\nThe following are equivalent ways of creating a Matrix of size N x N with variables custom variables created with a JuMP extension using the Poly(X) positional argument to specify its variables:\n\n# Using the `@variable` macro\n@variable(model, x[1:N,1:N], Symmetric, Poly(X))\n# Without the `@variable` macro\nx = Matrix{JuMP.variabletype(model, Poly(X))}(N, N)\ninfo = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)\nfor i in 1:N, j in i:N\n    x[i,j] = x[j,i] = JuMP.addvariable(model, buildvariable(error, info, Poly(X)), \"x[$i,$j]\")\nend\n\n\n\n"
},

{
    "location": "variables.html#The-@variable-macro-1",
    "page": "Variables",
    "title": "The @variable macro",
    "category": "section",
    "text": "DRAFT: Describe the complete syntax of the @variable macro. Anonymous versus named variables. Describe the three possible container types returned and how to use them (Array, JuMPArray, and Dict).@variableHow to delete variables."
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
    "text": "@constraint(m::Model, expr)\n\nAdd a constraint described by the expression expr.\n\n@constraint(m::Model, ref[i=..., j=..., ...], expr)\n\nAdd a group of constraints described by the expression expr parametrized by i, j, ...\n\nThe expression expr can either be\n\nof the form func in set constraining the function func to belong to the set set, e.g. @constraint(m, [1, x-1, y-2] in MOI.SecondOrderCone(3)) constrains the norm of [x-1, y-2] be less than 1;\nof the form a sign b, where sign is one of ==, ≥, >=, ≤ and <= building the single constraint enforcing the comparison to hold for the expression a and b, e.g. @constraint(m, x^2 + y^2 == 1) constrains x and y to lie on the unit circle;\nof the form a ≤ b ≤ c or a ≥ b ≥ c (where ≤ and <= (resp. ≥ and >=) can be used interchangeably) constraining the paired the expression b to lie between a and c;\nof the forms @constraint(m, a .sign b) or @constraint(m, a .sign b .sign c) which broadcast the constraint creation to each element of the vectors.\n\nNote for extending the constraint macro\n\nEach constraint will be created using addconstraint(m, buildconstraint(_error, func, set)) where\n\n_error is an error function showing the constraint call in addition to the error message given as argument,\nfunc is the expression that is constrained\nand set is the set in which it is constrained to belong.\n\nFor expr of the first type (i.e. @constraint(m, func in set)), func and set are passed unchanged to buildconstraint but for the other types, they are determined from the expressions and signs. For instance, @constraint(m, x^2 + y^2 == 1) is transformed into addconstraint(m, buildconstraint(_error, x^2 + y^2, MOI.EqualTo(1.0))).\n\nTo extend JuMP to accept new constraints of this form, it is necessary to add the corresponding methods to buildconstraint. Note that this will likely mean that either func or set will be some custom type, rather than e.g. a Symbol, since we will likely want to dispatch on the type of the function or set appearing in the constraint.\n\n\n\n"
},

{
    "location": "constraints.html#JuMP.@SDconstraint",
    "page": "Constraints",
    "title": "JuMP.@SDconstraint",
    "category": "macro",
    "text": "@SDconstraint(m::Model, expr)\n\nAdd a semidefinite constraint described by the expression expr.\n\n@SDconstraint(m::Model, ref[i=..., j=..., ...], expr)\n\nAdd a group of semidefinite constraints described by the expression expr parametrized by i, j, ...\n\nThe expression expr needs to be of the form a sign b where sign is ⪰, ≥, >=, ⪯, ≤ or <= and a and b are square matrices. It constrains that a - b (or b - a if the sign is ⪯, ≤ or <=) is positive semidefinite.\n\n\n\n"
},

{
    "location": "constraints.html#Constraints-1",
    "page": "Constraints",
    "title": "Constraints",
    "category": "section",
    "text": "DRAFT: Describe how constraints are represented (link to MOI docs). Constraints are very similar to variables in (1) how names work (2) how attributes work, and (3) the macro syntax for constructing them. They\'re a bit different because they\'re parameterized by function-set type. Describe constraints vs. ConstraintRefs. Describe JuMP.constraintobject. How to delete constraints. How to modify constraints by setting attributes and MOI.modifyconstraint!. Describe semidefinite constraints and symmetry handling. Refer to NLP docs for nonlinear constraints.@constraint\n@SDconstraint"
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
    "text": "generatecontainer(T, indexvars, indexsets, requestedtype)\n\nReturn a tuple, the first element of which is code that generates a container for objects of type T given the index variables, index sets, and requestedtype. requestedtype may be one of :Array, :JuMPArray, :Dict, or :Auto. Return error-producing code if requested type is incompatible. For the case of :Auto, the following rules are used to determine the appropriate container:\n\nIf all index sets are either explicit 1:B objects for any B or symbols which refer to objects of type Base.OneTo, then an Array is generated of the appropriate size. Types of symbols/expressions are not known at compile time, so we defer to type-safe functions to check the Base.OneTo condition.\nIf condition (1) does not hold, and the index sets are independent (the index variable for one set does not appear in the definition of another), then an JuMPArray is generated of the appropriate size.\nOtherwise, generate an empty Dict{Any,T}.\n\nThe second element of the return tuple is a Bool, true if the container type automatically checks for duplicate terms in the index sets and false otherwise.\n\nExamples\n\ngeneratecontainer(VariableRef, [:i,:j], [:(1:N), :(1:T)], :Auto)\n# Returns code equivalent to:\n# :(Array{VariableRef}(length(1:N), length(1:T))\n\ngeneratecontainer(VariableRef, [:i,:j], [:(1:N), :(2:T)], :Auto)\n# Returns code equivalent to:\n# :(JuMPArray(Array{VariableRef}(length(1:N), length(2:T)), $indexvars...))\n\ngeneratecontainer(VariableRef, [:i,:j], [:(1:N), :(S)], :Auto)\n# Returns code that generates an Array if S is of type Base.OneTo,\n# otherwise an JuMPArray.\n\ngeneratecontainer(VariableRef, [:i,:j], [:(1:N), :(1:j)], :Auto)\n# Returns code equivalent to:\n# :(Dict{Any,VariableRef}())\n\n\n\n"
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
    "text": "There a two different naming aspects that need to be distinguished when creating variables/contraints (resp. a container of variables/constraints):The name of the local variable created (if any) holding the reference (resp. the container of references) which corresponds to the name that can be used to retrieve it using m[:name].\nThe name of the variable/constraint (resp. each variable/constraint in the container) used for printing. This corresponds to the MOI.VariableName/MOI.ConstraintName attribute.When creating a variable using the syntax @variable(m; kwargs...), creating a constraint using the syntax @constraint(m, expr) or when creating a container with the syntax [...] in a macro, we say that the variable or constraint is anonymous. For anonymous variables/constraints, no local variable is created holding the reference or container of references and it is not stored in the model, i.e. it is not possible to retrieve it using m[:name].Otherwise, when it is not anonymous, the name used both for the local variable created and the key for retrieving the reference or container of references in the model are determined from the macro expression. For instance, when creating a container with the syntax name[...] or when creating a constraint with @constraint(m, name, expr), the name used is name.The name of the variable/constraint used for printing is based on the base name which is specified by the basename keyword argument. When the basename keyword argument is not specified, the name depends on whether the variable is anonymous:if the variable/constraint is anonymous, then the MOI.VariableName/MOI.ConstraintName attribute is not set and the name used for printing is noname,\notherwise, the base name is set to the name used for the local variable created.The name of the variable/constraint set to the MOI.VariableName/MOI.ConstraintName attribute and used for printing is then basename for single variable/constraint and basename[i1,i2,...,in] for the reference at indices i1, i2, ..., in in a container."
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
    "text": "with_optimizer(constructor::Type, args...; kwargs...)\n\nReturn an OptimizerFactory that creates optimizers using the constructor constructor with positional arguments args and keyword arguments kwargs.\n\nExamples\n\nThe following returns an optimizer factory that creates IpoptOptimizers using the constructor call IpoptOptimizer(print_level=0):\n\nwith_optimizer(IpoptOptimizer, print_level=0)\n\n\n\n"
},

{
    "location": "solvers.html#JuMP.optimize",
    "page": "Solvers",
    "title": "JuMP.optimize",
    "category": "function",
    "text": "function optimize(model::Model,\n                  optimizer_factory::Union{Nothing, OptimizerFactory} = nothing;\n                  ignore_optimize_hook=(model.optimizehook===nothing))\n\nOptimize the model. If optimizer_factory is not nothing, it first set the optimizer to a new one created using the optimizer factory.\n\n\n\n"
},

{
    "location": "solvers.html#JuMP.Model-Tuple{}",
    "page": "Solvers",
    "title": "JuMP.Model",
    "category": "method",
    "text": "Model(; caching_mode::MOIU.CachingOptimizerMode=MOIU.Automatic,\n        bridge_constraints::Bool=true)\n\nReturn a new JuMP model without any optimizer; the model is stored the model in a cache. The mode of the CachingOptimizer storing this cache is caching_mode. The optimizer can be set later in the JuMP.optimize call. If bridge_constraints is true, constraints that are not supported by the optimizer are automatically bridged to equivalent supported constraints when an appropriate is defined in the MathOptInterface.Bridges module or is defined in another module and is explicitely added.\n\n\n\n"
},

{
    "location": "solvers.html#JuMP.Model-Tuple{JuMP.OptimizerFactory}",
    "page": "Solvers",
    "title": "JuMP.Model",
    "category": "method",
    "text": "Model(optimizer_factory::OptimizerFactory;\n      caching_mode::MOIU.CachingOptimizerMode=MOIU.Automatic,\n      bridge_constraints::Bool=true)\n\nReturn a new JuMP model using the optimizer factory optimizer_factory to create the optimizer. The optimizer factory can be created by the with_optimizer function.\n\nExamples\n\nThe following creates a model using the optimizer IpoptOptimizer(print_level=0):\n\nmodel = JuMP.Model(with_optimizer(IpoptOptimizer, print_level=0))\n\n\n\n"
},

{
    "location": "solvers.html#Automatic-and-Manual-modes-1",
    "page": "Solvers",
    "title": "Automatic and Manual modes",
    "category": "section",
    "text": "In Automatic and Manual modes, two MOI layers are automatically applied to the optimizer:CachingOptimizer: maintains a cache of the model so that when the optimizer does not support an incremental change to the model, the optimizer\'s internal model can be discarded and restored from the cache just before optimization. The CachingOptimizer has two different modes: Automatic and Manual corresponding to the two JuMP modes with the same names.\nLazyBridgeOptimizer (this can be disabled using the bridge_constraints keyword argument to Model constructor): when a constraint added is not supported by the optimizer, it tries transform the constraint into an equivalent form, possibly adding new variables and constraints that are supported by the optimizer. The applied transformations are selected among known recipes which are called bridges. A few default bridges are defined in MOI but new ones can be defined and added to the LazyBridgeOptimizer used by JuMP.See the MOI documentation for more details on these two MOI layers.To attach an optimizer to a JuMP model, JuMP needs to create a new empty optimizer instance. New optimizer instances can be obtained using an OptimizerFactory that can be created using the with_optimizer function:with_optimizerThe factory can be provided either at model construction time or at JuMP.optimize time:JuMP.optimizeNew JuMP models are created using the Model constructor:Model()\nModel(::JuMP.OptimizerFactory)TODO: how to control the caching optimizer states"
},

{
    "location": "solvers.html#JuMP.direct_model",
    "page": "Solvers",
    "title": "JuMP.direct_model",
    "category": "function",
    "text": "direct_model(backend::MOI.ModelLike)\n\nReturn a new JuMP model using backend to store the model and solve it. As opposed to the Model constructor, no cache of the model is stored outside of backend and no bridges are automatically applied to backend. The absence of cache reduces the memory footprint but it is important to bear in mind the following implications of creating models using this direct mode:\n\nWhen backend does not support an operation such as adding variables/constraints after solver or modifying constraints, an error is thrown. With models created using the Model constructor, such situations can be dealt with by storing the modifications in a cache and loading them into the optimizer when JuMP.optimize is called.\nNo constraint bridging is supported by default.\nThe optimizer used cannot be changed the model is constructed.\nThe model created cannot be copied.\n\n\n\n"
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
    "text": "This section describes the coding style rules that apply to JuMP code and that we recommend for JuMP models and surrounding Julia code. The motivations for a style guide include:conveying best practices for writing readable and maintainable code\nreducing the amount of time spent on bike-shedding by establishing basic naming and formatting conventions\nlowering the barrier for new contributors by codifying the existing practices (e.g., you can be more confident your code will pass review if you follow the style guide)In some cases, the JuMP style guide diverges from the Julia style guide. All such cases will be explicitly noted and justified.info: Info\nThe style guide is always a work in progress, and not all JuMP code follows the rules. When modifying JuMP, please fix the style violations of the surrounding code (i.e., leave the code tidier than when you started). If large changes are needed, consider separating them into another PR."
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
