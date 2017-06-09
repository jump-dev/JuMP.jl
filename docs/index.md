# JuMP
### Julia for Mathematical Programming

JuMP is a domain-specific modeling language for **[mathematical programming]** embedded in **[Julia]**. It currently supports a number of open-source and commercial solvers ([CPLEX], [COIN Clp], [COIN Cbc], [ECOS], [GLPK], [Gurobi], [Ipopt], [KNITRO], [MOSEK], [NLopt], [SCS]) for a variety of problem classes, including **[linear programming]**, **[(mixed) integer programming]**, **[second-order conic programming]**, **[semidefinite programming]**, and **[nonlinear programming]**.

JuMP makes it easy to specify and **solve optimization problems without expert knowledge**, yet at the same time allows experts to implement advanced algorithmic techniques such as exploiting efficient hot-starts in linear programming or using callbacks to interact with branch-and-bound solvers. JuMP is also **fast** - benchmarking has shown that it can create problems at similar speeds to special-purpose commercial tools such as AMPL while maintaining the expressiveness of a generic high-level programming language. JuMP can be easily embedded in complex work flows including simulations and web servers.

While neither Julia nor JuMP have reached version 1.0 yet, the releases are stable enough for everyday use and are being used in a number of research projects and neat applications by a growing community of users who are early adopters. JuMP remains under active development, and we welcome your feedback, suggestions, and bug reports.

Installing JuMP
---------------

If you are familiar with Julia you can get started quickly by using the
package manager to install JuMP::

    julia> Pkg.add("JuMP")

And a solver, e.g.::

    julia> Pkg.add("Clp")  # Will install Cbc as well

Then read the :ref:`quick-start` and/or see a :ref:`simple-example`.
The subsequent sections detail the complete functionality of JuMP.

-----------
Citing JuMP
-----------

If you find JuMP useful in your work, we kindly request that you cite the following `paper <http://dx.doi.org/10.1287/ijoc.2014.0623>`_:

.. code-block:: none

    @article{LubinDunningIJOC,
    author = {Miles Lubin and Iain Dunning},
    title = {Computing in Operations Research Using Julia},
    journal = {INFORMS Journal on Computing},
    volume = {27},
    number = {2},
    pages = {238-248},
    year = {2015},
    doi = {10.1287/ijoc.2014.0623},
    URL = {http://dx.doi.org/10.1287/ijoc.2014.0623}
    }

A preprint of this paper is freely available on `arXiv <http://arxiv.org/abs/1312.1431>`_.


[mathematical programming]: http://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: http://julialang.org/
[COIN Clp]: https://projects.coin-or.org/Clp
[COIN Cbc]: https://projects.coin-or.org/Cbc
[ECOS]: https://github.com/ifa-ethz/ecos
[GLPK]: http://www.gnu.org/software/glpk/
[Gurobi]: http://www.gurobi.com/
[MOSEK]: http://mosek.com/
[CPLEX]: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
[Ipopt]: https://projects.coin-or.org/Ipopt
[KNITRO]: http://www.ziena.com/knitro.htm
[NLopt]: http://ab-initio.mit.edu/wiki/index.php/NLopt
[SCS]: https://github.com/cvxgrp/scs
[linear programming]: http://en.wikipedia.org/wiki/Linear_programming
[(mixed) integer programming]: http://en.wikipedia.org/wiki/Integer_programming
[second-order conic programming]: http://en.wikipedia.org/wiki/Second-order_cone_programming
[semidefinite programming]: https://en.wikipedia.org/wiki/Semidefinite_programming
[nonlinear programming]: http://en.wikipedia.org/wiki/Nonlinear_programming