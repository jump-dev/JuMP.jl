.. _sdp:

---------------------
Semidefinite modeling
---------------------

JuMP has *experimental* support for general semidefinite optimization problems. 

Currently, `Mosek <http://mosek.com/>`_
is the only supported solver. To install Mosek, run::

    Pkg.add("Mosek")

Semidefinite Variables
^^^^^^^^^^^^^^^^^^^^^^

Semidefinite variables can be created with the ``@defSDPVar`` macro; the dimension of the matrix must be specified::
    
    @defSDPVar(m, X[3])

Variables are implicitly assumed to be positive semidefinite, but other bounds (with respect to the semidefinite partial order) can be explicitly stated::

    @defSDPVar(m, eye(5,5) <= Y[5,5] <= 5*ones(5,5))

Matrix Variables
^^^^^^^^^^^^^^^^

Matrix variables can be created with the ``@defMatrixVar`` macro; the dimension of the matrix must be specified. Matrix variables are allowed to be either one-dimensional (vectors) or two-dimensional (matrices)::

    @defMatrixVar(m, Z[2,3])

Component-wise variable bounds can be explicitly added::

    @defMatrixVar(m, ones(5,3) <= W[5,3] <= eye(5,3) + 2*ones(5,3))

Matrix Expressions
^^^^^^^^^^^^^^^^^^

Matrix expressions are a natural extension of affine expressions for matrix variables. Matrix variables can be multiplied by or added with scalar matrices, and concatenated to create larger structures::

    @defSDPVar(m, X[3])
    @defMatrixVar(m, Y[2,3])
    R = diagm(2ones(3)) +  diagm(ones(2),-1) + diagm(ones(2),1)
    matexpr = [3*ones(3,3)*X Y; Y' -R*Y] + eye(5,5)

SDP solvers are strongly reliant on data sparsity to solve large-scale problems; because of this, it is important to use sparse matrices whenever possible.

Matrix Function Variables
^^^^^^^^^^^^^^^^^^^^^^^^^

Matrix function variables are affine functions of matrix expressions. Supported operations include trace and element reference::

    trace(A*X)
    X[3,2]

Matrix Function Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Affine expressions can be constructed from a mixture of matrix function variables, scalar variables, and scalar constants::

    trace(A*X) + y + 3

Norm Expressions
^^^^^^^^^^^^^^^^

Primal Constraints
^^^^^^^^^^^^^^^^^^

Primal constraints are scalar linear constraints on both matrix and scalar variables: :math:``L \leq \sum_{i} f_i(X_i) + \sum_{j} c_jy_j \leq U``, with matrix variables :math:`X_i`, scalar variables :math:`y_i`, symmetric matrices :math:`A_i`, scalar constants :math:`c_j,L`, and :math:`U`, and matrix function variables :math:`f_i`. These constraints are the natural form for many modern SDP solvers.

Dual Constraints
^^^^^^^^^^^^^^^^

Dual constraints are semidefinite constraints on scalar variables: :math:``\sum_{i} y_iB_i \succcurlyeq C`` for (symmetric) constant matrices :math:`B_i` and :math:`C` and scalar variables :math:`y_i`. 

Matrix Constraints
^^^^^^^^^^^^^^^^^^

Matrix constraints are semidefinite constraints on matrix variables: :math:``\sum_{i} A_iX_iB_i \succcurlyeq C`` for matrices :math:`A_i,B_i`, and :math:`C` and matrix variables :math:`X_i`.

SDP Example
^^^^^^^^^^^

    m = Model()
    @defSDPVar(m, X[3] >= eye(3,3))
    @defVar(m, 0 <= y <= 1)

    addConstraint(m, X[1,2] + y == 5)
    addConstraint(m, eye(2,2)*y <= 2*ones(2,2))
    addConstraint(m, ones(3,3)*X <= 5*eye(3,3))
    setObjective(m, :Max, trace(X) + y)

    solve(m)
