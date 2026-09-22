# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Optimal power flow

# *This tutorial was originally contributed by James Foster (@jd-foster).*

# This tutorial formulates and solves an alternating current optimal power flow
# (AC-OPF) problem, a much-studied nonlinear problem from the field of
# electrical engineering.

# Once we've formulated and solved the nonlinear problem, we will turn our focus
# to obtaining a good estimate of the objective value at the global optimum
# through the use of semidefinite programming.

# One main purpose of this tutorial is to highlight JuMP's ability to directly
# formulate problems involving complex-valued decision variables and complex
# matrix cones such as the  [`HermitianPSDCone`](@ref) object.

# For another example of modeling with complex decision variables, see the
# [Quantum state discrimination](@ref) tutorial, and see the
# [Complex number support](@ref) section of the manual for more details.

# This tutorial takes a matrix-oriented approach focused on network nodes
# that simplifies the construction of semidefinite programs.
# Another approach is to formulate the problem focusing on network lines
# where it is easier to work with flow constraints. A general approach is
# provided by the Julia/JuMP package [PowerModels.jl](https://lanl-ansi.github.io/PowerModels.jl/stable/),
# an open-source framework to a broad range of power flow model formulations
# along with utilities for working with detailed network data.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import Clarabel
import DataFrames
import Ipopt
import LinearAlgebra
import SparseArrays
import Test

# ## Initial formulation

# Optimal power flow problems for electrical transmission typically pose the
# following question: what is the most cost-effective operation of electricity
# generators while meeting constraints on the safe limits of network components?

# We'll use the 9-_node_ network test case [`case9mod`](https://www.maths.ed.ac.uk/optenergy/LocalOpt/9busnetwork.html)
# to explore this problem.

# The graph of the network, shown here, has three nodes (or _buses_) each for
# the different purposes of generation ``G`` (nodes 1, 2, and 3), trans-shipment
# (nodes 4, 6, and 8), and demand ``D`` (nodes 5, 7, and 9).

# ![Nine Nodes](../../assets/case9mod.png)

# This example is a modified version of the [MATPOWER](https://matpower.org/) ([Zimmerman2011](@cite))
# test case `case9` ([archive](https://github.com/MATPOWER/matpower/tree/master/data))
# created by [Bukhsh2013](@cite) for their test case archive of
# optimal power flow problems with local optima. This test case is also
# extensively evaluated in [Krasko2017](@cite).

# Here *bus* and network *node* are taken as analogous terms, as are *branch*
# and transmission *line*.

# For future reference, let's name the number of nodes in the network:

N = 9;

# The network data can be summarised using a small number of arrays. Using the
# `sparsevec` function from the `SparseArrays` standard library package, we can
# give the indices and values of the non-zero data points:

## Real generation: lower (`lb`) and upper (`ub`) bounds
P_Gen_lb = SparseArrays.sparsevec([1, 2, 3], [10, 10, 10], N)
P_Gen_ub = SparseArrays.sparsevec([1, 2, 3], [250, 300, 270], N)
## Reactive generation: lower (`lb`) and upper (`ub`) bounds
Q_Gen_lb = SparseArrays.sparsevec([1, 2, 3], [-5, -5, -5], N)
Q_Gen_ub = SparseArrays.sparsevec([1, 2, 3], [300, 300, 300], N)
## Power demand levels (real, reactive, and complex form)
P_Demand = SparseArrays.sparsevec([5, 7, 9], [54, 60, 75], N)
Q_Demand = SparseArrays.sparsevec([5, 7, 9], [18, 21, 30], N)
S_Demand = P_Demand + im * Q_Demand

# The key decision variables are the real power injections ``P^G`` and reactive
# power injections ``Q^G``over the allowed range of the generators. All other
# buses must restrict their generation variables to 0. On the other hand, these
# non-generator nodes have a fixed  real and reactive power demand, denoted
# ``P^D`` and ``Q^D`` respectively (these are fixed at 0 in the case of
# trans-shipment and generator nodes).

# The cost of operating each generator is modeled as a quadratic function of
# its real power output; in our specific test case, the objective function to
# minimize is:
# ```math
# \begin{align}
#     \min      &    0.11 \;\; (P^G_1)^2 +   5 P^G_1 + 150  \\
#               & +  0.085 \; (P^G_2)^2 + 1.2 P^G_2 + 600  \\
#               & + 0.1225 \;  (P^G_3)^2 +     P^G_3 + 335 \\
# \end{align}
# ```
# Let's create an initial JuMP model with some of this data:

model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, P_Gen_lb[i] <= P_G[i in 1:N] <= P_Gen_ub[i])
@objective(
    model,
    Min,
    (0.11 * P_G[1]^2 + 5 * P_G[1] + 150) +
    (0.085 * P_G[2]^2 + 1.2 * P_G[2] + 600) +
    (0.1225 * P_G[3]^2 + P_G[3] + 335),
);

# Even before solving an optimization problem, we can estimate a lower bound on
# the best objective value by substituting the lower bound on each generator's
# real power range (all 10, as it turns out in this case):

basic_lower_bound = value(lower_bound, objective_function(model));
Test.@test isapprox(basic_lower_bound, 1188.75; atol = 1e-2)  #src
println("Objective value (basic lower bound) : $basic_lower_bound")

# to see that we can do no better than an objective cost of 1188.75.

# (Direct substitution works because a
# [quadratic function](https://en.wikipedia.org/wiki/Quadratic_function#Graph_of_the_univariate_function)
# of a single variable ``x`` with positive coefficients is strictly increasing
# for all ``x \geq 0``.)

# In fact, we can get a quick but even better estimate from the direct
# observation that the real power generated must meet or exceed the real power
# demand:

@constraint(model, sum(P_G) >= sum(P_Demand))
optimize!(model)
@assert is_solved_and_feasible(model)
better_lower_bound = round(objective_value(model); digits = 2)
println("Objective value (better lower bound): $better_lower_bound")

# However, there are additional power flow constraints that must be satisfied.

# Power must flow from one or more generation nodes through the transmission
# lines and end up at a demand node. The state variables of our steady-state
# alternating current (AC) electrical network are *complex-valued* voltage
# variables ``V_1, \ldots, V_N``. Voltages capture both a magnitude and phase of
# the node's electrical state in relation to the rest of the system. An AC power
# system also extends the notion of resistance in wires found in a direct
# current (DC) circuit to a  complex quantity, known as the *impedance*, of each
# transmission line. The reciprocal of impedance is known as *admittance*.
# Together, these complex quantities are used to express a complex version of
# *Ohm's law*: current flow through a line is proportional to the difference in
# voltages on each end of the line, multiplied by the admittance.

# ## Network data

# Let's assemble the data we need for writing the complex power flow constraints.
# The data for the problem consists of a list of the real and imaginary parts of
# the line impedance. We obtain the following data table from the `branch data`
# section of the `case9mod` MATPOWER format file:

branch_data = DataFrames.DataFrame([
    (1, 4, 0.0, 0.0576, 0.0),
    (4, 5, 0.017, 0.092, 0.158),
    (6, 5, 0.039, 0.17, 0.358),
    (3, 6, 0.0, 0.0586, 0.0),
    (6, 7, 0.0119, 0.1008, 0.209),
    (8, 7, 0.0085, 0.072, 0.149),
    (2, 8, 0.0, 0.0625, 0.0),
    (8, 9, 0.032, 0.161, 0.306),
    (4, 9, 0.01, 0.085, 0.176),
]);
DataFrames.rename!(branch_data, [:F_BUS, :T_BUS, :BR_R, :BR_X, :BR_Bc])

# The first two columns describe the network, supplying the *from* and *to*
# connection points of the lines. The last three columns give the branch
# resistance, branch reactance and *line-charging susceptance*.

# We will also need to reference the `base_MVA` number (used for re-scaling):

base_MVA = 100;

# and the number of lines:

M = size(branch_data, 1)

# From the first two columns of the branch data table, we can create a sparse
# [incidence matrix](https://en.wikipedia.org/wiki/Incidence_matrix) that
# simplifies handling of the network layout:

A =
    SparseArrays.sparse(branch_data.F_BUS, 1:M, 1, N, M) +
    SparseArrays.sparse(branch_data.T_BUS, 1:M, -1, N, M)

# We form the network impedance vector from the next two columns

z = (branch_data.BR_R .+ im * branch_data.BR_X) / base_MVA;

# and calculate it's corresponding *bus admittance* matrix as

Y_0 = A * SparseArrays.spdiagm(1 ./ z) * A';

# while the last column gives the branch line-charging susceptance

y_sh = 1 / 2 * (im * branch_data.BR_Bc) * base_MVA;

# and leads to the shunt admittance matrix

Y_sh = SparseArrays.spdiagm(
    LinearAlgebra.diag(A * SparseArrays.spdiagm(y_sh) * A'),
);

# (The construction of the shunt admittance matrix `Y_sh` looks
# somewhat more complicated than `Y_0` because we only want to add the
# diagonal elements in the calculation; the line-charging is used only in the
# nodal voltage terms and not the line voltage terms.)

# The full *bus admittance* matrix ``Y`` is then defined as

Y = Y_0 + Y_sh;

# ## JuMP model

# Now we're ready to write the complex power flow constraints we need to more
# accurately model the electricity system.

# We'll introduce a number of constraints that model both the physics and
# operational requirements.

# Let's start by initializing a new model:

model = Model(Ipopt.Optimizer)
set_silent(model)

# Then we'll create the nodal power generation variables:

@variable(
    model,
    S_G[i in 1:N] in ComplexPlane(),
    lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
    upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
)

# We need complex nodal voltages (the system state variables):

@variable(model, V[1:N] in ComplexPlane(), start = 1.0 + 0.0im)

# and operational constraints for maintaining voltage magnitude levels:
@constraint(model, [i in 1:N], 0.9^2 <= real(V[i])^2 + imag(V[i])^2 <= 1.1^2)

# We also need to fix an origin or _reference angle_ from which all other
# complex voltage angles (arguments) are determined.
# Here we will use node 1 as the nominated _reference bus_.
# Fixing the imaginary component of a reference bus to zero sets its complex
# voltage angle to 0, while constraining the real part to be non-negative
# disallows equivalent solutions that are just a reflection by 180 degrees:

@constraint(model, imag(V[1]) == 0);
@constraint(model, real(V[1]) >= 0);

# The power flow equations express a conservation of energy (power) principle,
# where power generated less the power consumed must balance the power exchanged
# with the network:

@constraint(model, S_G - S_Demand .== V .* conj(Y * V))

# As above, the objective function is a quadratic cost of real power:

P_G = real(S_G)
@objective(
    model,
    Min,
    (0.11 * P_G[1]^2 + 5 * P_G[1] + 150) +
    (0.085 * P_G[2]^2 + 1.2 * P_G[2] + 600) +
    (0.1225 * P_G[3]^2 + P_G[3] + 335),
);

# We're finally ready to solve our nonlinear AC-OPF problem:

optimize!(model)
@assert is_solved_and_feasible(model)
Test.@test isapprox(objective_value(model), 3087.84; atol = 1e-2)  #src
solution_summary(model)

#-

objval_solution = round(objective_value(model); digits = 2)
println("Objective value (feasible solution) : $(objval_solution)")

# The solution's power generation (in rectangular form)
# and complex voltage values (in polar form using degrees) are:

DataFrames.DataFrame(;
    Bus = 1:N,
    ComplexPowerGen = round.(value.(S_G); digits = 2),
    VoltageMagnitude = round.(abs.(value.(V)); digits = 2),
    VoltageAngle_Deg = round.(rad2deg.(angle.(value.(V))); digits = 2),
)

# ## Relaxations and better objective bounds

# The Ipopt solver uses an interior-point algorithm. It has local optimality
# guarantees, but is unable to certify whether the solution is globally optimal.
# The solution we found is indeed globally optimal. The work to verify this has
# been done in [Bukhsh2013](@cite) and [Krasko2017](@cite), and
# different solvers (such as Gurobi, SCIP and GLOMIQO) are also able to verify
# this.

# The techniques of *convex relaxations* can also be used to improve on our
# current best lower bound:

better_lower_bound

# To this end, observe that the nonlinear constraints in the AC-OPF formulation
# are quadratic equalities for power flow along with quadratic voltage
# inequalities.

# Let's linearize these constraints by first making the substitution
# ``W = V V^*``, where:
# ```math
# W = V V^* \quad \iff  \quad W_{ii} = | V_i |^2, \quad W_{ik} = V_i \; \overline{V_k}, \quad \forall i, \, k \in \{ 1, \ldots, N \}
# ```
# and where ``V^*`` is the [conjugate transpose](https://en.wikipedia.org/wiki/Conjugate_transpose) of ``V``.

# On the face of it, this turns a quadratic voltage bound constraint like:
# ```math
# v_L \leq |V_i |^2  \leq v_U, \quad i  \in \{ 1, \ldots, N \}
# ```
# for some real ``v_L`` and ``v_U`` into a simple two-sided bound:
# ```math
# v_L \leq W_{ii}  \leq v_U,
# ```
# while each quadratic expression for the nodal power term:
# ```math
#     S^{Node}_i = V_i \overline{(YV)_i}
# ```
# becomes the linear combination:
# ```math
#     S^{Node}_i = (E_{ii} Y^T) \bullet W.
# ```
# Here ``A \bullet B = \operatorname{tr}(A^* B)`` is the
# [Frobenius inner product](https://en.wikipedia.org/wiki/Frobenius_inner_product)
# of two complex matrices, while ``E_{kn}`` denotes the _matrix unit_ with a
# single nonzero entry of 1 in row ``k`` and column ``n``.

E(k, n) = SparseArrays.sparse([k], [n], 1, N, N);

# Of course, we've shifted the nonlinearity into the equality constraint
# ``W = V V^*``: it is this constraint we will now relax using a
# semidefinite programming approach.

# We will make use of complex voltages and relax
# ``W = V V^*`` to;
# ```math
# W \succeq V V^*,
# ```
# where the relation ``\succeq`` is the ordering in the Hermitian positive
# semidefinite cone.

# The above constraint is equivalent to:
# ```math
# \begin{bmatrix} 1 & V^* \\ V & W \\ \end{bmatrix} \succeq 0
# ```
# by the theory of the [Schur complement](https://en.wikipedia.org/wiki/Schur_complement).
# This matrix inequality implies a number of second-order cone constraints by
# taking certain ``2 \times 2`` minors of the matrix for each
# ``i  \in \{ 1, \ldots, N \}``:
# ```math
# \begin{bmatrix} 1 & V_i^* \\ V_i & W_{ii} \\ \end{bmatrix} \succeq 0,
# ```
# which is equivalent to the real second-order cone inequality:
# ```math
# \operatorname{real}(W_{ii}) \geq \operatorname{real}(V_i)^2 + \operatorname{imag}(V_i)^2.
# ```
# We include these implied constraints as well for demonstration purposes.

# Putting it all together we get the following semidefinite relaxation of the
# AC-OPF problem:

model = Model(Clarabel.Optimizer)
set_attribute(model, "tol_gap_rel", 1e-3)
set_attribute(model, "tol_feas", 1e-3)
set_attribute(model, "tol_ktratio", 1e-3)
@variable(
    model,
    S_G[i in 1:N] in ComplexPlane(),
    lower_bound = P_Gen_lb[i] + Q_Gen_lb[i] * im,
    upper_bound = P_Gen_ub[i] + Q_Gen_ub[i] * im,
)
@variable(model, W[1:N, 1:N] in HermitianPSDCone())
@variable(model, V[1:N] in ComplexPlane(), start = 1.0 + 0.0im)
@constraint(model, [i in 1:N], 0.9^2 <= real(W[i, i]) <= 1.1^2)
@constraint(model, real(V[1]) >= 0)
@constraint(model, imag(V[1]) == 0)
@constraint(model, 0.9 <= real(V[1]) <= 1.1)
@constraint(model, LinearAlgebra.Hermitian([1 V'; V W]) in HermitianPSDCone())
## 2 x 2 minor inequalities:
@constraint(
    model,
    [i in 1:N],
    [0.5, real(W[i, i]), real(V[i]), imag(V[i])] in RotatedSecondOrderCone()
)
@constraint(
    model,
    [i in 1:N],
    S_G[i] - S_Demand[i] == LinearAlgebra.tr((conj(Y) * E(i, i)) * W),
)
P_G = real(S_G)
@objective(
    model,
    Min,
    (0.11 * P_G[1]^2 + 5 * P_G[1] + 150) +
    (0.085 * P_G[2]^2 + 1.2 * P_G[2] + 600) +
    (0.1225 * P_G[3]^2 + P_G[3] + 335),
)
optimize!(model)

#-

Test.@test is_solved_and_feasible(model; allow_almost = true)
sdp_relaxation_lower_bound = round(objective_value(model); digits = 2)
Test.@test isapprox(sdp_relaxation_lower_bound, 2753.04; rtol = 1e-3)     #src
println(
    "Objective value (W & V relax. lower bound): $sdp_relaxation_lower_bound",
)

# We can more easily see solution values by rounding out noisy data:

W_1 = SparseArrays.sparse(round.(value.(W); digits = 2))

# and recover an approximation to the voltage variables as:

DataFrames.DataFrame(;
    Bus = 1:N,
    Magnitude = round.(abs.(value.(V)); digits = 2),
    AngleDeg = round.(rad2deg.(angle.(value.(V))); digits = 2),
)

# For further information on exploiting sparsity see [Jabr2012](@cite).

# This relaxation has the advantage that we can work directly with complex
# voltages to extend the formulation, strengthen the relaxation and gain
# additional approximate information about the voltage variables.
