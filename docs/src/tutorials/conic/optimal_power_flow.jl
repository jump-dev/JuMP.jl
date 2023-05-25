# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Optimal power flow

# *This tutorial was originally contributed by James Foster (@jd-foster).*

# This tutorial formulates and solves an alternating current optimal power flow (AC-OPF) problem,
# a much-studied nonlinear problem from the field of electrical engineering.
# Once we've formulated and solved the nonlinear problem, 
# we will turn our focus to obtaining a good estimate of the objective value
# at the global optimum through the use of semidefinite programming techniques.

# The tutorial highlights JuMP's ability to directly formulate problems involving
# complex-valued decision variables and complex matrix cones such as the 
# [`HermitianPSDCone`](@ref) object.
# See also the [Complex number support](@ref) section of the manual for more details.

# For another example of modeling with complex decision variables
# see the [Quantum state discrimination](@ref) tutorial.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import Ipopt
import SCS
import LinearAlgebra
import SparseArrays
import DataFrames
import Test

# ## Initial formulation

# Optimal power flow problems for electrical transmission typically pose the
# following question: what is the most cost-effective operation of electricity
# generators while meeting constraints on the safe limits of network components?

# The graph of the network shown here has nine nodes, with three each used for
# the different purposes of
# generation ``G`` (nodes 1, 2 and 3),
# trans-shipment (nodes 4, 6 and 8) 
# and demand ``D`` (5, 7, and 9).

# ![Nine Nodes](../../assets/case9mod.png)

# We're using the 9-_bus_ network test case `case9mod` to explore this problem.
# This test case is a version of the
# [MATPOWER](https://matpower.org/) test case `case9`
# ([archive](https://github.com/MATPOWER/matpower/tree/master/data))
# that has been modified by Bukhsh et al. (2013) for their
# test case archive of optimal power flow problems with local optima.
# This test case is also extensively evaluated in Krasko and Rebennack (2017).

# Here *bus* and network *node* are taken as analogous terms, as are *branch* and transmission *line*.

# For future reference, let's name the number of nodes in the network:
N = 9

# The network data can be summarised using a small number of arrays.
# With the `sparsevec` function from the `SparseArrays` standard library package
# we can just give the indices and values of the non-zero data points:

# _Generation power bounds: lower (`lb`) and upper (`ub`)_

# _Active_
P_Gen_lb = SparseArrays.sparsevec([1, 2, 3], [10, 10, 10], N)
P_Gen_ub = SparseArrays.sparsevec([1, 2, 3], [250, 300, 270], N)

# _Reactive_
Q_Gen_lb = SparseArrays.sparsevec([1, 2, 3], [-5, -5, -5], N)
Q_Gen_ub = SparseArrays.sparsevec([1, 2, 3], [300, 300, 300], N)

# _Power demand levels (real, reactive, and complex form)_
P_Demand = SparseArrays.sparsevec([5, 7, 9], [54, 60, 75], N)
Q_Demand = SparseArrays.sparsevec([5, 7, 9], [18, 21, 30], N)
S_Demand = P_Demand + im * Q_Demand

# The key decision variables here are the real power injections ``P^G`` and
# reactive power injections ``Q^G``over the allowed range of the generators.
# All other buses must restrict their generation variables to 0.
# On the other hand these non-generator nodes have a fixed 
# real and reactive power demand, denoted ``P^D`` and ``Q^D`` respectively
# (these are fixed at 0 in the case of trans-shipment and generator nodes).

# The cost of operating each generator is here modeled as a quadratic function of its real power output;
# in our specific test case, the objective function to minimize is
# ```math
# \begin{align}
#     \min \;\; & \; 0.11 \;\; (P^G_1)^2 +   5 P^G_1 + 150  \\
#        + \;\; & 0.085 \; (P^G_2)^2 + 1.2 P^G_2 + 600  \\
#          \;\; & 0.1225 \;  (P^G_3)^2 +     P^G_3 + 335 \\
# \end{align}
# ```
# Let's create an initial JuMP model with some of this data:
model = Model(Ipopt.Optimizer)
set_silent(model)

@variable(model, P_Gen_lb[i] <= P_G[i = 1:N] <= P_Gen_ub[i])

@objective(
    model,
    Min,
    0.11 * P_G[1]^2 +
    5 * P_G[1] +
    150 +
    0.085 * P_G[2]^2 +
    1.2 * P_G[2] +
    600 +
    0.1225 * P_G[3]^2 +
    P_G[3] +
    335
)

# Even before solving an optimization problem, we can estimate a lower bound on the
# best objective value by substituting the lower bound on each generator's real power range 
# (all 10, as it turns out in this case)

objval_basic_lb = value(lower_bound, objective_function(model));
Test.@test isapprox(objval_basic_lb, 1188.75, atol = 1e-2)  #src
println("Objective value (basic lower bound): $(objval_basic_lb)")

# to see that we can do no better than an objective cost of 1188.75.

# (Direct substitution works because a
# [quadratic function](https://en.wikipedia.org/wiki/Quadratic_function#Graph_of_the_univariate_function)
# of a single variable ``x`` with positive coefficients is strictly increasing for all ``x \geq 0``.)

# In fact, we can get a quick but even better estimate from the direct observation that the
# real power generated must meet or exceed the real power demand.

@constraint(model, sum(P_G) >= sum(P_Demand))
optimize!(model)

objval_better_lb = round(objective_value(model); digits = 2)
println("Objective value (better lower bound): $(objval_better_lb)")

# However, there are additional power flow constraints that must be satisfied.

# Power must flow from one or more generation nodes through the transmission lines
# and end up at a demand node. The state variables of our steady-state alternating current (AC)
# electrical network are *complex-valued* voltage variables ``V_1, \ldots, V_N``. Voltages capture both a magnitude and phase 
# of the node's electrical state in relation to the rest of the system.
# An AC power system also extends the notion of resistance in wires found in a direct current (DC) circuit to a 
# complex quantity, known as the *impedance*, of each transmission line.
# The reciprocal of impedance is known as *admittance*.
# Together these complex quantities are used to express a complex version of *Ohm's law*: current flow
# through a line is proportional to the difference in voltages on each end of the line
# multiplied by the admittance.

# ## Network data

# Let's assemble the data we need for writing the complex power flow constraints. 
# The data for the problem consists of a list of the real and imaginary parts of the line impedance.
# We obtain the following data table from the `branch data` section of the `case9mod` MATPOWER format file:
df_br = DataFrames.DataFrame([
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
# Let's make this into a more accessible `DataFrame`:
DataFrames.rename!(df_br, [:F_BUS, :T_BUS, :BR_R, :BR_X, :BR_Bc])

# The first two columns describe the network, supplying the *from* and *to* connection points of the lines.
# The last three columns give the branch resistance, branch reactance and *line-charging susceptance*.

# We will also need to reference the `baseMVA` number (used for re-scaling):
baseMVA = 100;
# and the number of lines:
M = size(df_br, 1)

# From the first two columns of the branch data table,
# we can create a sparse [incidence matrix](https://en.wikipedia.org/wiki/Incidence_matrix)
# that simplifies handling of the network layout:
A =
    SparseArrays.sparse(df_br.F_BUS, 1:M, 1, N, M) +
    SparseArrays.sparse(df_br.T_BUS, 1:M, -1, N, M)

# We form the network impedance vector
z = (df_br.BR_R .+ im * df_br.BR_X) / baseMVA
# and the branch line-charging susceptance
y_sh = 1 / 2 * (im * df_br.BR_Bc) * baseMVA
# and then the *bus admittance* matrix is defined as
Y =
    A * SparseArrays.spdiagm(1 ./ z) * A' + SparseArrays.spdiagm(
        LinearAlgebra.diag(A * SparseArrays.spdiagm(y_sh) * A'),
    )

# (The second term looks more complicated because we only want to add the diagonal elements in the calculation;
# the line-charging is used only in the nodal voltage terms and not the line voltage terms.)

# ## JuMP model

# Now we're ready to write the complex power flow constraints we need
# to more accurately model the electricity system.

# We'll introduce a number of constraints that model both the physics and operational requirements.

# Let's start by initializing a new model:
model = Model(Ipopt.Optimizer)
set_silent(model)

# **Generation**

# Create the nodal power generation variables:
@variable(model, S_G[1:N] in ComplexPlane())
P_G = real(S_G);
Q_G = imag(S_G);

# Generators should operate over a prescribed range:
@constraint(model, [i = 1:N], P_Gen_lb[i] <= P_G[i] <= P_Gen_ub[i])
@constraint(model, [i = 1:N], Q_Gen_lb[i] <= Q_G[i] <= Q_Gen_ub[i])

# **Demand**

# Create the nodal power demand variables:
@variable(model, S_D[1:N] in ComplexPlane())

# The loads in this model are assumed to be fixed and of constant-power type:
@constraint(model, S_D .== S_Demand)

# **Operational state** 

# We need  complex nodal voltages (the system state variables):
@variable(model, V[1:N] in ComplexPlane(), start = 1.0 + 0.0im)

# Operational constraints for maintaining voltage magnitude levels:
@constraint(model, [i = 1:N], 0.9^2 <= real(V[i])^2 + imag(V[i])^2 <= 1.1^2)

# Fixing the imaginary component of a _slack bus_ to zero sets its complex voltage angle to 0,
# which serves as an origin or reference value for all other complex voltage angles.
# Here we're using node 1 as the nominated slack bus:
@constraint(model, imag(V[1]) == 0)

# **Power flow constraints**

# The current at a node is an expression representing a generalised version of Ohm's law and
# [Kirchhoff's circuit laws](https://en.wikipedia.org/wiki/Kirchhoff%27s_circuit_laws):
I_Node = Y * V

# This next expression represents the power exchanged with the network from each node.
# It is the product of nodal voltage and current but in its generalised complex form here:
S_Node = LinearAlgebra.diagm(V) * conj(I_Node)

# The power flow equations express a conservation of energy (power) principle, where
# power generated less the power consumed must balance the power exchanged with the network:
@constraint(model, S_G - S_D .== S_Node)

# **Objective**

# Quadratic cost of real power (as above):
@objective(
    model,
    Min,
    0.11 * P_G[1]^2 +
    5 * P_G[1] +
    150 +
    0.085 * P_G[2]^2 +
    1.2 * P_G[2] +
    600 +
    0.1225 * P_G[3]^2 +
    P_G[3] +
    335
)

# We're finally ready to solve our AC-OPF problem.
optimize!(model)
Test.@test isapprox(objective_value(model), 3087.84, atol = 1e-2)  #src
solution_summary(model)
objval_solution = round(objective_value(model); digits = 2)
println("Objective value (feasible solution): $(objval_solution))")

# We can see the solution's voltages:
DataFrames.DataFrame(;
    Bus = 1:N,
    Magnitude = round.(abs.(value.(V)), digits = 2),
    AngleDeg = round.(rad2deg.(angle.(value.(V))), digits = 2),
)

# ## Relaxations and better objective bounds
# The Ipopt solver uses an interior-point algorithm. It has local optimality guarantees, but is unable to certify
# whether the solution is globally optimal. The solution we found is indeed globally optimal.
# The work to verify this has been done in Bukhsh et al. (2013) and Krasko and Rebennack (2017),
# and different solvers (such as Gurobi, SCIP and GLOMIQO) are also able to verify this. 

# The techniques of *convex relaxations* can also be used to improve on our current best lower bound:
println("Objective value (better lower bound): $(objval_better_lb)")

# ## References and further resources

# **Krasko**, V., & S. **Rebennack**.
# [_Chapter 15: Global Optimization: Optimal Power Flow Problem._](https://doi.org/10.1137/1.9781611974683.ch15)
# In Advances and Trends in Optimization with Engineering Applications, 187—205. MOS-SIAM Series on Optimization.
# Society for Industrial and Applied Mathematics, 2017. 

# **Bukhsh**, W. A., Grothey, A., McKinnon, K. I., & Trodden, P. A.
# [_Local solutions of the optimal power flow problem._](https://doi.org/10.1109/TPWRS.2013.2274577)
# IEEE Transactions on Power Systems, 28(4), 4780-4788 (2013).

# [**Test case `case9mod`**](https://www.maths.ed.ac.uk/optenergy/LocalOpt/9busnetwork.html):
# from the
# [Test Case Archive of Optimal Power Flow (OPF) Problems with Local Optima](https://www.maths.ed.ac.uk/optenergy/LocalOpt/).

# R. D. Zimmerman, C. E. Murillo-Sanchez, & R. J. Thomas, 
# [_MATPOWER: Steady-State Operations, Planning and Analysis Tools for Power Systems
# Research and Education,_](https://doi.org/10.1109/TPWRS.2010.2051168)
# IEEE Transactions on Power Systems, 26(1), 12-19 (2011).

# **MATPOWER data format**:
# [MATPOWER manual](https://matpower.org/docs/MATPOWER-manual.pdf); see especially Appendix B "Data File Format" and Table B-3.

# **PowerModels.jl**:
# the Julia/JuMP package [PowerModels.jl](https://lanl-ansi.github.io/PowerModels.jl/stable/) provides
# an interface to a wide range of power flow formulations along with utilities for working with detailed network data.

# **Ipopt solver**:
# Wächter, A., Biegler, L. 
# [_On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming._](https://doi.org/10.1007/s10107-004-0559-y)
# Math. Program. 106, 25—57 (2006).
