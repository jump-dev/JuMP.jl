# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Optimal power flow

# This tutorial formulates and solves an optimal power flow problem,
# a much-studied nonlinear problem from the field of electrical engineering.
# Our particular focus is on obtaining a good estimate of the objective value
# at the global optimum through the use of semidefinite programming techniques.

# The tutorial highlights JuMP's ability to formulate problems involving
# complex-valued decision variables and complex matrix cones such as the 
# [`HermitianPSDCone`](@ref) object.
# See the [Complex number support](@ref) section of the manual for more details.

# For another example of modeling with complex decision variables
# see the [Quantum state discrimination](@ref) tutorial.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import LinearAlgebra: diag
import SparseArrays: sparse, spdiagm, sparsevec
import DataFrames: DataFrame, rename
import SCS
import Test

# ## Formulation

# Optimal power flow problems for electrical transmission typically pose the
# following question: what is the most cost-effective operation of electricity
# generators while meeting constraints on the safe limits of network components?

# We will use the 9 _bus_ network test case `case9mod` to explore this problem
# following the example of Krasko and Rebennack (2017).

# Here *bus* and network *node* are taken as analogous terms, as are *branch* and transmission *line*.

# The graph of the network shows three nodes each used for the different purposes of
# generation ``G`` (nodes 1, 2 and 3),
# trans-shipment (nodes 4, 6 and 8) 
# and demand ``L`` (5, 7, and 9).

# ![Nine Nodes](../../assets/case9mod.png)

# The key decision variables here are the real power injections ``P^G`` and
# reactive power injections ``Q^G``over the allowed range of the generators.
# All other buses must restrict their generation variables to 0.
# On the other hand these non-generator nodes have a fixed 
# real and reactive power demand, denoted ``P^L`` and ``Q^L`` respectively
# (fixed at 0 in the case of trans-shipment and generator nodes).

# The cost of operating each generator is here modeled as a quadratic function of its real power output;
# in our specific test case, the objective function to minimize is
# ```math
# \begin{align}
#     \min \;\; & \; 0.11 \;\; (P^G_1)^2 +   5 P^G_1 + 150  \\
#        + \;\; & 0.085 \; (P^G_2)^2 + 1.2 P^G_2 + 600  \\
#          \;\; & 0.1225 \;  (P^G_3)^2 +     P^G_3 + 335 \\
# \end{align}
# ```
# Let's create a basic JuMP model with these settings:
N = 9 # number of network nodes

import Ipopt
model = Model(Ipopt.Optimizer)

P_G_lb = sparsevec([1,2,3], [10, 10, 10], N)
P_G_ub = sparsevec([1,2,3], [250, 300, 270], N)
@variable(model, P_G_lb[i] <= P_G[i=1:N] <= P_G_ub[i])


@objective(model, Min,
      0.11*P_G[1]^2 +   5*P_G[1] + 150
+    0.085*P_G[2]^2 + 1.2*P_G[2] + 600
+   0.1225*P_G[3]^2 +     P_G[3] + 335)


# Even before solving, we can substitute the lower bound on each generator's real power range 
# (all 10, as it turns out in this case)

value(z -> Dict(zip(P_G, P_G_lb))[z], objective_function(model))

# to see that we can do no better than an objective cost of 1188.75.

# However, we have additional power flow constraints to satisfy.
# In fact, we can get a better quick estimate from the direct observation that the
# active power generated must meet or exceed the active power demand.

P_L = sparsevec([5, 7, 9], [54, 60, 75], N)
@constraint(model, sum(P_G) >= sum(P_L))
optimize!(model)
objective_value(model)

# Power must flow from one or more generation nodes through the transmission lines
# and end up at a demand node. The state variables of our steady-state alternating current (AC)
# electrical network are *complex-valued* voltage variables ``V_1, \ldots, V_9``. Voltages capture both a magnitude and phase 
# of the node's electrical state in relation to the rest of the system.
# An AC power system also extends the notion of resistance in wires found in a direct current (DC) circuit to a 
# complex quantity, known as the *impedance*, of each transmission line.
# The reciprocal of impedance is known as *admittance*.
# Together these complex quantities are used to express a complex version of *Ohm's law*: current flow
# through a line is proportional to the difference in voltages on each end of the line
# multiplied by the admittance.

# ## Data

# Let's assemble the data we need for writing the complex power flow constraints. 
# The data for the problem consists of a list of the real and imaginary parts of the line impedance.
# We obtain from the `case9mod` MATPOWER test case `branch data` the following data table:
#! format: off 
ColName=[ :F_BUS, :T_BUS, :BR_R  ,:BR_X   ,:BR_Bc ]
  Lines=[(   1,      4,   0,      0.0576,  0,    )
         (   4,      5,   0.017,  0.092,   0.158 )
         (   6,      5,   0.039,  0.17,    0.358 )
         (   3,      6,   0,      0.0586,  0,    )
         (   6,      7,   0.0119, 0.1008,  0.209 )
         (   8,      7,   0.0085, 0.072,   0.149 )
         (   2,      8,   0,      0.0625,  0,    )
         (   8,      9,   0.032,  0.161,   0.306 )
         (   4,      9,   0.01,   0.085,   0.176 )

]
#! format: on
# Let's make this into a more accessible `DataFrame`:
df_br = rename(DataFrame(Lines), ColName)

# The first two columns describe the network, supplying the *from* and *to* connection points of the lines.
# The last three columns give the branch resistance, branch reactance and *line-charging susceptance*.

# We will also need to reference the `baseMVA` number (used for re-scaling):
baseMVA = 100;
# and the number of lines:
M = length(Lines)

# From the first two columns of the branch data table,
# we can create a sparse [incidence matrix](https://en.wikipedia.org/wiki/Incidence_matrix)
# that simplifies handling of the network layout:
A = sparse(df_br.F_BUS, 1:M, 1, N, M) + sparse(df_br.T_BUS, 1:M, -1, N, M)

# We form the network impedance vector
z = (df_br.BR_R .+ im * df_br.BR_X) / baseMVA
# and the branch line-charging susceptance
y_sh = 1/2 * (im * df_br.BR_Bc) * baseMVA
# and then the *bus admittance* matrix is defined as
Y = A * spdiagm(1 ./ z) * A' + spdiagm(diag( A * spdiagm(y_sh) * A' ))

# (The second term looks more complicated because we only want to add the diagonal elements in the calculation;
# the line-charging is used only in the nodal voltage terms and not the line voltage terms.)

# ## JuMP model

# Now we're ready to write the complex power flow constraints we need
# to more accurately model the physics of the electricity system.
@variable(model, V[1:N] in ComplexPlane())

@constraint(model, [i=1:N], 0.9 <= real(V[i])^2 + imag(V[i])^2 <= 1.1)

@constraint(model, [i=1:N], P_G[i] - P_L[i] == real(V[i]'*(Y*V)[i]))
@constraint(model, [i=1:N], Q_G[i] - Q_L[i] == imag(V[i]'*(Y*V)[i]))

@constraint(model, [i=1:N], S_G[i] - S_L[i] == V[i]'*(Y*V)[i])

optimize!(model)

# ## Results
# ## Conclusion

# ## References and further resources

# Krasko, Vitaliy, and Steffen Rebennack. “Chapter 15: Global Optimization: Optimal Power Flow Problem.” 
# In Advances and Trends in Optimization with Engineering Applications, 187–205. MOS-SIAM Series on Optimization.
# Society for Industrial and Applied Mathematics, 2017. https://doi.org/10.1137/1.9781611974683.ch15.

# [Test case `case9mod`](https://www.maths.ed.ac.uk/optenergy/LocalOpt/9busnetwork.html) from the
# [Test Case Archive of Optimal Power Flow (OPF) Problems with Local Optima](https://www.maths.ed.ac.uk/optenergy/LocalOpt/)

# [MATPOWER manual](https://matpower.org/docs/MATPOWER-manual.pdf) (see especially Appendix B "Data File Format" and Table B-3)

# The Julia/JuMP package [PowerModels.jl](https://lanl-ansi.github.io/PowerModels.jl/stable/) provides
# an interface to a wide range of power flow formulations along with utilities for working with detailed network data.
