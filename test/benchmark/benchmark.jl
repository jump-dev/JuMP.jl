# Particularly interesting graph coloring instance derived from an AC power flow model.
using HDF5, JLD
using ReverseDiffSparse
#using ProfileView

d = load("benchmark.jld");
I = d["I"];
J = d["J"];
n = d["n"];

g = gen_adjlist(zip(I,J), n)
color, num_colors = acyclic_coloring(g)
@time color, num_colors = acyclic_coloring(g)
#Profile.clear_malloc_data()
#@profile color, num_colors = acyclic_coloring(g)

#ProfileView.view()
