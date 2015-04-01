# Particularly interesting graph coloring instance derived from an AC power flow model.
using HDF5, JLD
using ReverseDiffSparse
#using ProfileView

d = load("benchmark2.jld");
I = d["I"];
J = d["J"];
n = d["n"];

g = gen_adjlist(I,J, n)
color, num_colors = acyclic_coloring(g)
@time color, num_colors = acyclic_coloring(g)
#@profile color, num_colors = acyclic_coloring(g)

rinfo = ReverseDiffSparse.recovery_preprocess(g, color, num_colors)
Profile.clear_malloc_data()
@time rinfo = ReverseDiffSparse.recovery_preprocess(g, color, num_colors)

#ProfileView.view()
