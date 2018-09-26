#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite
# These files are JuMP implementations of the model described in
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems,
#  Springer, No, 187, 1981
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# This file will run all the HS test problems, allowing for easy testing
#############################################################################

if VERSION >= v"0.7-"
    import Pkg
    ipt = "Ipopt" in keys(Pkg.installed())
else
    ipt = Pkg.installed("Ipopt") !== nothing
end
if !ipt
    println("Cannot run NLP tests because Ipopt is not installed.")
    exit(1)
end
using Ipopt, JuMP, Compat, Compat.Test
nlp_solver = IpoptSolver(print_level=0)

HS_path = dirname(@__FILE__)

files = readdir(HS_path)

println("Running HS Tests")
for file in files
    occursin("runhs", file) && continue
    include(file)
end
