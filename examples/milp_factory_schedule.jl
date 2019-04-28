#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################


#=
julia> VERSION
v"1.1.0"

pkg> status JuMP Cbc GLPK
  [9961bab8] Cbc v0.6.0
  [60bf3e95] GLPK v0.9.1
  [4076af6c] JuMP v0.19.0
=#

"""
This is a Julia translation of part 5 from "Introduction to to Linear Programming with Python" 
available at https://github.com/benalexkeen/Introduction-to-linear-programming

"""

using JuMP,GLPK,Test
const MOI = JuMP.MathOptInterface

function example_factory_schedule()

   d_max_cap = Dict(
      (9, :A)=>210000,
      (11, :A)=>80000,
      (6, :B)=>70000,
      (3, :A)=>120000,
      (8, :A)=>200000,
      (1, :A)=>100000,
      (7, :A)=>155000,
      (10, :B)=>100000,
      (5, :B)=>0,
      (7, :B)=>60000,
      (6, :A)=>140000,
      (12, :B)=>150000,
      (2, :A)=>110000,
      (5, :A)=>160000,
      (8, :B)=>100000,
      (9, :B)=>100000,
      (4, :B)=>100000,
      (11, :B)=>120000,
      (2, :B)=>55000,
      (3, :B)=>60000,
      (4, :A)=>145000,
      (1, :B)=>50000,
      (12, :A)=>150000,
      (10, :A)=>197000)

   d_min_cap = Dict(
      (9, :A)=>20000,
      (11, :A)=>20000,
      (6, :B)=>20000,
      (3, :A)=>20000,
      (8, :A)=>20000,
      (1, :A)=>20000,
      (7, :A)=>20000,
      (10, :B)=>20000,
      (5, :B)=>0,
      (7, :B)=>20000,
      (6, :A)=>20000,
      (12, :B)=>20000,
      (2, :A)=>20000,
      (5, :A)=>20000,
      (8, :B)=>20000,
      (9, :B)=>20000,
      (4, :B)=>20000,
      (11, :B)=>20000,
      (2, :B)=>20000,
      (3, :B)=>20000,
      (4, :A)=>20000,
      (1, :B)=>20000,
      (12, :A)=>20000,
      (10, :A)=>20000)

   d_var_cost = Dict(
      (9, :A)=>9,
      (11, :A)=>8,
      (6, :B)=>6,
      (3, :A)=>12,
      (8, :A)=>7,
      (1, :A)=>10,
      (7, :A)=>5,
      (10, :B)=>11,
      (5, :B)=>0,
      (7, :B)=>4,
      (6, :A)=>8,
      (12, :B)=>12,
      (2, :A)=>11,
      (5, :A)=>8,
      (8, :B)=>6,
      (9, :B)=>8,
      (4, :B)=>5,
      (11, :B)=>10,
      (2, :B)=>4,
      (3, :B)=>3,
      (4, :A)=>9,
      (1, :B)=>5,
      (12, :A)=>8,
      (10, :A)=>10)

   d_fixed_cost = Dict(
      (9, :A)=>500,
      (11, :A)=>500,
      (6, :B)=>600,
      (3, :A)=>500,
      (8, :A)=>500,
      (1, :A)=>500,
      (7, :A)=>500,
      (10, :B)=>600,
      (5, :B)=>0,
      (7, :B)=>600,
      (6, :A)=>500,
      (12, :B)=>600,
      (2, :A)=>500,
      (5, :A)=>500,
      (8, :B)=>600,
      (9, :B)=>600,
      (4, :B)=>600,
      (11, :B)=>600,
      (2, :B)=>600,
      (3, :B)=>600,
      (4, :A)=>500,
      (1, :B)=>600,
      (12, :A)=>500,
      (10, :A)=>500)

   d_demand = Dict(
      2=>100000,
      11=>140000,
      7=>150000,
      9=>200000,
      10=>190000,
      8=>170000,
      6=>130000,
      4=>130000,
      3=>130000,
      5=>140000,
      12=>100000,
      1=>120000)
   
   months = 1:12
   factories = [:A, :B]

   model = Model(with_optimizer(GLPK.Optimizer))

   @variables(model, begin
      status[m in months,f in factories],Bin
      production[ m in months,f in factories], Int
   end)

   #production cannot be less than minimum capacity
   @constraint(model, [m in months, f in factories],
      production[m,f] >= d_min_cap[m,f] * status[m,f])

   #production cannot be more tha maximum capacity
   @constraint(model,[m in months, f in factories],
      production[m,f] <= d_max_cap[m,f] * status[m,f])

   #production must equal demand
   @constraint(model,[m in months],
      production[m,:A] + production[m,:B] == d_demand[m])

   #factory B is shut down during month 5
   @constraints(model, begin
      status[5,:B] == 0
      production[5,:B]==0
   end)

   #minimize the cost of production 
   @objective(model, Min, sum(
      d_fixed_cost[m,f] * status[m,f] + d_var_cost[m,f] * production[m,f] 
      for m in months, f in factories))

   optimize!(model)

   @test JuMP.termination_status(model) == MOI.OPTIMAL
   @test JuMP.objective_value(model) == 12906400.0

   #spot check individual values 
   @test value.(production)[1,:A] == 70000
   @test value.(status)[1,:A] == 1
   @test value.(status)[5,:B] == 0
   @test value.(production)[5,:B] == 0
end

example_factory_schedule()