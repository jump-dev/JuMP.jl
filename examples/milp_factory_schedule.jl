#=
julia> VERSION
v"1.1.0"

pkg> status JuMP Cbc GLPK
  [9961bab8] Cbc v0.6.0
  [60bf3e95] GLPK v0.9.1
  [4076af6c] JuMP v0.19.0

This is a Julia translation of part 5 from "Introduction to to Linear Programming with Python" 
available at https://github.com/benalexkeen/Introduction-to-linear-programming

=#

using DelimitedFiles
vars=readdlm("data/factory_variables.csv",',',header=true)
dem = readdlm("data/monthly_demand.csv",',',header=true)

max_cap = Dict()
 for i in 1:size(vars[1],1)
    max_cap[(vars[1][i,:1],Symbol(vars[1][i,:2]))] = vars[1][i,:3]
 end

min_cap = Dict()
for i in 1:size(vars[1],1)
   min_cap[(vars[1][i,:1],Symbol(vars[1][i,:2]))] = vars[1][i,:4]
end

var_cost = Dict()
for i in 1:size(vars[1],1)
    var_cost[(vars[1][i,:1],Symbol(vars[1][i,:2]))] = vars[1][i,:5]
end

#fixed costs dictionary
fixed_cost = Dict()
for i in 1:size(vars[1],1)
    fixed_cost[(vars[1][i,:1],Symbol(vars[1][i,:2]))] = vars[1][i,:6]
end

#demand dictionary
demand = Dict()
for i in 1:size(dem[1],1)
    demand[(dem[1][i,:1])] = dem[1][i,:2]
end

#hardcoded, easier to reference?
months = 1:12
facts = [:A, :B]

using JuMP

#using Cbc
#model = Model(with_optimizer(Cbc.Optimizer))

using GLPK
model = Model(with_optimizer(GLPK.Optimizer))

@variables model begin
   status[m in months,f in facts],Bin
   production[ m in months,f in facts], Int
end

@constraint(model,[m in months, f in facts],
production[m,f] >= min_cap[m,f] * status[m,f]
)

@constraint(model,[m in months, f in facts],
production[m,f] <= max_cap[m,f] * status[m,f]
)

@constraint(model,[m in months],
production[m,:A] + production[m,:B] == demand[m]
)

@constraints(model, begin
status[5,:B] == 0
production[5,:B]==0
end)

@objective(model, Min, sum(
    fixed_cost[m,f] * status[m,f] + var_cost[m,f] * production[m,f] 
    for m in months, f in facts))

optimize!(model)

using Test

@test JuMP.termination_status(model) == MOI.OPTIMAL
@test JuMP.objective_value(model) == 12906400.0

#spot check individual values 
@test value.(production)[1,:A] == 70000
@test value.(status)[1,:A] == 1

@test value.(status)[5,:B] == 0
@test value.(production)[5,:B] == 0
