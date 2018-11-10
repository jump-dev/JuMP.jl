using JuMP
using Clp
using Compat

solver = Clp.Optimizer

# Nutrition guidelines
numCategories = 4
categories = ["calories", "protein", "fat", "sodium"]
minNutrition = [1800, 91, 0, 0]
maxNutrition = [2200, Inf, 65, 1779]

# Foods
numFoods = 9
foods = ["hamburger", "chicken", "hot dog", "fries",
                 "macaroni", "pizza", "salad", "milk", "ice cream"]
cost = [2.49, 2.89, 1.50, 1.89, 2.09, 1.99, 2.49, 0.89, 1.59]
nutritionValues = [410 24 26 730;
                   420 32 10 1190;
                   560 20 32 1800;
                   380  4 19 270;
                   320 12 10 930;
                   320 15 12 820;
                   320 31 12 1230;
                   100  8 2.5 125;
                   330  8 10 180]

# Build model
m = Model(with_optimizer(solver))

# Variables for nutrition info
@variable(m, minNutrition[i] <= nutrition[i=1:numCategories] <= maxNutrition[i])
# Variables for which foods to buy
@variable(m, buy[i=1:numFoods] >= 0)

# Objective - minimize cost
@objective(m, Min, Compat.dot(cost, buy))

# Nutrition constraints
for j = 1:numCategories
    @constraint(m, sum(nutritionValues[i,j]*buy[i] for i=1:numFoods) == nutrition[j])
end

@constraint(m, buy[8] + buy[9] <= 6)

# Solve
JuMP.optimize!(m)
