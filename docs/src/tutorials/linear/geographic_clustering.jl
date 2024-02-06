# Copyright (c) 2019 Matthew Help, Mathieu Tanneau and contributors              #src
#                                                                                #src
# Permission is hereby granted, free of charge, to any person obtaining a copy   #src
# of this software and associated documentation files (the "Software"), to deal  #src
# in the Software without restriction, including without limitation the rights   #src
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      #src
# copies of the Software, and to permit persons to whom the Software is          #src
# furnished to do so, subject to the following conditions:                       #src
#                                                                                #src
# The above copyright notice and this permission notice shall be included in all #src
# copies or substantial portions of the Software.                                #src
#                                                                                #src
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #src
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       #src
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    #src
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         #src
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  #src
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  #src
# SOFTWARE.                                                                      #src

# # Geographical clustering

# **This tutorial was originally contributed by Matthew Helm and Mathieu Tanneau.**

# The goal of this exercise is to cluster $n$ cities into $k$ groups, minimizing
# the total pairwise distance between cities *and* ensuring that the variance in
# the total populations of each group is relatively small.

# This tutorial uses the following packages:
using JuMP
import DataFrames
import HiGHS
import LinearAlgebra

# For this example, we'll use the 20 most populous cities in the United States.
cities = DataFrames.DataFrame(
    Union{String,Float64}[
        "New York, NY" 8.405 40.7127 -74.0059
        "Los Angeles, CA" 3.884 34.0522 -118.2436
        "Chicago, IL" 2.718 41.8781 -87.6297
        "Houston, TX" 2.195 29.7604 -95.3698
        "Philadelphia, PA" 1.553 39.9525 -75.1652
        "Phoenix, AZ" 1.513 33.4483 -112.0740
        "San Antonio, TX" 1.409 29.4241 -98.4936
        "San Diego, CA" 1.355 32.7157 -117.1610
        "Dallas, TX" 1.257 32.7766 -96.7969
        "San Jose, CA" 0.998 37.3382 -121.8863
        "Austin, TX" 0.885 30.2671 -97.7430
        "Indianapolis, IN" 0.843 39.7684 -86.1580
        "Jacksonville, FL" 0.842 30.3321 -81.6556
        "San Francisco, CA" 0.837 37.7749 -122.4194
        "Columbus, OH" 0.822 39.9611 -82.9987
        "Charlotte, NC" 0.792 35.2270 -80.8431
        "Fort Worth, TX" 0.792 32.7554 -97.3307
        "Detroit, MI" 0.688 42.3314 -83.0457
        "El Paso, TX" 0.674 31.7775 -106.4424
        "Memphis, TN" 0.653 35.1495 -90.0489
    ],
    ["city", "population", "lat", "lon"],
)

# ### Model Specifics

# We will cluster these 20 cities into 3 different groups and we will assume
# that the ideal or target population $P$ for a group is simply the total
# population of the 20 cities divided by 3:

n = size(cities, 1)
k = 3
P = sum(cities.population) / k

# ### Obtaining the distances between each city

# Let's compute the pairwise Haversine distance between each of the cities in
# our data set and store the result in a variable we'll call `dm`:

"""
    haversine(lat1, long1, lat2, long2, r = 6372.8)

Compute the haversine distance between two points on a sphere of radius `r`,
where the points are given by the latitude/longitude pairs `lat1/long1` and
`lat2/long2` (in degrees).
"""
function haversine(lat1, long1, lat2, long2, r = 6372.8)
    lat1, long1 = deg2rad(lat1), deg2rad(long1)
    lat2, long2 = deg2rad(lat2), deg2rad(long2)
    hav(a, b) = sin((b - a) / 2)^2
    inner_term = hav(lat1, lat2) + cos(lat1) * cos(lat2) * hav(long1, long2)
    d = 2 * r * asin(sqrt(inner_term))
    ## Round distance to nearest kilometer.
    return round(Int, d)
end

# Our distance matrix is symmetric so we'll convert it to a `LowerTriangular`
# matrix so that we can better interpret the objective value of our model:

dm = LinearAlgebra.LowerTriangular([
    haversine(cities.lat[i], cities.lon[i], cities.lat[j], cities.lon[j])
    for i in 1:n, j in 1:n
])

# ### Build the model

# Now that we have the basics taken  care of, we can set up our model, create
# decision variables, add constraints, and then solve.

# First, we'll set up a model that leverages the Cbc solver. Next, we'll set up
# a binary variable $x_{i,k}$ that takes the value $1$ if city $i$ is in group
# $k$ and $0$ otherwise. Each city must be in a group, so we'll add the
# constraint $\sum_k x_{i,k} = 1$ for every $i$.

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:n, 1:k], Bin)
@constraint(model, [i = 1:n], sum(x[i, :]) == 1);
## To reduce symmetry, we fix the first city to belong to the first group.
fix(x[1, 1], 1; force = true)

# The total population of a group $k$ is $Q_k = \sum_ix_{i,k}q_i$ where $q_i$ is
# simply the $i$-th value from the `population` column in our `cities` DataFrame.
# Let's add constraints so that $\alpha \leq (Q_k - P) \leq \beta$. We'll set
# $\alpha$ equal to $-3$ million and $\beta$ equal to $3$. By adjusting
# these thresholds you'll find that there is a tradeoff between having
# relatively even populations between groups and having geographically close
# cities within each group. In other words, the larger the absolute values of
# $\alpha$ and $\beta$, the closer together the cities in a group will be but
# the variance between the group populations will be higher.

@variable(model, -3 <= population_diff[1:k] <= 3)
@constraint(model, population_diff .== x' * cities.population .- P)

# Now we need to add one last binary variable $z_{i,j}$ to our model that we'll
# use to compute the total distance between the cities in our groups, defined as
# $\sum_{i,j}d_{i,j}z_{i,j}$. Variable $z_{i,j}$ will equal $1$ if cities $i$
# and $j$ are in the same group, and $0$ if they are not in the same group.

# To ensure that $z_{i,j} = 1$ if and only if cities $i$ and $j$ are in the same
# group, we add the constraints $z_{i,j} \geq x_{i,k} + x_{j,k} - 1$ for every
# pair $i,j$ and every $k$:

@variable(model, z[i = 1:n, j = 1:i], Bin)
for k in 1:k, i in 1:n, j in 1:i
    @constraint(model, z[i, j] >= x[i, k] + x[j, k] - 1)
end

# We can now add an objective to our model which will simply be to minimize the
# dot product of $z$ and our distance matrix, `dm`.

@objective(model, Min, sum(dm[i, j] * z[i, j] for i in 1:n, j in 1:i));

# We can then call `optimize!` and review the results.

optimize!(model)
@assert has_optimal_solution(model)

# ### Reviewing the Results

# Now that we have results, we can add a column to our `cities` DataFrame for
# the group and then loop through our $x$ variable to assign each city to its
# group. Once we have that, we can look at the total population for each group
# and also look at the cities in each group to verify that they are grouped by
# geographic proximity.

cities.group = zeros(n)

for i in 1:n, j in 1:k
    if round(Int, value(x[i, j])) == 1
        cities.group[i] = j
    end
end

for group in DataFrames.groupby(cities, :group)
    @show group
    println("")
    @show sum(group.population)
    println("")
end
