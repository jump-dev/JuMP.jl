#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# tsp.jl
#
# Solves the travelling salesman problem using integer programming and
# lazy generation of the subtour elimination constraints.
#############################################################################

using JuMP
using Gurobi
using Base.Test

# extractTour
# Given a n-by-n matrix representing the solution to an undirected TSP,
# extract the tour as a vector
# Input:
#  n        Number of cities
#  sol      n-by-n 0-1 symmetric matrix representing solution
# Output:
#  tour     n+1 length vector of tour, starting and ending at 1
function extractTour(n, sol)
    tour = [1]  # Start at city 1 always
    cur_city = 1
    while true
        # Look for first arc out of current city
        for j = 1:n
            if sol[cur_city,j] >= 1-1e-6
                # Found next city
                push!(tour, j)
                # Don't ever use this arc again
                sol[cur_city, j] = 0.0
                sol[j, cur_city] = 0.0
                # Move to next city
                cur_city = j
                break
            end
        end
        # If we have come back to 1, stop
        if cur_city == 1
            break
        end
    end  # end while
    return tour
end

# findSubtour
# Given a n-by-n matrix representing solution to the relaxed
# undirected TSP problem, find a set of nodes belonging to a subtour
# Input:
#  n        Number of cities
#  sol      n-by-n 0-1 symmetric matrix representing solution
# Outputs:
#  subtour  n length vector of booleans, true iff in a particular subtour
#  subtour_length   Number of cities in subtour (if n, no subtour found)
function findSubtour(n, sol)
    # Initialize to no subtour
    subtour = fill(false,n)
    # Always start looking at city 1
    cur_city = 1
    subtour[cur_city] = true
    subtour_length = 1
    while true
        # Find next node that we haven't yet visited
        found_city = false
        for j = 1:n
            if !subtour[j]
                if sol[cur_city, j] >= 1 - 1e-6
                    # Arc to unvisited city, follow it
                    cur_city = j
                    subtour[j] = true
                    found_city = true
                    subtour_length += 1
                    break  # Move on to next city
                end
            end
        end
        if !found_city
            # We are done
            break
        end
    end
    return subtour, subtour_length
end

# solveTSP
# Given a matrix of city locations, solve the TSP
# Inputs:
#   n       Number of cities
#   cities  n-by-2 matrix of (x,y) city locations
# Output:
#   path    Vector with order to cities are visited in
function solveTSP(n, cities)

    # Calculate pairwise distance matrix
    dist = zeros(n, n)
    for i = 1:n
        for j = i:n
            d = norm(cities[i,1:2] - cities[j,1:2])
            dist[i,j] = d
            dist[j,i] = d
        end
    end

    # Create a model that will use Gurobi to solve
    m = Model(solver=GurobiSolver())

    # x[i,j] is 1 iff we travel between i and j, 0 otherwise
    # Although we define all n^2 variables, we will only use
    # the upper triangle
    @variable(m, x[1:n,1:n], Bin)

    # Minimize length of tour
    @objective(m, Min, sum(dist[i,j]*x[i,j] for i=1:n for j=i:n))

    # Make x_ij and x_ji be the same thing (undirectional)
    # Don't allow self-arcs
    for i = 1:n
        @constraint(m, x[i,i] == 0)
        for j = (i+1):n
            @constraint(m, x[i,j] == x[j,i])
        end
    end

    # We must enter and leave every city once and only once
    for i = 1:n
        @constraint(m, sum(x[i,j] for j=1:n) == 2)
    end

    function subtour(cb)
        # Optional: display tour starting at city 1
        println("----\nInside subtour callback")
        println("Current tour starting at city 1:")
        print(extractTour(n, getvalue(x)))

        # Find any set of cities in a subtour
        subtour, subtour_length = findSubtour(n, getvalue(x))

        if subtour_length == n
            # This "subtour" is actually all cities, so we are done
            println("Solution visits all cities")
            println("----")
            return
        end

        # Subtour found - add lazy constraint
        # We will build it up piece-by-piece
        arcs_from_subtour = zero(AffExpr)

        for i = 1:n
            if !subtour[i]
                # If this city isn't in subtour, skip it
                continue
            end
            # Want to include all arcs from this city, which is in
            # the subtour, to all cities not in the subtour
            for j = 1:n
                if i == j
                    # Self-arc
                    continue
                elseif subtour[j]
                    # Both ends in same subtour
                    continue
                else
                    # j isn't in subtour
                    arcs_from_subtour += x[i,j]
                end
            end
        end

        # Add the new subtour elimination constraint we built
        println("Adding subtour elimination cut")
        println("----")
        @lazyconstraint(cb, arcs_from_subtour >= 2)
    end  # End function subtour

    # Solve the problem with our cut generator
    addlazycallback(m, subtour)
    solve(m)

    # Return best tour
    return extractTour(n, getvalue(x))
end  # end solveTSP


# Create a simple instance that looks like
#       +           +
#   +                   +
#       +           +
# The optimal tour is obvious, but the initial solution will be
#    /--+           +--\
#   +               |   +
#    \--+           +--/
n = 6
cities = [ 50 200;
                    100 100;
                    100 300;
                    500 100;
                    500 300;
                    550 200]
tour = solveTSP(n, cities)
println("Solution: ")
println(tour)
