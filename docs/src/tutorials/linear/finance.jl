# Copyright (c) 2019 Arpit Bhatia and contributors                               #src
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

# # Financial modeling problems

# **This tutorial was originally contributed by Arpit Bhatia.**

# Optimization models play an increasingly important role in financial
# decisions. Many computational finance problems can be solved efficiently using
# modern optimization techniques.
#
# In this tutorial we will discuss two such examples taken from the book
# [Optimization Methods in Finance](https://doi.org/10.1017/9781107297340).

# This tutorial uses the following packages
using JuMP
import HiGHS

# ## Short-term financing

# Corporations routinely face the problem of financing short term cash
# commitments such as the following:

# **Month**        |Jan    |Feb    |Mar    |Apr    |May    |Jun
# :-----:          |:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
# **Net Cash Flow**|-150   |-100   |200    |-200   |50     |300

# Net cash flow requirements are given in thousands of dollars. The company has
# the following sources of funds:

# - A line of credit of up to \$100K at an interest rate of 1% per month,
# - In any one of the first three months, it can issue 90-day commercial paper
#   bearing a total interest of 2% for the 3-month period,
# - Excess funds can be invested at an interest rate of 0.3% per month.

# Our task is to find out the most economical way to use these 3 sources such
# that we end up with the most amount of money at the end of June.
#
# We model this problem in the following manner:

# We will use the following decision variables:
# - the amount $u_{i}$ drawn from the line of credit in month $i$
# - the amount $v_{i}$ of commercial paper issued in month $i$
# - the excess funds $w_{i}$ in month $i$

# Here we have three types of constraints:
# 1. for every month, cash inflow = cash outflow for each month
# 2. upper bounds on $u_{i}$
# 3. nonnegativity of the decision variables $u_{i}$, $v_{i}$ and $w_{i}$.

# Our objective will be to simply maximize the company's wealth in June, which
# say we represent with the variable $m$.

financing = Model(HiGHS.Optimizer)

@variables(financing, begin
    0 <= u[1:5] <= 100
    0 <= v[1:3]
    0 <= w[1:5]
    m
end)

@objective(financing, Max, m)

@constraints(
    financing,
    begin
        u[1] + v[1] - w[1] == 150 # January
        u[2] + v[2] - w[2] - 1.01u[1] + 1.003w[1] == 100 # February
        u[3] + v[3] - w[3] - 1.01u[2] + 1.003w[2] == -200 # March
        u[4] - w[4] - 1.02v[1] - 1.01u[3] + 1.003w[3] == 200 # April
        u[5] - w[5] - 1.02v[2] - 1.01u[4] + 1.003w[4] == -50 # May
        -m - 1.02v[3] - 1.01u[5] + 1.003w[5] == -300 # June
    end
)

optimize!(financing)

objective_value(financing)

# ## Combinatorial auctions

# In many auctions, the value that a bidder has for a set of items may not be
# the sum of the values that he has for individual items.
#
# Examples are equity trading, electricity markets, pollution right auctions and
# auctions for airport landing slots.
#
# To take this into account, combinatorial auctions allow the bidders to submit
# bids on combinations of items.

# Let $M=\{1,2, \ldots, m\}$ be the set of items that the auctioneer has to
# sell. A bid is a pair $B_{j}=\left(S_{j}, p_{j}\right)$ where
# $S_{j} \subseteq M$ is a nonempty set of items and $p_{j}$ is the price offer
# for this set.
#
# Suppose that the auctioneer has received $n$ bids
# $B_{1}, B_{2}, \ldots, B_{n}.$ The goal of this problem is to help an
# auctioneer determine the winners in order to maximize his revenue.

# We model this problem by taking a decision variable $y_{j}$ for every bid. We
# add a constraint that each item $i$ is sold at most once. This gives us the
# following model:

# ```math
# \begin{aligned}
# \max && \sum_{i=1}^{n} p_{j} y_{j} \\
# \text { s.t. }  && \sum_{j : i \in S_{j}} y_{j} \leq 1 && \forall i=\{1,2 \ldots m\} \\
# && y_{j} \in\{0,1\} && \forall j \in\{1,2 \ldots n\}
# \end{aligned}
# ```

bid_values = [6 3 12 12 8 16]
bid_items = [[1], [2], [3 4], [1 3], [2 4], [1 3 4]]

auction = Model(HiGHS.Optimizer)
@variable(auction, y[1:6], Bin)
@objective(auction, Max, sum(y' .* bid_values))
for i in 1:6
    @constraint(auction, sum(y[j] for j in 1:6 if i in bid_items[j]) <= 1)
end

optimize!(auction)

objective_value(auction)

#-

value.(y)
