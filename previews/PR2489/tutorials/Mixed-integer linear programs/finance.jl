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

# # Finance

# **Originally Contributed by**: Arpit Bhatia

# Optimization models play an increasingly important role in financial
# decisions. Many computational finance problems can be solved efficiently using
# modern optimization techniques.
#
# In this tutorial we will discuss 3 such examples.
#
# The first two are taken from the book [Optimization Methods in Finance](https://doi.org/10.1017/9781107297340),
# while the third is the famous Markowitz Portfolio Optimization problem with
# data from [lecture notes from a course taught at Georgia Tech by Shabir Ahmed](https://www2.isye.gatech.edu/~sahmed/isye6669/).

# This tutorial uses the following packages
using JuMP
import GLPK
import Ipopt
import Statistics

# ## Short Term Financing

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

# Our objective will be to simply maximimse the company's wealth in June, which
# say we represent with the variable $m$.

financing = Model(GLPK.Optimizer)

@variables(financing, begin
    0 <= u[1:5] <= 100
    0 <= v[1:3]
    0 <= w[1:5]
    m
end)

@objective(financing, Max, m)

@constraints(financing, begin
    u[1] + v[1] - w[1] == 150 # January
    u[2] + v[2] - w[2] - 1.01u[1] + 1.003w[1] == 100 # February
    u[3] + v[3] - w[3] - 1.01u[2] + 1.003w[2] == -200 # March
    u[4] - w[4] - 1.02v[1] - 1.01u[3] + 1.003w[3] == 200 # April
    u[5] - w[5] - 1.02v[2] - 1.01u[4] + 1.003w[4] == -50 # May
    -m - 1.02v[3] - 1.01u[5] + 1.003w[5] == -300 # June
end)

optimize!(financing)

objective_value(financing)

# ## Combinatorial Auctions

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

auction = Model(GLPK.Optimizer)
@variable(auction, y[1:6], Bin)
@objective(auction, Max, sum(y' .* bid_values))
for i in 1:6
    @constraint(auction, sum(y[j] for j in 1:6 if i in bid_items[j]) <= 1)
end

optimize!(auction)

objective_value(auction)

#-

value.(y)

# ## Portfolio Optimization

# Suppose we are considering investing 1000 dollars in three non-dividend paying
# stocks, IBM (IBM), Walmart (WMT), and Southern Electric (SEHI), for a
# one-month period.

# This means we will use the money to buy shares of the three stocks at the
# current market prices, hold these for one month, and sell the shares off at
# the prevailing market prices at the end of the month.
#
# As a rational investor, we hope to make some profit out of this endeavor,
# i.e., the return on our investment should be positive.
#
# Suppose we bought a stock at $p$ dollars per share in the beginning of the
# month, and sold it off at $s$ dollars per share at the end of the month. Then
# the one-month return on a share of the stock is $ \frac{s-p}{p} $.

# Since the stock prices are quite uncertain, so is the end-of-month return on
# our investment. Our goal is to invest in such a way that the expected
# end-of-month return is at least \$50 or 5%. Furthermore, we want to make
# sure that the “risk” of not achieving our desired return is minimum.

# Note that we are solving the problem under the following assumptions:
# 1. We can trade any continuum of shares.
# 2. No short-selling is allowed.
# 3. There are no transaction costs.

# We model this problem by taking decision variables $x_{i}, i=1,2,3,$ denoting
# the dollars invested in each of the 3 stocks.

# Let us denote by $\tilde{r}_{i}$ the random variable corresponding to the
# monthly return (increase in the stock price) per dollar for stock $i$.

# Then, the return (or profit) on $x_{i}$ dollars invested in stock $i$ is
# $\tilde{r}_{i} x_{i},$ and the total (random) return on our investment is
# $\sum_{i=1}^{3} \tilde{r}_{i} x_{i}.$ The expected return on our investment is
# then $\mathbb{E}\left[\sum_{i=1}^{3} \tilde{r}_{i} x_{i}\right]=\sum_{i=1}^{3} \overline{r}_{i} x_{i},$
# where $\overline{r}_{i}$ is the expected value of the $\tilde{r}_{i}.$

# Now we need to quantify the notion of “risk” in our investment.

# Markowitz, in his Nobel prize winning work, showed that a rational investor’s
# notion of minimizing risk can be closely approximated by minimizing the
# variance of the return of the investment portfolio. This variance is given by:

# $$
# \operatorname{Var}\left[\sum_{i=1}^{3} \tilde{r}_{i} x_{i}\right] = \sum_{i=1}^{3} \sum_{j=1}^{3} x_{i} x_{j} \sigma_{i j}
# $$

# where $\sigma_{i j}$ is the covariance of the return of stock $i$ with stock $j$.

# Note that the right hand side of the equation is the most reduced form of the
# expression and we have not shown the intermediate steps involved in getting to
# this form. We can also write this equation as:

# ```math
# \operatorname{Var}\left[\sum_{i=1}^{3} \tilde{r}_{i} x_{i}\right] =x^{T} Q x
# ```

# Where $Q$ is the covariance matrix for the random vector $\tilde{r}$.

# Finally, we can write the model as:

# $$
# \begin{align*}
# \min x^{T} Q x \\
# \text { s.t. } \sum_{i=1}^{3} x_{i} \leq 1000.00 \\
# \overline{r}^{T} x \geq 50.00 \\
# x \geq 0
# \end{align*}
# $$

# After that long discussion, lets now use JuMP to solve the portfolio
# optimization problem for the data given below.

# | Month        |  IBM     |  WMT    |  SEHI  |
# |--------------|----------|---------|--------|
# | November-00  |  93.043  |  51.826 |  1.063 |
# | December-00  |  84.585  |  52.823 |  0.938 |
# | January-01   |  111.453 |  56.477 |  1.000 |
# | February-01  |  99.525  |  49.805 |  0.938 |
# | March-01     |  95.819  |  50.287 |  1.438 |
# | April-01     |  114.708 |  51.521 |  1.700 |
# | May-01       |  111.515 |  51.531 |  2.540 |
# | June-01      |  113.211 |  48.664 |  2.390 |
# | July-01      |  104.942 |  55.744 |  3.120 |
# | August-01    |  99.827  |  47.916 |  2.980 |
# | September-01 |  91.607  |  49.438 |  1.900 |
# | October-01   |  107.937 |  51.336 |  1.750 |
# | November-01  |  115.590 |  55.081 |  1.800 |

stock_data = [
    93.043 51.826 1.063;
    84.585 52.823 0.938;
    111.453 56.477 1.000;
    99.525 49.805 0.938;
    95.819 50.287 1.438;
    114.708 51.521 1.700;
    111.515 51.531 2.540;
    113.211 48.664 2.390;
    104.942 55.744 3.120;
    99.827 47.916 2.980;
    91.607 49.438 1.900;
    107.937 51.336 1.750;
    115.590 55.081 1.800;
]

# Calculating stock returns

stock_returns = Array{Float64}(undef, 12, 3)
for i in 1:12
    stock_returns[i, :] = (stock_data[i + 1, :] .- stock_data[i, :]) ./ stock_data[i, :]
end
stock_returns

# Calculating the expected value of monthly return:

r = Statistics.mean(stock_returns, dims = 1)

# Calculating the covariance matrix Q

Q = Statistics.cov(stock_returns)

# JuMP Model

portfolio = Model(Ipopt.Optimizer)
set_silent(portfolio)
@variable(portfolio, x[1:3] >= 0)
@objective(portfolio, Min, x' * Q * x)
@constraint(portfolio, sum(x) <= 1000)
@constraint(portfolio, sum(r[i] * x[i] for i = 1:3) >= 50)
optimize!(portfolio)

objective_value(portfolio)

#-

value.(x)
