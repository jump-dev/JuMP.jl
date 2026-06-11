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

# # Getting started with data and plotting

# In this tutorial we will learn how to read tabular data into Julia, and some
# of the basics of plotting.

# If you're new to Julia, start by reading [Getting started with Julia](@ref)
# and [Getting started with JuMP](@ref) first.

# !!! note
#     There are multiple ways to read the same kind of data into Julia. This
#     tutorial focuses on DataFrames.jl because it provides the ecosystem to
#     work with most of the required file types in a straightforward manner.

# Before we get started, we need this constant to point to where the data files
# are.

const DATA_DIR = joinpath(@__DIR__, "data")

# ## Where to get help

# Read the documentation
# * Plots.jl: [http://docs.juliaplots.org/latest/](http://docs.juliaplots.org/latest/)
# * CSV.jl: [http://csv.juliadata.org/stable](http://csv.juliadata.org/stable)
# * DataFrames.jl: [https://dataframes.juliadata.org/stable/](https://dataframes.juliadata.org/stable/)

# ## Preliminaries

# To get started, we need to install some packages.

# ### DataFrames.jl

# The `DataFrames` package provides a set of tools for working with tabular
# data. It is available through the Julia package manager.
# ```julia
# using Pkg
# Pkg.add("DataFrames")
# ```

import DataFrames

# !!! info "What is a DataFrame?"
#     A DataFrame is a data structure like a table or spreadsheet. You can use
#     it for storing and exploring a set of related data values. Think of it as
#     a smarter array for holding tabular data.

# ### Plots.jl

# The `Plots` package provides a set of tools for plotting. It is available
# through the Julia package manager.
# ```julia
# using Pkg
# Pkg.add("Plots")
# ```

import Plots

# ### CSV .jl

# CSV and other delimited text files can be read by the CSV.jl package.

# ```julia
# Pkg.add("CSV")
# ```

import CSV

# ## DataFrame basics

# To read a CSV file into a DataFrame, we use the `CSV.read` function.

csv_df = CSV.read(joinpath(DATA_DIR, "StarWars.csv"), DataFrames.DataFrame)

# Let's try plotting some of this data

Plots.scatter(
    csv_df.Weight,
    csv_df.Height;
    xlabel = "Weight",
    ylabel = "Height",
)

# That doesn't look right. What happened? If you look at the dataframe above, it
# read `Weight` in as a `String` column because there are "NA" fields. Let's
# correct that, by telling CSV to consider "NA" as `missing`.

csv_df = CSV.read(
    joinpath(DATA_DIR, "StarWars.csv"),
    DataFrames.DataFrame;
    missingstring = "NA",
)

# Then let's re-plot our data

Plots.scatter(
    csv_df.Weight,
    csv_df.Height;
    title = "Height vs Weight of StarWars characters",
    xlabel = "Weight",
    ylabel = "Height",
    label = false,
    ylims = (0, 3),
)

# Better!

# !!! tip
#     Read the [CSV documentation](https://csv.juliadata.org/stable/) for other
#     parsing options.

# DataFrames.jl supports manipulation using functions similar to pandas. For
# example, split the dataframe into groups based on eye-color:

by_eyecolor = DataFrames.groupby(csv_df, :Eyecolor)

# Then recombine into a single dataframe based on a function operating over the
# split dataframes:

eyecolor_count = DataFrames.combine(by_eyecolor) do df
    return DataFrames.nrow(df)
end

# We can rename columns:

DataFrames.rename!(eyecolor_count, :x1 => :count)

# Drop some missing rows:

DataFrames.dropmissing!(eyecolor_count, :Eyecolor)

# Then we can visualize the data:

sort!(eyecolor_count, :count; rev = true)
Plots.bar(
    eyecolor_count.Eyecolor,
    eyecolor_count.count;
    xlabel = "Eye color",
    ylabel = "Number of characters",
    label = false,
)

# ### Other Delimited Files

# We can also use the `CSV.jl` package to read any other delimited text file
# format.

# By default, CSV.File will try to detect a file's delimiter from the first 10
# lines of the file.

# Candidate delimiters include `','`, `'\t'`, `' '`, `'|'`, `';'`, and `':'`. If
# it can't auto-detect the delimiter, it will assume `','`.

# Let's take the example of space separated data.

ss_df = CSV.read(joinpath(DATA_DIR, "Cereal.txt"), DataFrames.DataFrame)

# We can also specify the delimiter as follows:

delim_df = CSV.read(
    joinpath(DATA_DIR, "Soccer.txt"),
    DataFrames.DataFrame;
    delim = "::",
)

# ## Working with DataFrames

# Now that we have read the required data into a DataFrame, let us look at some
# basic operations we can perform on it.

# ### Querying Basic Information

# The `size` function gets us the dimensions of the DataFrame:

DataFrames.size(ss_df)

# We can also use the `nrow` and `ncol` functions to get the number of rows and
# columns respectively:

DataFrames.nrow(ss_df), DataFrames.ncol(ss_df)

# The `describe` function gives basic summary statistics of data in a DataFrame:

DataFrames.describe(ss_df)

# Names of every column can be obtained by the `names` function:

DataFrames.names(ss_df)

# Corresponding data types are obtained using the broadcasted `eltype` function:

eltype.(ss_df)

# ### Accessing the Data

# Similar to regular arrays, we use numerical indexing to access elements of a
# DataFrame:

csv_df[1, 1]

# The following are different ways to access a column:

csv_df[!, 1]

#-

csv_df[!, :Name]

#-

csv_df.Name

#-

csv_df[:, 1] # Note that this creates a copy.

# The following are different ways to access a row:

csv_df[1:1, :]

#-

csv_df[1, :] # This produces a DataFrameRow.

# We can change the values just as we normally assign values.

# Assign a range to scalar:

csv_df[1:3, :Height] .= 1.83

# Assign a vector:

csv_df[4:6, :Height] = [1.8, 1.6, 1.8]

#-

csv_df

# !!! tip
#     There are a lot more things which can be done with a DataFrame. Read the
#     [docs](https://juliadata.github.io/DataFrames.jl/stable/) for more
#     information.

# For information on dplyr-type syntax:
# * Read the [DataFrames.jl documentation](https://dataframes.juliadata.org/stable/man/querying_frameworks/)
# * Check out [DataFramesMeta.jl](https://github.com/JuliaData/DataFramesMeta.jl)

# ## Example: the passport problem

# Let's now apply what we have learned to solve a real problem.

# ### Data manipulation

# The [Passport Index Dataset](https://github.com/ilyankou/passport-index-dataset)
# lists travel visa requirements for 199 countries, in `.csv` format. Our task
# is to find the minimum number of passports required to visit all countries.

passport_data = CSV.read(
    joinpath(DATA_DIR, "passport-index-matrix.csv"),
    DataFrames.DataFrame,
);

# In this dataset, the first column represents a passport (=from) and each
# remaining column represents a foreign country (=to).

# The values in each cell are as follows:
# * 3 = visa-free travel
# * 2 = eTA is required
# * 1 = visa can be obtained on arrival
# * 0 = visa is required
# * -1 is for all instances where passport and destination are the same

# Our task is to find out the minimum number of passports needed to visit every
# country without requiring a visa.

# The values we are interested in are -1 and 3. Let's modify the dataframe so
# that the -1 and 3 are `1` (true), and all others are `0` (false):

function modifier(x)
    if x == -1 || x == 3
        return 1
    else
        return 0
    end
end

for country in passport_data.Passport
    passport_data[!, country] = modifier.(passport_data[!, country])
end

# The values in the cells now represent:
# * 1 = no visa required for travel
# * 0 = visa required for travel

# ### JuMP Modeling

# To model the problem as a mixed-integer linear program, we need a binary
# decision variable $x_c$ for each country $c$. $x_c$ is $1$ if we select
# passport $c$ and $0$ otherwise. Our objective is to minimize the sum
# $\sum x_c$ over all countries.

# Since we wish to visit all the countries, for every country, we must own at
# least one passport that lets us travel to that country visa free. For one
# destination, this can be mathematically represented as
# $\sum_{c \in C} a_{c,d} \cdot x_{d} \geq 1$, where $a$ is the `passport_data`
# dataframe.

# Thus, we can represent this problem using the following model:

# ```math
# \begin{aligned}
# \min && \sum_{c \in C} x_c \\
# \text{s.t.} && \sum_{c \in C} a_{c,d} x_c \geq 1 && \forall d \in C \\
# && x_c \in \{0,1\} && \forall c \in C.
# \end{aligned}
# ```

# We'll now solve the problem using JuMP:

using JuMP
import HiGHS

# First, create the set of countries:

C = passport_data.Passport

# Then, create the model and initialize the decision variables:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[C], Bin)
@objective(model, Min, sum(x))
@constraint(model, [d in C], passport_data[!, d]' * x >= 1)
model

# Now optimize!

optimize!(model)

# We can use the [`solution_summary`](@ref) function to get an overview of the
# solution:

solution_summary(model)

# ### Solution

# Let's  have a look at the solution in more detail:

println("Minimum number of passports needed: ", objective_value(model))

#-

println("Optimal passports:")
for c in C
    if value(x[c]) > 0.5
        println(" * ", c)
    end
end

# Interesting! We need some passports, like New Zealand and the United States,
# which have widespread access to a large number of countries. However, we also
# need passports like North Korea which only have visa-free access to a very
# limited number of countries.

# !!! note
#     We use `value(x[c]) > 0.5` rather than `value(x[c]) == 1` to avoid
#     excluding solutions like `x[c] = 0.99999` that are "1" to some tolerance.
