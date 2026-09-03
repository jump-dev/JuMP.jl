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

# # Working with Data Files

# **Originally Contributed by**: Arpit Bhatia

# In many cases we might need to read data available in an external file rather
# than type it into Julia ourselves.

# This tutorial is concerned with reading tabular data into Julia and using it
# for a JuMP model.

# We'll be reading data using the [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl)
# package and some other packages specific to file types.

# The data are stored in the [`/docs/src/tutorials/Getting started/data`](https://github.com/jump-dev/JuMP.jl/tree/master/docs/src//tutorials/Getting started/data)
# directory of the JuMP source code.

const DATA_DIR = joinpath(@__DIR__, "data");

# !!! note
#     There are multiple ways to read the same kind of data into Julia. This
#     tutorial focuses on DataFrames.jl because it provides the ecosystem to
#     work with most of the required file types in a straightforward manner.

# ### DataFrames.jl

# The `DataFrames` package provides a set of tools for working with tabular
# data. It is available through the Julia package system.
# ```julia
# using Pkg
# Pkg.add("DataFrames")
# ```

import DataFrames

# ### What is a DataFrame?

# A DataFrame is a data structure like a table or spreadsheet. You can use it
# for storing and exploring a set of related data values. Think of it as a
# smarter array for holding tabular data.

# ## Reading Tabular Data into a DataFrame

# We will begin by reading data from different file formats into a DataFrame
# object.

# ### Excel Sheets

# Excel files can be read using the [XLSX.jl](https://github.com/filepenoris/XLSX.jl)
# package.

# ```julia
# Pkg.add("XLSX")
# ```

import XLSX

# To read a Excel file into a DataFrame, we use the following julia code. The
# first argument to the `readtable` function is the file to be read and the
# second argument is the name of the sheet.

excel_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "SalesData.xlsx"), "SalesOrders")...
)

# ### CSV files

# CSV and other delimited text files can be read by the CSV.jl package.

# ```julia
# Pkg.add("CSV")
# ```

import CSV

# To read a CSV file into a DataFrame, we use the `CSV.read` function.

csv_df = CSV.read(joinpath(DATA_DIR, "StarWars.csv"), DataFrames.DataFrame)

# ### Other Delimited Files

# We can also use the `CSV.jl` package to read any other delimited text file
# format.

# By default, CSV.File will try to detect a file's delimiter from the first 10
# lines of the file.

# Candidate delimiters include `','`, `'\t'`, `' '`, `'|'`, `';'`, and `':'`. If
# it can't auto-detect the delimiter, it will assume `','`.

# Let's take the example of space separated data.

ss_df = CSV.read(joinpath(DATA_DIR, "Cereal.txt"), DataFrames.DataFrame)

# We can also specify the delimiter by passing the `delim` argument.

delim_df = CSV.read(
    joinpath(DATA_DIR, "Soccer.txt"), DataFrames.DataFrame, delim = "::"
)

# Note that by default, are read-only. If we wish to make changes to the data
# read, we pass the `copycols = true` argument in the function call.

ss_df = CSV.read(
    joinpath(DATA_DIR, "Cereal.txt"), DataFrames.DataFrame, copycols = true
)

# ## Working with DataFrames

# Now that we have read the required data into a DataFrame, let us look at some
# basic operations we can perform on it.

# ### Querying Basic Information

# The `size` function gets us the dimensions of the DataFrame.

DataFrames.size(ss_df)

# We can also us the `nrow` and `ncol` functions to get the number of rows and
# columns respectively.

DataFrames.nrow(ss_df), DataFrames.ncol(ss_df)

# The `describe` function gives basic summary statistics of data in a DataFrame.

DataFrames.describe(ss_df)

# Names of every column can be obtained by the `names` function.

DataFrames.names(ss_df)

# Corresponding data types are obtained using the broadcasted `eltype` function.

eltype.(ss_df)

# ### Accessing the Data

# Similar to regular arrays, we use numerical indexing to access elements of a
# DataFrame.

csv_df[1, 1]

# The following are different ways to access a column.

csv_df[!, 1]

#-

csv_df[!, :Name]

#-

csv_df.Name

#-

csv_df[:, 1] # Note that this creates a copy.

# The following are different ways to access a row.

csv_df[1:1, :]

#-

csv_df[1, :] # This produces a DataFrameRow.

# We can change the values just as we normally assign values.

# Assign a range to scalar.

excel_df[1:3, 5] .= 1

# Vector to equal length vector.

excel_df[4:6, 5] = [4, 5, 6]

# Subset of the DataFrame to another data frame of matching size.

excel_df[1:2, 6:7] =  DataFrames.DataFrame(
    [-2 -2; -2 -2], [Symbol("Unit Cost"), :Total]
)

#-

excel_df

# !!! tip
#     There are a lot more things which can be done with a DataFrame. Read the
#     [docs](https://juliadata.github.io/DataFrames.jl/stable/) for more
#     information.

# ## A Complete Modelling Example - Passport Problem

# Let's now apply what we have learnt to solve a real modelling problem.

# The [Passport Index Dataset](https://github.com/ilyankou/passport-index-dataset)
# lists travel visa requirements for 199 countries, in `.csv` format. Our task
# is to find out the minimum number of passports required to visit all
# countries.

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

# Thus, the values we are interested in are -1 and 3. Let us modify the data in
# the following manner:

passport_data = CSV.read(
    joinpath(DATA_DIR, "passport-index-matrix.csv"),
    DataFrames.DataFrame;
    copycols = true,
)

for i in 1:DataFrames.nrow(passport_data)
    for j in 2:DataFrames.ncol(passport_data)
        if passport_data[i, j] == -1 || passport_data[i, j] == 3
            passport_data[i, j] = 1
        else
            passport_data[i, j] = 0
        end
    end
end

# The values in the cells now represent:
# * 1 = no visa required for travel
# * 0 = visa required for travel

# Let us associate each passport with a decision variable $pass_{cntr}$ for
# each country. We want to minimize the sum $\sum pass_{cntr}$ over all countries.

# Since we wish to visit all the countries, for every country, we should own at
# least one passport that lets us travel to that country visa free. For one
# destination, this can be mathematically represented as
# $\sum_{cntr \in world} passportdata_{cntr,dest} \cdot pass_{cntr} \geq 1$.

# Thus, we can represent this problem using the following model:

# ```math
# \begin{aligned}
# \min && \sum_{cntr \in World} pass_{cntr} \\
# \text{s.t.} && \sum_{cntr \in World} passportdata_{cntr,dest} \cdot pass_{cntr} \geq 1 && \forall dest \in World \\
# && pass_{cntr} \in \{0,1\} && \forall cntr \in World
# \end{aligned}
# ```

# We'll now solve the problem using JuMP.

using JuMP
import GLPK

# First, create the set of countries:

World = names(passport_data)[2:end]

# Then, create the model and initialize the decision variables:
model = Model(GLPK.Optimizer)
@variable(model, pass[cntr in World], Bin)

# Define the objective function

@objective(model, Min, sum(pass[cntr] for cntr in World))

#-

@constraint(model, [dest in World], passport_data[:, dest]' * pass >= 1)

# Now optimize!

optimize!(model)
println("Minimum number of passports needed: ", objective_value(model))

#-

optimal_passports = [cntr for cntr in World if value(pass[cntr]) > 0.5]
println("Countries:")
for p in optimal_passports
    println(" ", p)
end

# !!! note
#     We use `value(pass[i]) > 0.5` rather than `value(pass[i]) == 1` to avoid
#     excluding solutions like `pass[i] = 0.99999` that are "1" to some
#     tolerance.
