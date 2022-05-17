# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Serving web apps

# This tutorial demonstrates how to setup and serve JuMP models via a REST API.

# First, we need the standard JuMP packages:

using JuMP
import HiGHS

# We also need [HTTP.jl](https://github.com/JuliaWeb/HTTP.jl) to act as our
# REST server, and [JSON.jl](https://github.com/JuliaIO/JSON.jl) to marshall
# data.

import HTTP
import JSON

# ## The server side

# The core components of our REST server are _endpoints_. These are functions
# which accept a `Dict{String,Any}` of input parameters, and return a
# `Dict{String,Any}` as output. The types are `Dict{String,Any}` because we're
# going to read these to and from JSON.

# Here's a very simple endpoint: it accepts `params` as input, formulates and
# solves a trivial mixed-integer program, and then returns a dictionary with the
# result.

function endpoint_solve(params::Dict{String,Any})
    if !haskey(params, "lower_bound")
        return Dict{String,Any}(
          "status" => "failure",
          "reason" => "missing lower_bound param",
        )
    end
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x >= params["lower_bound"], Int)
    optimize!(model)
    return Dict{String,Any}(
        "status" => "okay",
        "terminaton_status" => termination_status(model),
        "primal_status" => primal_status(model),
        "x" => value(x),
    )
end

# When we call this, we get:

endpoint_solve(Dict{String,Any}("lower_bound" => 1.2))

#-

endpoint_solve(Dict{String,Any}())

# For a second function, we need a function that accepts an `HTTP.Request`
# object and returns an `HTTP.Reponse` object.

function serve_solve(request::HTTP.Request)
    data = JSON.parse(String(request.body))
    solution = endpoint_solve(data)
    return HTTP.Response(200, JSON.json(solution))
end

# Finally, we need an HTTP server:

router = HTTP.Router()
HTTP.@register(router, "/solve", serve_solve)
@async HTTP.serve(router, HTTP.ip"127.0.0.1", 8080)

# !!! info
#     `@async` run the server in a background process. If you omit `@async`, the
#     server will block the current Julia process. The reason we've used
#     `@async` is so that we can demonstrate the client-side code in this
#     tutorial without starting multiple instances of Julia.

# ## The client side

# Now that we have a server, we can send it requests via this function:

function send_request(data::Dict)
    ret = HTTP.request(
        "POST",
        ## This should match the URL and endpoint we defined for our server.
        "http://127.0.0.1:8080/solve",
        ["Content-Type" => "application/json"],
        JSON.json(data),
    )
    if ret.status != 200
        ## This could happen if there are time-outs, network errors, etc.
        return Dict(
            "status" => "failure",
            "reason" => "HTTP status = $(ret.status)",
        )
    end
    return JSON.parse(String(ret.body))
end

# Let's see what happens:

send_request(Dict("lower_bound" => 0))

#-

send_request(Dict("lower_bound" => 1.2))

#-

send_request(Dict("invalid_param" => 1.2))

# ## Next steps

# This tutorial has been deliberately kept simple. For more complicated
# examples, consult the [HTTP.jl documentation](https://juliaweb.github.io/HTTP.jl/stable/).
