# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Serving web apps

# This tutorial demonstrates how to setup and serve JuMP models via a REST API.

# In the example app we are building, we solve a trivial mixed-integer program,
# which is parameterized by the lower bound of a variable. To call the service,
# users send an HTTP POST request with JSON contents indicating the lower bound.
# The returned value is the solution of the mixed-integer program as JSON.

# First, we need JuMP and a solver:

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
    elseif !(params["lower_bound"] isa Real)
        return Dict{String,Any}(
            "status" => "failure",
            "reason" => "lower_bound is not a number",
        )
    end
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x >= params["lower_bound"], Int)
    optimize!(model)
    ret = Dict{String,Any}(
        "status" => "okay",
        "terminaton_status" => termination_status(model),
        "primal_status" => primal_status(model),
    )
    ## Only include the `x` key if it has a value.
    if has_values(model)
        ret["x"] = value(x)
    end
    return ret
end

# When we call this, we get:

endpoint_solve(Dict{String,Any}("lower_bound" => 1.2))

#-

endpoint_solve(Dict{String,Any}())

# For a second function, we need a function that accepts an `HTTP.Request`
# object and returns an `HTTP.Response` object.

function serve_solve(request::HTTP.Request)
    data = JSON.parse(String(request.body))
    solution = endpoint_solve(data)
    return HTTP.Response(200, JSON.json(solution))
end

# Finally, we need an HTTP server. There are a variety of ways you can do this
# in HTTP.jl. We use an explicit `Sockets.listen` so we have manual control of
# when we shutdown the server.

function setup_server(host, port)
    server = HTTP.Sockets.listen(host, port)
    HTTP.serve!(host, port; server = server) do request
        try
            ## Extend the server by adding other endpoints here.
            if request.target == "/api/solve"
                return serve_solve(request)
            else
                return HTTP.Response(404, "target $(request.target) not found")
            end
        catch err
            ## Log details about the exception server-side
            @info "Unhandled exception: $err"
            ## Return a response to the client
            return HTTP.Response(500, "internal error")
        end
    end
    return server
end

# !!! warning
#     HTTP.jl does not serve requests on a separate thread. Therefore, a
#     long-running job will block the main thread, preventing concurrent users from
#     submitting requests. To work-around this, read [HTTP.jl issue 798](https://github.com/JuliaWeb/HTTP.jl/issues/798)
#     or watch [Building Microservices and Applications in Julia](https://www.youtube.com/watch?v=uLhXgt_gKJc&t=9543s)
#     from JuliaCon 2020.

server = setup_server(HTTP.ip"127.0.0.1", 8080)

# ## The client side

# Now that we have a server, we can send it requests via this function:

function send_request(data::Dict; endpoint::String = "solve")
    ret = HTTP.request(
        "POST",
        ## This should match the URL and endpoint we defined for our server.
        "http://127.0.0.1:8080/api/$endpoint",
        ["Content-Type" => "application/json"],
        JSON.json(data),
    )
    if ret.status != 200
        ## This could happen if there are time-outs, network errors, etc.
        return Dict(
            "status" => "failure",
            "code" => ret.status,
            "body" => String(ret.body),
        )
    end
    return JSON.parse(String(ret.body))
end

# Let's see what happens:

send_request(Dict("lower_bound" => 0))

#-

send_request(Dict("lower_bound" => 1.2))

# If we don't send a `lower_bound`, we get:

send_request(Dict("invalid_param" => 1.2))

# If we don't send a `lower_bound` that is a number, we get:

send_request(Dict("lower_bound" => "1.2"))

# Finally, we can shutdown our HTTP server:

close(server)

# ## Next steps

# For more complicated examples relating to HTTP servers, consult the
# [HTTP.jl documentation](https://juliaweb.github.io/HTTP.jl/stable/).

# To see how you can integrate this with a larger JuMP model, read
# [Design patterns for larger models](@ref).
