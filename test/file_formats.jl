#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestFileFormats

import JSON
using JuMP
using Test

function test_mof_file()
    model = Model()
    @variable(model, x)
    @constraint(model, my_c, 3 * x >= 1)
    @objective(model, Min, 2 * x^2 + x + 1)
    write_to_file(model, "my_model.mof.json")
    model_2 = read_from_file("my_model.mof.json")
    @test sprint(print, model) == sprint(print, model_2)
    return rm("my_model.mof.json")
end

function test_mof_io()
    model = Model()
    @variable(model, x)
    @constraint(model, my_c, 3 * x >= 1)
    @objective(model, Min, 2 * x^2 + x + 1)
    io = IOBuffer()
    @test_throws(
        ErrorException("Unable to infer the file format from an IO stream."),
        write(io, model; format = MOI.FileFormats.FORMAT_AUTOMATIC)
    )
    write(io, model; format = MOI.FileFormats.FORMAT_MOF)
    seekstart(io)
    @test_throws(
        ErrorException("Unable to infer the file format from an IO stream."),
        read(io, Model; format = MOI.FileFormats.FORMAT_AUTOMATIC)
    )
    seekstart(io)
    model_2 = read(io, Model; format = MOI.FileFormats.FORMAT_MOF)
    @test sprint(print, model) == sprint(print, model_2)
end

function test_mof_nlp()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    @NLconstraint(model, x^2 + y^2 <= 100.0)
    @constraint(model, x + y == 10)
    io = IOBuffer()
    write(io, model; format = MOI.FileFormats.FORMAT_MOF)
    seekstart(io)
    file = JSON.parse(io)
    @test file["name"] == "MathOptFormat Model"
    @test length(file["constraints"]) == 2
end

function test_nl_nlp()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    @NLconstraint(model, x^2 + y^2 <= 100.0)
    @constraint(model, x + y == 10)
    io = IOBuffer()
    write(io, model; format = MOI.FileFormats.FORMAT_NL)
    seekstart(io)
    file = read(io, String)
    # Check that the constant 100 occurs in the file
    @test occursin("n100", file)
    return
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if !startswith("$(name)", "test_")
            continue
        end
        @testset "$(name)" begin
            getfield(@__MODULE__, name)()
        end
    end
end

end

TestFileFormats.runtests()
