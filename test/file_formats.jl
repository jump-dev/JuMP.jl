#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using JuMP
using Test

@testset "File formats" begin
    @testset "MOF" begin
        model = Model()
        @variable(model, x)
        @constraint(model, my_c, 3 * x >= 1)
        @objective(model, Min, 2 * x^2 + x + 1)
        write_to_file(model, "my_model.mof.json")
        model_2 = read_from_file("my_model.mof.json")
        @test sprint(print, model) == sprint(print, model_2)
        rm("my_model.mof.json")
    end
    @testset "MPS" begin
        model = Model()
        @variable(model, x >= 0)
        @constraint(model, my_c, 3 * x >= 1)
        @objective(model, Min, 2 * x)
        write_to_file(model, "my_model.mps")
        model_2 = read_from_file("my_model.mps")
        @test sprint(print, model) == sprint(print, model_2)
        rm("my_model.mps")
    end
    @testset "LP" begin
        model = Model()
        @variable(model, x >= 0)
        @constraint(model, my_c, 3 * x >= 1)
        @objective(model, Min, 2 * x)
        write_to_file(model, "my_model.lp")
        @test read("my_model.lp", String) ==
            "minimize\nobj: 2 x\nsubject to\nmy_c: 3 x >= 1\nBounds\nx >= 0\nEnd\n"
        @test_throws(
            ErrorException("read! is not implemented for LP files."),
            read_from_file("my_model.lp")
        )
        rm("my_model.lp")
    end
    @testset "CBF" begin
        model = Model()
        @variable(model, X[1:2, 1:2], PSD)
        @constraint(model, my_c, sum(X) >= 1)
        @objective(model, Min, sum(X))
        write_to_file(model, "my_model.cbf")
        @test read("my_model.cbf", String) ==
            "VER\n3\n\nOBJSENSE\nMIN\n\nVAR\n3 1\nF 3\n\nOBJACOORD\n3\n0 1.0\n1 2.0\n2 1.0\n\nCON\n1 1\nL+ 1\n\nACOORD\n3\n0 0 1.0\n0 1 2.0\n0 2 1.0\n\nBCOORD\n1\n0 -1.0\n\nPSDCON\n1\n2\n\nHCOORD\n3\n0 0 0 0 1.0\n0 1 1 0 1.0\n0 2 1 1 1.0\n\n"
        model_2 = read_from_file("my_model.cbf")
        # Note: we replace ' in ' => ' ∈ ' because the unicode doesn't print on
        # Windows systems for some reason.
        @test replace(sprint(print, model_2), " in " => " ∈ ") ==
            "Min noname + 2 noname + noname\nSubject to\n [noname + 2 noname + noname - 1] ∈ MathOptInterface.Nonnegatives(1)\n [noname, noname, noname] ∈ MathOptInterface.PositiveSemidefiniteConeTriangle(2)\n"
        rm("my_model.cbf")
    end
    @testset "Base read/write via io" begin
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
end
