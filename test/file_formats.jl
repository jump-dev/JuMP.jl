#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

@testset "File formats" begin
    @testset "MPS" begin
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
            ErrorException("Read from file is not implemented for LP files."),
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
        @test sprint(print, model_2) ==
            "Min noname + 2 noname + noname\nSubject to\n [noname + 2 noname + noname - 1] ∈ MathOptInterface.Nonnegatives(1)\n [noname, noname, noname] ∈ MathOptInterface.PositiveSemidefiniteConeTriangle(2)\n"
        rm("my_model.cbf")
    end
end
