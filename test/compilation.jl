#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/compilation.jl
# Testing ahead-of-time compilation of JuMP using PackageCompiler.
# Intended to be run separately from runtests.jl.
#############################################################################

using Compat.Test, PackageCompiler

@testset "Compilation" begin
    @testset "Snooping simulation" begin
        #= Executes snooping commands as they will be executed in PackageCompiler, but captures any errors. Errors are elided
            when running PackageCompiler itself. Snooping commands demonstrate illustrative uses of JuMP (LP and MILP creation
            and solving). This code is adapted from PackageCompiler.snoop_vanilla(). =#

            snoop_commands_file = joinpath(@__DIR__, "snoop_commands.jl")

            if VERSION ≤ v"0.7-"
                code_object = """
                while !eof(STDIN)
                    eval(Main, deserialize(STDIN))
                end
                """
                julia_cmd = PackageCompiler.build_julia_cmd(
                    PackageCompiler.get_backup!(false, nothing), nothing, nothing, nothing, nothing, nothing,
                    nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
                )
                in, io = open(`$julia_cmd --eval $code_object`, "w", STDOUT)

                serialize(in, quote
                    try
                        include($snoop_commands_file)
                    catch ex
                        exit(100)
                    end
                    exit(0)
                end)
                close(in)
                wait(io)

                @test io.exitcode == 0
            else
                using Serialization

                code_object = """
                        using Serialization
                        while !eof(stdin)
                            Core.eval(Main, deserialize(stdin))
                        end
                        """
                julia_cmd = PackageCompiler.build_julia_cmd(
                    PackageCompiler.get_backup!(false, nothing), nothing, nothing, nothing, nothing, nothing,
                    nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
                )
                process = open(Cmd(`$julia_cmd --eval $code_object`, dir=dirname(@__FILE__)), stdout, write=true)
                serialize(process, quote
                    try
                        include($snoop_commands_file)
                    catch
                        exit(100)
                    end
                    exit(0)
                end)

                wait(process)

                @test process.exitcode == 0
            end
    end  # @testset "Snooping simulation"

    @testset "Compilation of system image" begin
        if VERSION > v"0.7-"
            using Libdl
        end

        # Generate a new system image including JuMP, but do not overwrite current system image for Julia installation.
        PackageCompiler.compile_package(("JuMP", joinpath(dirname(@__FILE__), "snoop_commands.jl")))

        sys_image_file = joinpath(PackageCompiler.sysimg_folder(), "sys.$(Libdl.dlext)")
        julia_exe = Base.julia_cmd().exec[1]

        @test isfile(sys_image_file)
        @test success(`$julia_exe -J $sys_image_file`)
    end  # @testset "Compilation of system image"

end  # @testset "Compilation"
