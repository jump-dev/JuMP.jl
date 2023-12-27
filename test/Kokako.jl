# Copyright 2023, Oscar Dowson
#
# Use of this source code is governed by an MIT-style license that can be found
# at https://opensource.org/licenses/MIT.

# !!! note "Info for JuMP developers"
#     If you remove this file from JuMP.jl, be sure to update the documentation
#     in `/docs/src/developers/extensions.md`.

module Kokako

import Test

"""
    run_tests(
        module_to_test::Module,
        args...;
        test_prefix::String = "test_",
        include_names::Vector{String} = String[""],
        exclude_names::Vector{String} = String[],
        kwargs...
    )

Test a single module `module_to_test` by looping through every function
whose function name begins with `test_prefix`, contains a fragment from the
vector `include_names` and does not contain a fragment in the vector
`exclude_names`.

Each function meeting the description is called in a new `Test.@testset`,
passing `args` and `kwargs` as positional and keyword arguments.
"""
function run_tests(
    module_to_test::Module,
    args...;
    test_prefix::String = "test_",
    include_names::Vector{String} = String[""],
    exclude_names::Vector{String} = String[],
    kwargs...,
)
    for name in names(module_to_test; all = true)
        test_function, test_name = getfield(module_to_test, name), string(name)
        if !(test_function isa Function)
            continue
        elseif !startswith(test_name, test_prefix)
            continue
        elseif !any(needle -> occursin(needle, test_name), include_names)
            continue
        elseif any(needle -> occursin(needle, test_name), exclude_names)
            continue
        end
        Test.@testset "$name" begin
            test_function(args...; kwargs...)
        end
    end
    return
end

"""
    run_tests(
        modules_to_test::Vector{Pair{String,Module}},
        args...;
        test_prefix::String = "test"_,
        include_names::Dict{String,Vector{String}} = Dict{String,Vector{String}}(),
        exclude_names::Dict{String,Vector{String}} = Dict{String,Vector{String}}(),
        kwargs...,
    )

Test all the modules in `modules_to_test`. See the singular `run_tests` function
for a description of the arguments.

`include_names` and `exclude_names` are dictionaries which map the first
element in the `Pair{String,Module}` pair to the value that is passed to
`run_tests` for that module.
"""
function run_tests(
    modules_to_test::Vector{Pair{String,Module}},
    args...;
    test_prefix::String = "test_",
    include_names::Dict{String,Vector{String}} = Dict{String,Vector{String}}(),
    exclude_names::Dict{String,Vector{String}} = Dict{String,Vector{String}}(),
    kwargs...,
)
    for (file, module_to_test) in modules_to_test
        Test.@testset "$file:$test_prefix" begin
            run_tests(
                module_to_test,
                args...;
                test_prefix = test_prefix,
                include_names = get(include_names, file, String[""]),
                exclude_names = get(exclude_names, file, String[]),
                kwargs...,
            )
        end
    end
    return
end

"""
    include_modules_to_test(dir::String[, files::Vector{String}])

Find and include all modules in Julia source files beginning with `test_` that
are contained within `dir` (including subdirectories), or in the explicit list
of filenames `files` that are relative to `dir`.

Returns a `Vector{Pair{String,Module}}` mapping the filename to the included
module that can be passed to `run_tests`.

!!! warning
    Because this function calls `Base.include`, this function _must_ be called
    from the top-level of the runtests.jl script. It _must not_ be called from
    within a function or other local scope.
"""
function include_modules_to_test(dir::String)
    if !isempty(ARGS)
        return include_modules_to_test(dir, ARGS)
    end
    test_files = String[]
    for (root, dirs, files) in walkdir(dir)
        for file in files
            if startswith(file, "test_") && endswith(file, ".jl")
                filename = joinpath(root, file)
                filename = replace(filename, joinpath(dir, "") => "")
                push!(test_files, filename)
            end
        end
    end
    sort!(test_files)
    return include_modules_to_test(dir, test_files)
end

function include_modules_to_test(dir::String, files::Vector{String})
    modules = Pair{String,Module}[]
    for file in files
        if isfile(joinpath(dir, file))
            @info "Loading $file"
            push!(modules, file => Base.include(Main, joinpath(dir, file)))
        end
    end
    return modules
end

"""
    include_modules_to_test(package::Module)

Calls `include_modules_to_test(dir)` where `dir` is the `/test` directory of the
package `package`.
"""
function include_modules_to_test(package::Module)
    root = dirname(dirname(pathof(package)))
    modules = include_modules_to_test(joinpath(root, "test"))
    if VERSION >= v"1.9"
        append!(modules, include_modules_to_test(joinpath(root, "ext")))
    end
    return modules
end

end  # module
