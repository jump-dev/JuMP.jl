# How to use a custom binary

Many solvers are not written in Julia, but instead in languages like C or C++.
For many open-source solvers, we automatically install the appropriate binary
when you run `Pkg.add("Solver")`. For example, `Pkg.add("ECOS")` will also
install the ECOS binary.

This page explains how this installation works, and how you can use a custom
binary.

## Background

Each solver that JuMP supports is structured as a Julia package. For example,
the interface for the [ECOS](https://github.com/embotech/ecos) solver is
provided by the [ECOS.jl](https://github.com/jump-dev/ECOS.jl) package.

!!! tip
    This page uses the example of ECOS.jl because it is simple to compile. Other
    solvers follow similar conventions. For example, the interface to the Clp
    solver is provided by Clp.jl.

The ECOS.jl package provides an interface between the C API of ECOS and
MathOptInterface. However, it does not handle the installation of the solver
binary; that is the job for a JLL package.

A JLL is a Julia package that wraps a pre-compiled binary.
The binaries are built using [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil)
(for example, [ECOS](https://github.com/JuliaPackaging/Yggdrasil/blob/master/E/ECOS/build_tarballs.jl))
and hosted in the [JuliaBinaryWrappers](https://github.com/JuliaBinaryWrappers)
GitHub repository (for example, [ECOS_jll.jl](https://github.com/JuliaBinaryWrappers/ECOS_jll.jl)).

JLL packages contain little code. Their only job is to `dlopen` a dynamic
library, along with any dependencies.

JLL packages manage their binary dependencies using [Julia's artifact system](https://pkgdocs.julialang.org/v1/artifacts/).
Each JLL package has an `Artifacts.toml` file which describes where to find each
binary artifact for each different platform that it might be installed on. Here
is the [Artifacts.toml file for ECOS_jll.jl](https://github.com/JuliaBinaryWrappers/ECOS_jll.jl/blob/main/Artifacts.toml).

The binaries installed by the JLL package should be sufficient for most users.
In a rare case however, you may require a custom binary. The two main reasons to
use a custom binary are:

 * You want a binary with custom compilation settings (for example, debugging)
 * You want a binary with a different set of dependencies that is available on
   Yggdrasil (for example, a commerial solver like Gurobi or CPLEX).

The following section explains how to replace the binaries provided by a JLL
package with the custom ones you have compiled. As a reminder, we use ECOS as an
example for simplicity, but the steps are the same for other solvers.

## Using ECOS with a custom binary

The first step is to find where Julia stores the binaries:
```julia
julia> using ECOS_jll

julia> ECOS_jll.artifact_dir
"/Users/oscar/.julia/artifacts/2addb75332eff5a1657b46bb6bf30d2410bc7ecf"
```
The name of the last folder is important, and we will need it later.

!!! tip
    This path may be different on other machines.

In order to use a custom installation of ECOS, we need to reproduce the
structure of this directory. Here is what it contains:
```julia
julia> readdir(ECOS_jll.artifact_dir)
4-element Vector{String}:
 "include"
 "lib"
 "logs"
 "share"

julia> readdir(joinpath(ECOS_jll.artifact_dir, "lib"))
1-element Vector{String}:
 "libecos.dylib"
```

In most cases you need only reproduce the `include`, `lib`, and `bin`
directories (if they exist). You can safely ignore any `logs` or `share`
directories. Take careful note of what files each directory contains and what
they are called.

### Compile a custom version of the solver

The next step is to compile a custom version of ECOS. Because ECOS is written in
C with no dependencies, this is easy to do if you have a C compiler:
```julia
oscar@Oscars-MBP jll_example % git clone https://github.com/embotech/ecos.git
[... lines omitted ...]
oscar@Oscars-MBP jll_example % cd ecos
oscar@Oscars-MBP ecos % make shared
[... many lines omitted...]
oscar@Oscars-MBP ecos % mkdir lib
oscar@Oscars-MBP ecos % cp libecos.dylib lib
```

!!! warning
    Compiling custom solver binaries is an advanced operation. Due to the
    complexities of compiling various solvers, the JuMP community is unable to
    help you diagnose and fix compilation issues.

After this compilation step, we now have a folder `/tmp/jll_example/ecos`
that contains `lib` and `include` directories with the same files as `ECOS_jll`:
```julia
julia> readdir(joinpath("ecos", "lib"))
1-element Vector{String}:
 "libecos.dylib"
```

### Override the artifact location

As a final step, we need to tell Julia to use our custom installation instead of
the default. We can do this by making an over-ride file at
`~/.julia/artifacts/Overrides.toml`.

`Overrides.toml` has the following content:
```julia
# Override for ECOS_jll
2addb75332eff5a1657b46bb6bf30d2410bc7ecf = "/tmp/jll_example/ecos"
```
Where `2addb75332eff5a1657b46bb6bf30d2410bc7ecf` is the folder from the original
`ECOS_jll.artifact_dir` and `"/tmp/jll_example/ecos"` is the location of our new
installation. Replace these as appropriate for your system.

If you restart Julia after creating the override file, you will see:
```julia
julia> using ECOS_jll

julia> ECOS_jll.artifact_dir
"/tmp/jll_example/ecos"
```
Now when we use ECOS it will use our custom binary.
