# How to use a custom binary

Many solvers are not written in Julia, but instead in languages like C or C++.
JuMP interacts with these solvers through binary dependencies.

For many open-source solvers, we automatically install the appropriate binary
when you run `Pkg.add("Solver")`. For example, `Pkg.add("ECOS")` will also
install the ECOS binary.

This page explains how this installation works, and how you can use a custom
binary.

!!! compat
    These instructions require Julia 1.6 or later.

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
Binaries are built using [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil)
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
In rare cases, however, you may require a custom binary. The two main reasons to
use a custom binary are:

 * You want a binary with custom compilation settings (for example, debugging)
 * You want a binary with a set of dependencies that are not available on
   Yggdrasil (for example, a commercial solver like Gurobi or CPLEX).

The following sections explain how to replace the binaries provided by a JLL
package with the custom ones you have compiled. As a reminder, we use ECOS as an
example for simplicity, but the steps are the same for other solvers.

## [Explore the JLL you want to override](@id jll_structure)

The first step is to explore the structure and filenames of the JLL package we
want to override.

Find the location of the files using `.artifact_dir`:
```julia
julia> using ECOS_jll

julia> ECOS_jll.artifact_dir
"/Users/oscar/.julia/artifacts/2addb75332eff5a1657b46bb6bf30d2410bc7ecf"
```

!!! tip
    This path may be different on other machines.

Here is what it contains:
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

Other solvers may have a `bin` directory containing executables. To use a custom
binary of ECOS, we need to replace `/lib/libecos.dylib` with our custom binary.

## [Compile a custom binary](@id compile_ecos)

The next step is to compile a custom binary. Because ECOS is written in C with
no dependencies, this is easy to do if you have a C compiler:
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

## Overriding a single library

To override the `libecos` library, we need to know what `ECOS_jll` calls it. (In
most cases, it will also be `libecos`, but not always.)

There are two ways you can check.

 1. Check the bottom of the JLL's GitHub README. For example,
    [ECOS_jll](https://github.com/JuliaBinaryWrappers/ECOS_jll.jl#products) has
    a single `LibraryProduct` called `libecos`.
 2. Type `ECOS_jll.` and the press the `[TAB]` key twice to auto-complete
    available options:
    ```julia
    julia> ECOS_jll.
    LIBPATH           PATH_list          best_wrapper       get_libecos_path   libecos_handle
    LIBPATH_list      __init__           dev_jll            is_available       libecos_path
    PATH              artifact_dir       find_artifact_dir  libecos
    ```
    Here you can see there is `libecos`, and more usefully for us,
    `libecos_path`.

Once you know the name of the variable to override (the one that ends in
`_path`), use  [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl)
to specify a new path:
```julia
using Preferences
set_preferences!(
    "LocalPreferences.toml",
    "ECOS_jll",
    "libecos_path" => "/tmp/jll_example/ecos/lib/libecos"
)
```

This will create a file in your current directory called `LocalPreferences.toml`
with the contents:
```julia
[ECOS_jll]
libecos_path = "/tmp/jll_example/ecos/lib/libecos"
```

Now if you restart Julia, you will see:
```julia
julia> using ECOS_jll

julia> ECOS_jll.libecos
"/tmp/jll_example/ecos/lib/libecos"
```

To go back to using the default library, just delete the `LocalPreferences.toml`
file.

## Overriding an entire artifact

Sometimes a solver may provide a number of libraries and executables, and
specifying the path for each of the becomes tedious. In this case, we can use
Julia's `Override.toml` to replace an entire artifact.

Overriding an entire artifact requires you to replicate the structure and
contents of the JLL package that we [explored above](@ref jll_structure).

In most cases you need only reproduce the `include`, `lib`, and `bin`
directories (if they exist). You can safely ignore any `logs` or `share`
directories. Take careful note of what files each directory contains and what
they are called.

For our ECOS example, we already reproduced the structure when we
[compiled ECOS](@ref compile_ecos).

So, now we need to tell Julia to use our custom installation instead of
the default. We can do this by making an override file at
`~/.julia/artifacts/Overrides.toml`.

`Overrides.toml` has the following content:
```julia
# Override for ECOS_jll
2addb75332eff5a1657b46bb6bf30d2410bc7ecf = "/tmp/jll_example/ecos"
```
where `2addb75332eff5a1657b46bb6bf30d2410bc7ecf` is the folder from the original
`ECOS_jll.artifact_dir` and `"/tmp/jll_example/ecos"` is the location of our new
installation. Replace these as appropriate for your system.

If you restart Julia after creating the override file, you will see:
```julia
julia> using ECOS_jll

julia> ECOS_jll.artifact_dir
"/tmp/jll_example/ecos"
```
Now when we use ECOS it will use our custom binary.

## Using Cbc with a custom binary

As a second example, we demonstrate how to use
[Cbc.jl](https://github.com/jump-dev/Cbc.jl) with a custom binary.

### Explore the JLL you want to override

First, let's check where `Cbc_jll` is installed:
```julia
julia> using Cbc_jll

julia> Cbc_jll.artifact_dir
"/Users/oscar/.julia/artifacts/e481bc81db5e229ba1f52b2b4bd57484204b1b06"

julia> readdir(Cbc_jll.artifact_dir)
5-element Vector{String}:
 "bin"
 "include"
 "lib"
 "logs"
 "share"

julia> readdir(joinpath(Cbc_jll.artifact_dir, "bin"))
1-element Vector{String}:
 "cbc"

julia> readdir(joinpath(Cbc_jll.artifact_dir, "lib"))
10-element Vector{String}:
 "libCbc.3.10.5.dylib"
 "libCbc.3.dylib"
 "libCbc.dylib"
 "libCbcSolver.3.10.5.dylib"
 "libCbcSolver.3.dylib"
 "libCbcSolver.dylib"
 "libOsiCbc.3.10.5.dylib"
 "libOsiCbc.3.dylib"
 "libOsiCbc.dylib"
 "pkgconfig"
```

### Compile a custom binary

Next, we need to compile Cbc. Cbc can be difficult to compile (it has a lot of
dependencies), but for macOS users there is a homebrew recipe:
```
(base) oscar@Oscars-MBP jll_example % brew install cbc
[ ... lines omitted ... ]
(base) oscar@Oscars-MBP jll_example % brew list cbc
/usr/local/Cellar/cbc/2.10.5/bin/cbc
/usr/local/Cellar/cbc/2.10.5/include/cbc/ (76 files)
/usr/local/Cellar/cbc/2.10.5/lib/libCbc.3.10.5.dylib
/usr/local/Cellar/cbc/2.10.5/lib/libCbcSolver.3.10.5.dylib
/usr/local/Cellar/cbc/2.10.5/lib/libOsiCbc.3.10.5.dylib
/usr/local/Cellar/cbc/2.10.5/lib/pkgconfig/ (2 files)
/usr/local/Cellar/cbc/2.10.5/lib/ (6 other files)
/usr/local/Cellar/cbc/2.10.5/share/cbc/ (59 files)
/usr/local/Cellar/cbc/2.10.5/share/coin/ (4 files)
```

### Override single libraries

To use `Preferences.jl` to override specific libraries we first check the names
of each library in `Cbc_jll`:
```julia
julia> Cbc_jll.
LIBPATH               cbc                    get_libcbcsolver_path  libOsiCbc_path
LIBPATH_list          cbc_path               is_available           libcbcsolver
PATH                  dev_jll                libCbc                 libcbcsolver_handle
PATH_list             find_artifact_dir      libCbc_handle          libcbcsolver_path
__init__              get_cbc_path           libCbc_path
artifact_dir          get_libCbc_path        libOsiCbc
best_wrapper          get_libOsiCbc_path     libOsiCbc_handle
```

Then we add the following to `LocalPreferences.toml`:
```julia
[Cbc_jll]
cbc_path = "/usr/local/Cellar/cbc/2.10.5/bin/cbc"
libCbc_path = "/usr/local/Cellar/cbc/2.10.5/lib/libCbc.3.10.5"
libOsiCbc_path = "/usr/local/Cellar/cbc/2.10.5/lib/libOsiCbc.3.10.5"
libcbcsolver_path = "/usr/local/Cellar/cbc/2.10.5/lib/libCbcSolver.3.10.5"
```

!!! info
    Note that capitalization matters, so `libcbcsolver_path` corresponds to
    `libCbcSolver.3.10.5`.

### Override entire artifact

To use the homebrew install as our custom binary we add the following to
`~/.julia/artifacts/Overrides.toml`:
```julia
# Override for Cbc_jll
e481bc81db5e229ba1f52b2b4bd57484204b1b06 = "/usr/local/Cellar/cbc/2.10.5"
```
