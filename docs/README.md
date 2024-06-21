# JuMP Documentation README

**The documentation currently requires Julia 1.9 to build.**

JuMP's documentation is written with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).

## Initial setup

To build the documentation, you need to do a series of initialization steps.
However, you only need to do this once.

First, you will need a local copy of JuMP. If you don't have one already, run:
```
$ julia -e 'import Pkg; Pkg.develop("JuMP")'
```

This will create a copy of JuMP at `~/.julia/dev/JuMP`. (On Windows, this will
be located at `C:\\Users\\<your_user_name>\\.julia\\dev\\JuMP`.) Open a terminal,
and `cd` to that directory:
```
$ cd ~/.julia/dev/JuMP
```

The next step is to setup the `docs` environment.
```
$ julia --project=docs -e 'import Pkg; Pkg.instantiate(); Pkg.develop(Pkg.PackageSpec(path=".")); Pkg.develop("MathOptInterface")'
```

Now you're ready to build the documentation.

## Building the docs

Build the docs as follows:
```
$ cd ~/.julia/dev/JuMP
$ julia --project=docs docs/make.jl
```

The compiled documents can be viewed at `~/.julia/dev/JuMP/docs/build/index.html`.

## Updating Project.toml

Project.toml fixes the versions of JuMP-related packages such as MOI and the
solvers. This is to prevent minor changes in the upstream solvers (for example,
numerical differences or changes to their raw solver statuses) from causing the
documentation builds to fail on an un-related PR.

These versions should be periodically updated.

## Merging pull requests

If you merge two pull requests which modify the documentation in close succession,
it may happen that the second pull request finishes CI before the first. Therefore,
the docs from the second pull request will deploy and be updated in `/dev`. However,
when the first pull request finishes and deploys the docs, these changes will be
overwritten because they were not in the repository when the first pull request was
merged.

For an example of this, see https://github.com/jump-dev/JuMP.jl/issues/2947.

The work-around is to ensure that you never merge two pull requests which modify
the documentation in close succession---wait 10 minutes between merges.
