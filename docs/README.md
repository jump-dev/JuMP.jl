JuMP.jl Documentation README
================================

JuMP.jl's documentation is written with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). To install it, run the following command in a Julia session:

```julia
Pkg.add("Documenter")
```


Building the documentation
--------------------------

**This documentation currently requires Julia 1.0 to build.** 

If you are hacking on JuMP, and want to build the documentation, do the following:

```julia

(base) Oscars-MBP:JuMP oscar$ julia --project=docs
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.0.5 (2019-09-09)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(docs) pkg> st
    Status `~/.julia/dev/JuMP/docs/Project.toml`
  [e30172f5] Documenter v0.23.4

julia> pwd()
"/Users/oscar/.julia/dev/JuMP"

(docs) pkg> dev .
 Resolving package versions...
  Updating `~/.julia/dev/JuMP/docs/Project.toml`
  [4076af6c] + JuMP v0.21.3 [`..`]
  Updating `~/.julia/dev/JuMP/docs/Manifest.toml`
  [6e4b80f9] + BenchmarkTools v0.5.0
  [9e28174c] + BinDeps v1.0.1
  [b99e7846] + BinaryProvider v0.5.10
  [49dc2e85] + Calculus v0.5.1
  [523fee87] + CodecBzip2 v0.6.0
  [944b1d66] + CodecZlib v0.6.0
  [bbf7d656] + CommonSubexpressions v0.3.0
  [864edb3b] + DataStructures v0.17.19
  [163ba53b] + DiffResults v1.0.2
  [b552c78f] + DiffRules v1.0.1
  [f6369f11] + ForwardDiff v0.10.12
  [cd3eb016] + HTTP v0.8.16
  [83e8ac13] + IniFile v0.5.0
  [7d188eb4] + JSONSchema v0.3.0
  [4076af6c] + JuMP v0.21.3 [`..`]
  [1914dd2f] + MacroTools v0.5.5
  [b8f27783] + MathOptInterface v0.9.14
  [739be429] + MbedTLS v0.6.8
  [d8a4904e] + MutableArithmetics v0.2.10
  [77ba4419] + NaNMath v0.3.3
  [bac558e1] + OrderedCollections v1.3.0
  [276daf66] + SpecialFunctions v0.8.0
  [90137ffa] + StaticArrays v0.12.3
  [3bb67fe8] + TranscodingStreams v0.9.5
  [30578b45] + URIParser v0.4.1
  [a5390f91] + ZipFile v0.8.4
  [2f01184e] + SparseArrays 
  [10745b16] + Statistics 

julia> include("docs/make.jl")
[ Info: SetupBuildDirectory: setting up build directory.
[ Info: Doctest: running doctests.
[ Info: ExpandTemplates: expanding markdown templates.
[ Info: CrossReferences: building cross-references.
[ Info: CheckDocument: running document checks.
[ Info: Populate: populating indices.
[ Info: RenderDocument: rendering document.
[ Info: HTMLWriter: rendering HTML pages.
┌ Info: Deployment criteria:
│ - ✔ ENV["TRAVIS_REPO_SLUG"]="" occurs in repo="github.com/jump-dev/JuMP.jl.git"
│ - ✘ ENV["TRAVIS_PULL_REQUEST"]="" is "false"
│ - ✔ ENV["TRAVIS_TAG"]="" is (i) empty or (ii) a valid VersionNumber
│ - ✘ ENV["TRAVIS_BRANCH"]="" matches devbranch="master" (if tag is empty)
│ - ✘ ENV["DOCUMENTER_KEY"] exists
│ - ✔ ENV["TRAVIS_EVENT_TYPE"]="" is not "cron"
└ Deploying: ✘

```

Note that by adding JuMP to the environment, we have modified `/docs/Project.toml`, which 
we _don't_ want to commit in our PR to JuMP. Therefore, before commiting, run:
```julia
(docs) pkg> rm JuMP
  Updating `~/.julia/dev/JuMP/docs/Project.toml`
  [4076af6c] - JuMP v0.21.3 [`..`]
  Updating `~/.julia/dev/JuMP/docs/Manifest.toml`
  [6e4b80f9] - BenchmarkTools v0.5.0
  [9e28174c] - BinDeps v1.0.1
  [b99e7846] - BinaryProvider v0.5.10
  [49dc2e85] - Calculus v0.5.1
  [523fee87] - CodecBzip2 v0.6.0
  [944b1d66] - CodecZlib v0.6.0
  [bbf7d656] - CommonSubexpressions v0.3.0
  [864edb3b] - DataStructures v0.17.19
  [163ba53b] - DiffResults v1.0.2
  [b552c78f] - DiffRules v1.0.1
  [f6369f11] - ForwardDiff v0.10.12
  [cd3eb016] - HTTP v0.8.16
  [83e8ac13] - IniFile v0.5.0
  [7d188eb4] - JSONSchema v0.3.0
  [4076af6c] - JuMP v0.21.3 [`..`]
  [1914dd2f] - MacroTools v0.5.5
  [b8f27783] - MathOptInterface v0.9.14
  [739be429] - MbedTLS v0.6.8
  [d8a4904e] - MutableArithmetics v0.2.10
  [77ba4419] - NaNMath v0.3.3
  [bac558e1] - OrderedCollections v1.3.0
  [276daf66] - SpecialFunctions v0.8.0
  [90137ffa] - StaticArrays v0.12.3
  [3bb67fe8] - TranscodingStreams v0.9.5
  [30578b45] - URIParser v0.4.1
  [a5390f91] - ZipFile v0.8.4
  [2f01184e] - SparseArrays 
  [10745b16] - Statistics 
```

The compiled documents can be viewed at `build/index.html`.
