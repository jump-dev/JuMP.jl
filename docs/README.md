JuMP.jl Documentation README
================================

JuMP.jl's documentation is written with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). To install it, run the following command in a Julia session:

```julia
Pkg.add("Documenter")
```


Building the documentation
--------------------------

The documentation is built using the following command:

```julia
julia --project=. --color=yes make.jl
```

The compiled documents can be viewed at `build/index.html`.
