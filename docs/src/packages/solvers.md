# Introduction

This section of the documentation contains brief documentation for some of the
solvers that JuMP supports. The list of solvers is not exhaustive, but instead
is intended to help you discover commonly used solvers.

## Affiliation

Packages beginning with `jump-dev/` are developed and maintained by the
JuMP developers. In many cases, these packages wrap external solvers that are
not developed by the JuMP developers and, while the Julia packages are all
open-source, in some cases the solvers themselves are closed source commercial
products.

Packages that do not begin with `jump-dev/` are developed independently, and the
JuMP developers have no affiliation with their developers or control over their
contents. The README files from each solver are taken directly from the GitHub
repository of each package. They are included here for the benefit of users
reading a single unified set of documentation.

## Adding new solvers

Written a solver? Add it to this section of the JuMP documentation by making
a pull request to the [`docs/packages.toml`](https://github.com/jump-dev/JuMP.jl/blob/master/docs/packages.toml)
file.
