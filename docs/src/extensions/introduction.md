# Introduction

This section of the documentation contains brief documentation for some popular
JuMP extensions. The list of extensions is not exhaustive, but instead is
intended to help you discover popular JuMP extensions, and to give you an
overview of the types of extensions that are possible to write with JuMP.

## Affiliation

Packages beginning with `jump-dev/` are developed and maintained by the
JuMP developers.

Packages that do not begin with `jump-dev/` are developed independently. The
developers of these packages requested or consented to the inclusion of their
README contents in the JuMP documentation for the benefit of users.

## Adding new extensions

Written an extension? Add it to this section of the JuMP documentation by making
a pull request to the [`docs/packages.toml`](https://github.com/jump-dev/JuMP.jl/blob/master/docs/packages.toml)
file.

## Weak dependencies

Some extensions listed in this section are implemented using the [weak dependency](https://pkgdocs.julialang.org/v1/creating-packages/#Weak-dependencies)
feature added to Julia in v1.9. These extensions are activated if and only if
you have `JuMP` and the other package loaded into your current scope with
`using` or `import`.

!!! compat
    Using a weak dependency requires Julia v1.9 or later.
