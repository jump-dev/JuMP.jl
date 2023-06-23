![JuMP logo](https://jump.dev/JuMP.jl/dev/assets/logo-with-text-background.svg "JuMP logo")
---

[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://numfocus.org)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jump.dev/JuMP.jl/stable/)
[![In Development](https://img.shields.io/badge/docs-dev-blue.svg)](https://jump.dev/JuMP.jl/dev/)

JuMP is a domain-specific modeling language for [mathematical optimization](https://en.wikipedia.org/wiki/Mathematical_optimization)
embedded in [Julia](https://julialang.org/). You can find out more about us by
visiting [jump.dev](https://jump.dev).


**Latest Release**: [![version](https://juliahub.com/docs/JuMP/DmXqY/1.12.0/version.svg)](https://juliahub.com/ui/Packages/JuMP/DmXqY/1.12.0) (`release-1.0` branch):
  * Installation via the Julia package manager:
    * `import Pkg; Pkg.add("JuMP")`
  * Get help:
    * Read the [Documentation](https://jump.dev/JuMP.jl/stable/): [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jump.dev/JuMP.jl/stable/)
    * Ask a question on the [Community forum](https://jump.dev/forum)
  * Testing status:
    * Github Actions: [![Build Status](https://github.com/jump-dev/JuMP.jl/workflows/CI/badge.svg?branch=release-1.0)](https://github.com/jump-dev/JuMP.jl/actions?query=workflow%3ACI)
  * [![deps](https://juliahub.com/docs/JuMP/deps.svg)](https://juliahub.com/ui/Packages/JuMP/DmXqY?t=2)

**Development version** (`master` branch):
  * Installation via the Julia package manager:
    * `import Pkg; Pkg.add(Pkg.PackageSpec(name="JuMP", rev="master"))`
  * Get help:
    * Read the [Documentation](https://jump.dev/JuMP.jl/dev/): [![In Development](https://img.shields.io/badge/docs-dev-blue.svg)](https://jump.dev/JuMP.jl/dev/)
    * Join the [Developer chatroom](https://jump.dev/chatroom)
    * Read the [NEWS](https://github.com/jump-dev/JuMP.jl/tree/master/NEWS.md)
  * Testing status:
    * Github Actions: [![Build Status](https://github.com/jump-dev/JuMP.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/JuMP.jl/actions?query=workflow%3ACI)
    * Test coverage: [![codecov](https://codecov.io/gh/jump-dev/JuMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/JuMP.jl)

## Need help?

Use the [Community forum](https://jump.dev/forum) to search for answers to
previously asked questions, or ask a new question.

The post [Please read: make it easier to help you](https://discourse.julialang.org/t/please-read-make-it-easier-to-help-you/14757),
describes the best practices for asking a question.

## Bug reports

Please report any issues via the Github [issue tracker]. All types of issues
are welcome and encouraged; this includes bug reports, documentation typos,
feature requests, etc.

[issue tracker]: https://github.com/jump-dev/JuMP.jl/issues

## Citing JuMP

If you find JuMP useful in your work, we kindly request that you cite the
following paper ([preprint](https://arxiv.org/abs/2206.03866)):

```bibtex
@article{Lubin2023,
    author = {Miles Lubin and Oscar Dowson and Joaquim {Dias Garcia} and Joey Huchette and Beno{\^i}t Legat and Juan Pablo Vielma},
    title = {{JuMP} 1.0: {R}ecent improvements to a modeling language for mathematical optimization},
    journal = {Mathematical Programming Computation},
    year = {2023},
    doi = {10.1007/s12532-023-00239-3}
}
```

For earlier works, see:

 * Our paper in SIAM Review ([journal](https://dx.doi.org/10.1137/15M1020575), [pdf](https://mlubin.github.io/pdf/jump-sirev.pdf)):
   ```bibtex
   @article{DunningHuchetteLubin2017,
       author = {Iain Dunning and Joey Huchette and Miles Lubin},
       title = {{JuMP}: {A} {M}odeling {L}anguage for {M}athematical {O}ptimization},
       journal = {SIAM Review},
       volume = {59},
       number = {2},
       pages = {295-320},
       year = {2017},
       doi = {10.1137/15M1020575},
   }
   ```

 * Our paper in IJOC ([journal](https://dx.doi.org/10.1287/ijoc.2014.0623), [preprint](https://arxiv.org/abs/1312.1431)):
   ```bibtex
   @article{LubinDunningIJOC,
       author = {Miles Lubin and Iain Dunning},
       title = {{C}omputing in {O}perations {R}esearch {U}sing {J}ulia},
       journal = {INFORMS Journal on Computing},
       volume = {27},
       number = {2},
       pages = {238-248},
       year = {2015},
       doi = {10.1287/ijoc.2014.0623},
   }
   ```

---

![NumFOCUS logo](https://jump.dev/JuMP.jl/dev/assets/numfocus-logo.png)

JuMP is a Sponsored Project of NumFOCUS, a 501(c)(3) nonprofit charity in the
United States. NumFOCUS provides JuMP with fiscal, legal, and administrative
support to help ensure the health and sustainability of the project. Visit
[numfocus.org](https://numfocus.org) for more information.

You can support JuMP by [donating](https://numfocus.salsalabs.org/donate-to-jump/index.html).

Donations to JuMP are managed by NumFOCUS. For donors in the United States,
your gift is tax-deductible to the extent provided by law. As with any donation,
you should consult with your tax adviser about your particular tax situation.

JuMP's largest expense is the annual JuMP-dev workshop. Donations will help us
provide travel support for JuMP-dev attendees and take advantage of other
opportunities that arise to support JuMP development.
