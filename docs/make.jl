#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

import Distributed

Distributed.@everywhere include(joinpath(@__DIR__, "make_utilities.jl"))

# Needed to make Documenter think that there is a PDF in the right place when
# link checking. In production we replace this by running the LaTeX build.
write(joinpath(@__DIR__, "src", "JuMP.pdf"), "")

@time literate_tutorials()

@time begin
    add_tutorial_overview()
    add_solver_readmes()
    add_changelog()
    add_moi_pages()
end

@time if _PDF || _IS_GITHUB_ACTIONS
    # Copy the src into a new directory for the latex build
    cp(joinpath(@__DIR__, "src"), joinpath(@__DIR__, "latex_src"); force = true)
    f_latex = Distributed.@spawnat :any make_latex()
    make_html()
    fetch(f_latex)
    # Hack for deploying: copy the pdf (and only the PDF) into the HTML build
    # directory! We don't want to copy everything in `latex_build` because it
    # includes lots of extraneous LaTeX files.
    cp(
        joinpath(@__DIR__, "latex_build", "JuMP.pdf"),
        joinpath(@__DIR__, "build", "JuMP.pdf");
        force = true,
    )
    rm(joinpath(@__DIR__, "latex_src"); force = true, recursive = true)
else
    @time make_html()
end

if _IS_GITHUB_ACTIONS
    Documenter.deploydocs(;
        repo = "github.com/jump-dev/JuMP.jl.git",
        push_preview = true,
    )
end
