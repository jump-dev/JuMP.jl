# Release process

This note serves as developer documentation for the JuMP release process, which
has accumulated a number of manual steps.

## Pre-release checks and updates

Before making a release:

 * Check that the pinned packages in `docs/Project.toml` are updated. We pin the
   versions so that changes in the solvers (changes in printing, small numeric
   changes) do not break the printing of the JuMP docs in arbitrary commits.

 * Check that the `rev` fields in `docs/packages.toml` are updated. We pin the
   versions of solvers and extensions to ensure that changes to their READMEs do
   not break the JuMP docs in arbitrary commits, and to ensure that the versions
   are compatible with the latest JuMP and MathOptInterface releases.

 * Update `docs/src/changelog.md`.

 * Run [`.github/workflows/extension-tests.yml`](https://github.com/jump-dev/JuMP.jl/actions/workflows/extension-tests.yml)
   using a `workflow_dispatch` trigger to check for any changes in JuMP that
   broke extensions.

 * Change the version number in `Project.toml`, and update the links in
   README.md.

These steps can be done in the same commit, or separately. The last commit
should have the message "Prep for vX.Y.Z." It must not include the text
`[ci skip]`, or a tag will not be created, see [Documenter#965](https://github.com/JuliaDocs/Documenter.jl/issues/965).

## Making a release

In the GitHub commit, comment `@JuliaRegistrator register`. This should
automatically publish a new version to the Julia registry, as well as create a
tag, and rebuild the documentation for this tag.

These steps can take quite a bit of time (1 hour or more), so don't be surprised
if the new documentation takes a while to appear. In addition, the links in the
README will be broken until JuliaHub fetches the new version on their servers.

## Post-release

Once the tag is created, you can update the relevant `release-` branch. The
latest release branch at the time of writing is `release-1.0` (we haven't
back-ported any patches that needed to create a `release-1.Y` branch), so to
update the release branch with the v1.10.0 tag, do:
```
git checkout release-1.0
git pull
git merge v1.10.0
git push
```
