# Checklists

The purpose of this page is to collate a series of checklists for commonly
performed changes to the source code of JuMP.

In each case, copy the checklist into the description of the pull request.

## Making a release

In preparation for a release, use the following checklist. These steps can be
done in the same commit, or separately. The last commit should have the message
"Prep for vX.Y.Z."

````
## Pre-release

 - [ ] Check that the pinned packages in `docs/Project.toml` are updated. We pin
       the versions so that changes in the solvers (changes in printing, small
       numeric changes) do not break the printing of the JuMP docs in arbitrary
       commits.
 - [ ] Check that the `rev` fields in `docs/packages.toml` are updated. We pin
       the versions of solvers and extensions to ensure that changes to their
       READMEs do not break the JuMP docs in arbitrary commits, and to ensure
       that the versions are compatible with the latest JuMP and
       MathOptInterface releases.
 - [ ] Check compat of `DimensionalData` in `Project.toml`
 - [ ] Check compat of `MacroTools` in `Project.toml`
 - [ ] Update `docs/src/changelog.md`
 - [ ] Run https://github.com/jump-dev/JuMP.jl/actions/workflows/extension-tests.yml
       using a `workflow_dispatch` trigger to check for any changes in JuMP that
       broke extensions.
 - [ ] Change the version number in `Project.toml`
 - [ ] The commit messages in this PR do not contain `[ci skip]`

## The release

 - [ ] After merging this pull request, comment `[at]JuliaRegistrator register` in
       the GitHub commit. This should automatically publish a new version to the
       Julia registry, as well as create a tag, and rebuild the documentation
       for this tag.

       These steps can take quite a bit of time (1 hour or more), so don't be
       surprised if the new documentation takes a while to appear. In addition,
       the links in the README will be broken until JuliaHub fetches the new
       version on their servers.

## Post-release

 - [ ] Once the tag is created, update the relevant `release-` branch. The latest
       release branch at the time of writing is `release-1.0` (we haven't
       back-ported any patches that needed to create a `release-1.Y` branch). To
       to update the release branch with the v1.10.0 tag, do:
       ```
       git checkout release-1.0
       git pull
       git merge v1.10.0
       git push
       ```
````

## Adding a new solver to the documentation

Use the following checklist when adding a new solver to the JuMP documentation.

````
## Basic

 - [ ] Check that the solver is a registered Julia package
 - [ ] Check that the solver supports the long-term support release of Julia
 - [ ] Check that the solver has a MathOptInterface wrapper
 - [ ] Check that the tests call `MOI.Test.runtests`. Some test excludes are
       permissible, but the reason for skipping a particular test should be
       documented.
 - [ ] Check that the README and/or documentation provides an example of how to
       use the solver with JuMP

## Documentation

 - [ ] Add a new row to the table in `docs/src/installation.md`

## Optional

 - [ ] Add package metadata to `docs/packages.toml`
````

## Adding a new shape

Use the following checklist when adding a new `AbstractShape`

````
## Basic

 - [ ] Add a new subtype of `AbstractShape`
 - [ ] Implement `vectorize(data, ::NewShape)::Vector`
 - [ ] Implement `reshape_vector(vector, ::NewShape)`
 - [ ] Implement `dual_shape`, or verify that the shape is self-dual
 - [ ] Add the tests from https://github.com/jump-dev/JuMP.jl/pull/3816
````
