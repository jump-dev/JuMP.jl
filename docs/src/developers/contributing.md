# How to contribute to JuMP

Welcome, this document explains some ways you can contribute to JuMP.

## Code of Conduct

This project and everyone participating in it is governed by the
[JuMP Code of Conduct](https://github.com/jump-dev/JuMP.jl/blob/master/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code.

## Join the community forum

First up, join the [community forum](https://jump.dev/forum).

The forum is a good place to ask questions about how to use JuMP. You can also
use the forum to discuss possible feature requests and bugs before raising a
GitHub issue (more on this below).

Aside from asking questions, the easiest way you can contribute to JuMP is to
help answer questions on the forum.

## Join the developer chatroom

If you're interested in contributing code to JuMP, the next place to join is the
[developer chatroom](https://jump.dev/chatroom). Let us know what you
have in mind, and we can point you in the right direction.

## Improve the documentation

Chances are, if you asked (or answered) a question on the community forum, then
it is a sign that the [documentation](https://jump.dev/JuMP.jl/dev/) could be
improved. Moreover, since it is your question, you are probably the best-placed
person to improve it.

The docs are written in Markdown and are built using
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).
You can find the source of all the docs
[here](https://github.com/jump-dev/JuMP.jl/tree/master/docs).

If your change is small (like fixing typos, or one or two sentence corrections),
the easiest way to do this is via GitHub's online editor. (GitHub has
[help](https://help.github.com/articles/editing-files-in-another-user-s-repository/)
on how to do this.)

If your change is larger, or touches multiple files, you will need to make the
change locally and then use Git to submit a pull request. (See
[Contribute code to JuMP](@ref) below for more on this.)

!!! tip
    If you need any help, come join the
    [developer chatroom](https://jump.dev/chatroom) and we will walk
    you through the process.

## File a bug report

Another way to contribute to JuMP is to file
[bug reports](https://github.com/jump-dev/JuMP.jl/issues/new?template=bug_report.md).

Make sure you read the info in the box where you write the body of the issue
before posting. You can also find a copy of that info
[here](https://github.com/jump-dev/JuMP.jl/blob/master/.github/ISSUE_TEMPLATE/bug_report.md).

!!! tip
    If you're unsure whether you have a real bug, post on the
    [community forum](https://jump.dev/forum)
    first. Someone will either help you fix the problem, or let you know the
    most appropriate place to open a bug report.

## Contribute code to JuMP

Finally, you can also contribute code to JuMP.

!!! warning
    If you do not have experience with Git, GitHub, and Julia development, the
    first steps can be a little daunting. However, there are lots of tutorials
    available online, including these for:
     * [GitHub](https://guides.github.com/activities/hello-world/)
     * [Git and GitHub](https://try.github.io/)
     * [Git](https://git-scm.com/book/en/v2)
     * [Julia package development](https://docs.julialang.org/en/v1/stdlib/Pkg/#Developing-packages-1)
    If you need any help, come join the [developer chatroom](https://jump.dev/chatroom)
    and we will walk you through the process.

Once you are familiar with Git and GitHub, the workflow for contributing code to
JuMP is similar to the following:

**Step 1: decide what to work on**

The first step is to find an [open issue](https://github.com/jump-dev/JuMP.jl/issues)
(or open a new one) for the problem you want to solve. Then, _before_ spending
too much time on it, discuss what you are planning to do in the issue to see if
other contributors are fine with your proposed changes. Getting feedback early can
improve code quality, and avoid time spent writing code that does not get merged into
JuMP.

!!! tip
    At this point, remember to be patient and polite; you may get a _lot_ of
    comments on your issue. However, do not be afraid. Comments mean that people are
    willing to help you improve the code that you are contributing to JuMP.

**Step 2: fork JuMP**

Go to [https://github.com/jump-dev/JuMP.jl](https://github.com/jump-dev/JuMP.jl)
and click the "Fork" button in the top-right corner. This will create a copy of
JuMP under your GitHub account.

**Step 3: install JuMP locally**

Open Julia and run:
```julia
] dev JuMP
```
This will download the JuMP Git repository to `~/.julia/dev/JuMP`. If you're on
Windows, this will be `C:\\Users\\<my_name>\\.julia\\dev\\JuMP`.

!!! warning
    `] command` means "first type `]` to enter the Julia pkg mode, then type the
    rest. Don't copy-paste the code directly.

**Step 4: checkout a new branch**

!!! note
    In the following, replace any instance of `GITHUB_ACCOUNT` with your GitHub
    user name.

The next step is to checkout a development branch. In a terminal (or command
prompt on Windows), run:
```
$ cd ~/.julia/dev/JuMP

$ git remote add GITHUB_ACCOUNT https://github.com/GITHUB_ACCOUNT/JuMP.jl.git

$ git checkout master

$ git pull

$ git checkout -b my_new_branch
```

!!! tip
    Lines starting with `$` mean "run these in a terminal (command prompt on
    Windows)."

**Step 5: make changes**

Now make any changes to the source code inside the `~/.julia/dev/JuMP`
directory.

Make sure you:
 * Follow the [Style guide](@ref) and run [JuliaFormatter](@ref)
 * Add tests and documentation for any changes or new features

!!! tip
    When you change the source code, you'll need to restart Julia for the
    changes to take effect. This is a pain, so install
    [Revise.jl](https://github.com/timholy/Revise.jl).

**Step 6a: test your code changes**

To test that your changes work, run the JuMP test-suite by opening Julia and
running:
```julia
cd("~/.julia/dev/JuMP")
] activate .
] test
```

!!! warning
    Running the tests might take a long time (~10--15 minutes).

!!! tip
    If you're using Revise.jl, you can also run the tests by calling `include`:
    ```julia
    include("test/runtests.jl")
    ```
    This can be faster if you want to re-run the tests multiple times.

**Step 6b: test your documentation changes**

Open Julia, then run:
```julia
cd("~/.julia/dev/JuMP/docs")
] activate .
include("src/make.jl")
```

!!! warning
    Building the documentation might take a long time (~10 minutes).

!!! tip
    If there's a problem with the tests that you don't know how to fix, don't
    worry. Continue to step 5, and one of the JuMP contributors will comment
    on your pull request telling you how to fix things.

**Step 7: make a pull request**

Once you've made changes, you're ready to push the changes to GitHub. Run:
```
$ cd ~/.julia/dev/JuMP

$ git add .

$ git commit -m "A descriptive message of the changes"

$ git push -u GITHUB_ACCOUNT my_new_branch
```

Then go to [https://github.com/jump-dev/JuMP.jl](https://github.com/jump-dev/JuMP.jl)
and follow the instructions that pop up to open a pull request.

**Step 8: respond to comments**

At this point, remember to be patient and polite; you may get a _lot_ of
comments on your pull request. However, do not be afraid. A lot of comments
means that people are willing to help you improve the code that you are
contributing to JuMP.

To respond to the comments, go back to step 5, make any changes, test the
changes in step 6, and then make a new commit in step 7. Your PR will
automatically update.

**Step 9: cleaning up**

Once the PR is merged, clean-up your Git repository ready for the
next contribution.
```
$ cd ~/.julia/dev/JuMP

$ git checkout master

$ git pull
```

!!! note
    If you have suggestions to improve this guide, please make a pull request.
    It's particularly helpful if you do this after your first pull request
    because you'll know all the parts that could be explained better.
