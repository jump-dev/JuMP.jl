# How to contribute to JuMP

Welcome! This document explains some of the ways you can contribute to JuMP.

## Code of Conduct

This project and everyone participating in it is governed by the
[JuMP Code of Conduct](https://github.com/JuliaOpt/JuMP.jl/blob/master/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code.

## Join the Discourse forum

First up, join the [Discourse forum](https://discourse.julialang.org/c/domain/opt).

The forum is a good place to ask questions about how to use JuMP. You can also
use the forum to discuss possible feature requests and bugs before raising a
Github issue (more on this below).

Aside from asking questions, the easiest way you can contribute to JuMP is to
help answer questions on Discourse!

## Improve the documentation

Chances are, if you asked (or answered) a question on Discourse, then it is
usually a sign that the [documentation](http://www.juliaopt.org/JuMP.jl/latest/)
could be improved. Moreover, since it is your question, you are probably the
best-placed person to improve it!

The docs are written in Markdown and are built using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).
You can find the source of all the docs [here](https://github.com/JuliaOpt/JuMP.jl/tree/master/docs).

If your change is small (like fixing typos, or one or two sentence corrections),
the easiest way to do this is via Github's online editor. (Github has
[help](https://help.github.com/articles/editing-files-in-another-user-s-repository/)
on how to do this.)

If your change is larger, or touches multiple files, you will need to make the
change locally and then use Git to submit a pull request. (See [Contribute code to JuMP](#Contribute code to JuMP)
below for more on this.)

If you need any help, come join our [Gitter](https://gitter.im/JuliaOpt/JuMP-dev)
channel and we will walk you through the process.

## File a bug report

A third way to contribute to JuMP is to file [bug reports](https://github.com/JuliaOpt/JuMP.jl/issues/new?template=bug_report.md).

(Make sure you read the info in the box where you write the body of the issue
before posting. You can also find a copy of that info [here](https://github.com/JuliaOpt/JuMP.jl/blob/master/.github/ISSUE_TEMPLATE/bug_report.md).)

## Contribute code to JuMP

Finally, you can also contribute code to JuMP!

If you don't have experience with Git, Github, and Julia development, the first
steps can be a little daunting. However, there are lots of tutorials available
online (such as [this](http://try.github.io/), [this](https://guides.github.com/activities/hello-world/),
and [this](https://docs.julialang.org/en/stable/manual/packages/#Making-changes-to-an-existing-package-1)).
If you need any help, come join our [Gitter](https://gitter.im/JuliaOpt/JuMP-dev)
channel and we will walk you through the process.

Once you are familiar with Git, the workflow for contributing code to JuMP is
along the lines of the following:
1. Find an [open issue](https://github.com/JuliaOpt/JuMP.jl/issues) (or open a
    new one) for the problem you want to solve
2. Discuss (in the issue, or on [Gitter](https://gitter.im/JuliaOpt/JuMP-dev))
_before_ spending too much time on it to test the waters first and see if other
contributors are fine with these changes or not
3. Make your changes locally. (The Julia manual has a [guide](https://docs.julialang.org/en/stable/manual/packages/#Making-changes-to-an-existing-package-1)
on how to do this.) Make sure you:
    - Follow the [style guide](http://www.juliaopt.org/JuMP.jl/latest/style.html)
    - Add tests and documentation for any changes or new features
4. Submit a pull request to Github
5. Update your pull request, responding to any comments

Remember to be patient and polite; you may get a _lot_ of comments on your pull
request! However, don't be afraid! A lot of comments means that people are
willing to help you improve the code that you are contributing to JuMP.

Thanks for contributing to JuMP!
