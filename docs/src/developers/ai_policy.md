# AI Policy

This document describes how the JuMP developers will manage and respond to
AI-assisted contributions to the JuMP ecosystem.

It was inspired by similar policies in [SymPy](https://docs.sympy.org/dev/contributing/ai-generated-code-policy.html)
and [SciPy](https://scipy.github.io/devdocs/dev/conduct/ai_policy.html).

## Communication

This part of the policy applies to all communication in our[community forum](https://jump.dev/forum),
[developer chatroom](https://jump.dev/chatroom), and [GitHub repositories](https://github.com/jump-dev).

Do not use AI to speak for you by copy-pasting a conversation in-and-out of a
chatbot.

We would rather speak to you directly, errors and grammatical mistakes included.

If English is not your first language, consider providing both a machine
generated translation into English and the original text in your preferred
language.

## Code contributions

This part of the policy applies to all pull requests and issues in our
[GitHub repositories](https://github.com/jump-dev).

We expect and encourage users to use AI assistance when developing. You do not
need to disclose what parts you used AI assistance for, but you are responsible
for all code that you submit to JuMP. However, do not use AI to author the text
descriptions of pull requests and issues; see our [Communication](@ref) policy.

When authoring new code in JuMP, keep in mind that the JuMP developers' two
biggest bottlenecks are:

1. capacity for code review of new pull requests
2. on-going support and maintenance of existing features.

For these reasons, please open an issue to discuss what you want to change
_before_ opening a pull request. In the issue, describe what you want to change
and why it matters to you.

Large pull requests that add new features not previously discussed in an issue
may be closed without review, even if they are correct and demonstrably useful.

## New Contributors

The JuMP developers welcome and encourage new contributors.

The best ways you can get involved are:

1. **Find and report bugs**: we can't fix things that we don't know about. There
   are always new bugs (or inconsistencies that we should better document) to
   find. Report a bug by [opening a GitHub issue](https://github.com/jump-dev/JuMP.jl/issues).

   Using AI tools to find bugs is acceptable, but you are responsible for
   understanding and explaining _why_ the issue is a bug (see our
   [Communication](@ref) policy). Rather than copy-pasting only the content of
   the AI analysis into the issue, first explain what you did to find and
   verify the bug, and then copy-paste the AI tool's output.

   Also note that there are many repositories in the JuMP ecosystem. Don't worry
   if you open an issue in the "wrong" one; we can easily transfer it to the
   correct repository.

2. **Tell us about confusing parts of the documentation**: if you get stuck
   trying to do something in JuMP, it means we didn't document things well
   enough. If you have suggestions for new tutorials we could add or how we
   could improve the documentation please leave a comment at
   ["Suggestions for documentation improvements"](https://github.com/jump-dev/JuMP.jl/issues/2348).

3. **Make your own packages**: write your own solver or JuMP extension under
   your personal account. There are no rules here. Vibe code as much as you
   like. Tell us about the things you have created by posting on the
   [community forum](https://jump.dev/forum), or give a talk at a
   [JuMP-dev workshop](https://jump.dev/categories/#jump-dev).
