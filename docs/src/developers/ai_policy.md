# AI Policy

This document describes how the JuMP developers will manage and respond to
AI-assisted contributions to the JuMP ecosystem.

“AI” herein refers to generative AI tools like large language models (LLMs) that
can generate, edit, and review software code, create and manipulate images, or
generate human-like communication.

The intent of the policy is to balance the benefits of AI-assisted contributions
against the long-term maintenance requirement of the JuMP ecosystem.

It was inspired by similar policies in [SymPy](https://docs.sympy.org/dev/contributing/ai-generated-code-policy.html)
and [SciPy](https://scipy.github.io/devdocs/dev/conduct/ai_policy.html).

## Communication

This part of the policy applies to all communication in our [community forum](https://jump.dev/forum),
[developer chatroom](https://jump.dev/chatroom), and [GitHub repositories](https://github.com/jump-dev).

Do not use AI to generate written communication. Write in your own words.
Human-to-human communication is essential for an open source community to thrive.

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
