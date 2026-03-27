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

## Responsibility

You are responsible for any code you submit to JuMP's repositories. You must
understand and be able to explain the code you submit as well as the existing
related code. It is not acceptable to submit a patch that you cannot understand
and explain in your own words. In explaining your contribution, do not use AI to
automatically generate descriptions; see our [Communication](@ref) policy.

## Copyright

All code in JuMP-related repositories is released under an open source license
(the exact license depends on the repository). Contributors license their code
under the same license. That means contributors must own the copyright of any
code you submit. It is your responsibility to not infringe on others copyright.
We will reject any pull requests where the copyright is in question.

## Disclosure

You must disclose whether AI has been used to assist in the development of your
pull request. If so, you must document which tools have been used, how they
were used, and specify what code or text is AI generated. We will reject any
pull request that does not include the disclosure.

## Code Quality

When authoring new code in JuMP, keep in mind that the JuMP developers' two
biggest bottlenecks are:

1. capacity for code review of new pull requests
2. on-going support and maintenance of existing features.

For these reasons, please open an issue to discuss what you want to change
_before_ opening a pull request. In the issue, describe what you want to change
and why it matters to you.

Large pull requests that add new features not previously discussed in an issue
may be closed without review, even if they are correct and demonstrably useful.

## AI Agents

The use of an AI agent that writes code and then submits a pull request
autonomously is not permitted. A human must check any generated code and submit
a pull request according to the [Responsibility](@ref) section above.
