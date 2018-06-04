Variables
=========

What is a JuMP VariableRef?
---------------------------

DRAFT: A JuMP variable is a reference to an index in a model. It's a thin
wrapper around MOI.VariableIndex. Variables have (1) names, and (2) attributes.
Describe the different scopes of a variable (e.g., as Julia variables, lookup by
string name, and lookup by symbol).

The `@variable` macro
---------------------

DRAFT: Describe the complete syntax of the `@variable` macro. Anonymous versus
named variables. Describe the three possible container types returned and how
to use them (`Array`, `JuMPArray`, and `Dict`).

How to delete variables.
