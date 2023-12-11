#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# test/utilities.jl
# Utility functions used across JuMP's tests.
#############################################################################

macro test_expression(expr)
    return esc(quote
        @test isequal_canonical(@expression(model, $expr), $expr)
    end)
end

macro test_expression_with_string(expr, str, inferrable = true)
    code = quote
        realized_expr = if $inferrable
            @inferred $expr
        else
            $expr
        end
        @test string(realized_expr) == $str
        @test isequal_canonical(@expression(model, $expr), realized_expr)
    end
    return esc(code)
end

function _strip_line_from_error(err::ErrorException)
    return ErrorException(replace(err.msg, r"^At.+\:[0-9]+\: `@" => "In `@"))
end

_strip_line_from_error(err::LoadError) = _strip_line_from_error(err.error)

_strip_line_from_error(err) = err

"""
    @test_throws_parsetime(error_type, expression)

Test that the macro call `expression` throws an `error_type` exception during
the parsing of the macro.

See https://discourse.julialang.org/t/test-throws-with-macros-after-pr-23533
"""
macro test_throws_parsetime(error_type, expression)
    return quote
        @test_throws(
            $(esc(_strip_line_from_error(error_type))),
            try
                @eval $expression
            catch err
                throw(_strip_line_from_error(err))
            end
        )
    end
end

"""
    @test_throws_runtime(error_type, expression)

Test that the macro call `expression` throws an `error_type` exception during
_runtime_.
"""
macro test_throws_runtime(error_type, expression)
    return quote
        @test_throws(
            $(esc(_strip_line_from_error(error_type))),
            try
                $(esc(expression))
            catch err
                throw(_strip_line_from_error(err))
            end
        )
    end
end
