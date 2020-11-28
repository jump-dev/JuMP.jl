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
    return esc(
        quote
            @test JuMP.isequal_canonical(@expression(model, $expr), $expr)
        end,
    )
end

macro test_expression_with_string(expr, str)
    return esc(
        quote
            realized_expr = @inferred $expr
            @test string(realized_expr) == $str
            @test JuMP.isequal_canonical(
                @expression(model, $expr),
                realized_expr,
            )
        end,
    )
end

function _strip_line_from_error(err::ErrorException)
    return ErrorException(replace(err.msg, r"^At.+\:[0-9]+\: `@" => "In `@"))
end
_strip_line_from_error(err::LoadError) = _strip_line_from_error(err.error)
_strip_line_from_error(err) = err

# Test that the macro call `m` throws an error exception during pre-compilation
macro test_macro_throws(errortype, m)
    # See https://discourse.julialang.org/t/test-throws-with-macros-after-pr-23533/5878
    quote
        @test_throws(
            $(esc(_strip_line_from_error(errortype))),
            try
                @eval $m
            catch err
                throw(_strip_line_from_error(err))
            end
        )
    end
end

# Test that the macro call `m` throws an error exception during _runtime_.
macro test_throws_strip(errortype, m)
    quote
        @test_throws(
            $(esc(_strip_line_from_error(errortype))),
            try
                $(esc(m))
            catch err
                throw(_strip_line_from_error(err))
            end
        )
    end
end
