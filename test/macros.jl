# macros.jl
# Tests for macros

# Check for changes in Julia's expression parsing
sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j = 1:M; i != j})
@test string(sumexpr) == ":(sum{\$(Expr(:parameters, :((i!=j)))),*(x[i,j],y[i,j]),i = 1:N,j = 1:M})"
@test sumexpr.head == :curly
@test length(sumexpr.args) == 5
@test sumexpr.args[1] == :sum
@test sumexpr.args[2].head == :parameters
@test sumexpr.args[3] == :(x[i,j] * y[i,j])
@test sumexpr.args[4].head == :(=)
@test sumexpr.args[5].head == :(=)


sumexpr = :(sum{x[i,j] * y[i,j], i = 1:N, j = 1:M})
@test string(sumexpr) == ":(sum{*(x[i,j],y[i,j]),i = 1:N,j = 1:M})"
@test sumexpr.head == :curly
@test length(sumexpr.args) == 4
@test sumexpr.args[1] == :sum
@test sumexpr.args[2] == :(x[i,j] * y[i,j])
@test sumexpr.args[3].head == :(=)
@test sumexpr.args[4].head == :(=)

# test JuMP's macros

let 
    m = Model(:Max)
    @defVar(m, w)
    @defVar(m, x)
    @defVar(m, y)
    @defVar(m, z)
    t = 10

    @addConstraint(m, 3x - y == 3.3(w + 2z) + 5) 
    @test conToStr(m.linconstr[end]) == "3.0 x - 1.0 y - 3.3 w - 6.6 z == 5.0"
    @addConstraint(m, (x+y)/2 == 1) 
    @test conToStr(m.linconstr[end]) == "0.5 x + 0.5 y == 1.0"
    @addConstraint(m, -1 <= x-y <= t) 
    @test conToStr(m.linconstr[end]) == "-1.0 <= 1.0 x - 1.0 y <= 10.0"
    @test_throws @addConstraint(m, x <= t <= y)
end

let
    m = Model(:Max)
    @defVar(m, x[1:3,1:3])
    @defVar(m, y)
    C = [1 2 3; 4 5 6; 7 8 9]
    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:2, j = 2:3 } <= 1)
    @test conToStr(m.linconstr[end]) == "2.0 _col2 + 3.0 _col3 + 5.0 _col5 + 6.0 _col6 <= 1.0"
    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:3; i != j} == y)
    @test conToStr(m.linconstr[end]) == "2.0 _col2 + 3.0 _col3 + 4.0 _col4 + 6.0 _col6 + 7.0 _col7 + 8.0 _col8 - 1.0 y == 0.0"

    @addConstraint(m, sum{ C[i,j]*x[i,j], i = 1:3, j = 1:i} == 0);
    @test conToStr(m.linconstr[end]) == "1.0 _col1 + 4.0 _col4 + 5.0 _col5 + 7.0 _col7 + 8.0 _col8 + 9.0 _col9 == 0.0"
end
