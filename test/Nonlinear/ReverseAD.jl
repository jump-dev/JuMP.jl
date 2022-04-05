module TestReverseAD

using Test
import JuMP: MOI
import JuMP: Nonlinear
import JuMP.Nonlinear.ReverseAD: Coloring

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_objective_quadratic_univariate()
    x = MOI.VariableIndex(1)
    data = Nonlinear.NonlinearData()
    Nonlinear.set_objective(data, :($x^2 + 1))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [x],
    )
    MOI.initialize(data, [:Grad, :Jac, :Hess])
    @test MOI.eval_objective(data, [1.2]) == 1.2^2 + 1
    g = [NaN]
    MOI.eval_objective_gradient(data, g, [1.2])
    @test g == [2.4]
    @test MOI.hessian_lagrangian_structure(data) == [(1, 1)]
    H = [NaN]
    MOI.eval_hessian_lagrangian(data, H, [1.2], 1.5, Float64[])
    @test H == 1.5 .* [2.0]
    MOI.eval_hessian_lagrangian_product(data, H, [1.2], [1.2], 1.5, Float64[])
    @test H == [1.5 * 2.0 * 1.2]
    return
end

function test_objective_quadratic_multivariate()
    x = MOI.VariableIndex(1)
    y = MOI.VariableIndex(2)
    data = Nonlinear.NonlinearData()
    Nonlinear.set_objective(data, :($x^2 + $x * $y + $y^2))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [x, y],
    )
    MOI.initialize(data, [:Grad, :Jac, :Hess])
    @test MOI.eval_objective(data, [1.2, 2.3]) == 1.2^2 + 1.2 * 2.3 + 2.3^2
    g = [NaN, NaN]
    MOI.eval_objective_gradient(data, g, [1.2, 2.3])
    @test g == [2 * 1.2 + 2.3, 1.2 + 2 * 2.3]
    @test MOI.hessian_lagrangian_structure(data) == [(1, 1), (2, 2), (2, 1)]
    H = [NaN, NaN, NaN]
    MOI.eval_hessian_lagrangian(data, H, [1.2, 2.3], 1.5, Float64[])
    @test H == 1.5 .* [2.0, 2.0, 1.0]
    v = [0.3, 0.4]
    hv = [NaN, NaN]
    MOI.eval_hessian_lagrangian_product(data, hv, [1.2, 2.3], v, 1.5, Float64[])
    @test hv ≈ 1.5 .* [2 1; 1 2] * v
    return
end

function test_objective_ifelse_comparison()
    x = MOI.VariableIndex(1)
    y = MOI.VariableIndex(2)
    data = Nonlinear.NonlinearData()
    Nonlinear.set_objective(data, :(ifelse(1 <= $x <= 2, $x^2, $y^2)))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [x, y],
    )
    MOI.initialize(data, [:Grad, :Jac, :Hess])
    @test MOI.eval_objective(data, [1.2, 2.3]) == 1.2^2
    @test MOI.eval_objective(data, [2.2, 2.3]) == 2.3^2
    g = [NaN, NaN]
    MOI.eval_objective_gradient(data, g, [1.2, 2.3])
    @test g == [2 * 1.2, 0.0]
    MOI.eval_objective_gradient(data, g, [2.2, 2.3])
    @test g == [0.0, 2 * 2.3]
    return
end

function test_objective_ifelse_logic()
    x = MOI.VariableIndex(1)
    y = MOI.VariableIndex(2)
    data = Nonlinear.NonlinearData()
    Nonlinear.set_objective(data, :(ifelse(1 <= $x && $x <= 2, $x^2, $y^2)))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [x, y],
    )
    MOI.initialize(data, [:Grad, :Jac, :Hess])
    @test MOI.eval_objective(data, [1.2, 2.3]) == 1.2^2
    @test MOI.eval_objective(data, [2.2, 2.3]) == 2.3^2
    g = [NaN, NaN]
    MOI.eval_objective_gradient(data, g, [1.2, 2.3])
    @test g == [2 * 1.2, 0.0]
    MOI.eval_objective_gradient(data, g, [2.2, 2.3])
    @test g == [0.0, 2 * 2.3]
    return
end

function test_objective_parameter()
    x = MOI.VariableIndex(1)
    data = Nonlinear.NonlinearData()
    p = Nonlinear.add_parameter(data, 1.2)
    Nonlinear.set_objective(data, :($p * $x))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [x],
    )
    MOI.initialize(data, [:Grad, :Jac, :Hess])
    @test MOI.eval_objective(data, [1.3]) == 1.2 * 1.3
    g = [NaN]
    MOI.eval_objective_gradient(data, g, [1.3])
    @test g == [1.2]
    return
end

function test_objective_subexpression()
    data = Nonlinear.NonlinearData()
    x = MOI.VariableIndex(1)
    input = :($x^2 + 1)
    expr = Nonlinear.add_expression(data, input)
    expr_2 = Nonlinear.add_expression(data, :($expr^2))
    Nonlinear.set_objective(data, :($expr_2^2))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [x],
    )
    MOI.initialize(data, [:Grad])
    @test MOI.eval_objective(data, [1.3]) == ((1.3^2 + 1)^2)^2
    g = [NaN]
    MOI.eval_objective_gradient(data, g, [1.3])
    @test g ≈ [2 * (1.3^2 + 1)^2 * (2 * (1.3^2 + 1)) * 2 * 1.3]
    return
end

function test_constraint_quadratic_univariate()
    x = MOI.VariableIndex(1)
    data = Nonlinear.NonlinearData()
    Nonlinear.add_constraint(data, :($x^2 <= 2.0))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [x],
    )
    MOI.initialize(data, [:Grad, :Jac, :Hess])
    g = [NaN]
    x_val = [1.2]
    MOI.eval_constraint(data, g, x_val)
    @test g == x_val .^ 2 .- 2
    @test MOI.jacobian_structure(data) == [(1, 1)]
    J = [NaN]
    MOI.eval_constraint_jacobian(data, J, x_val)
    @test J == 2 .* x_val
    @test MOI.hessian_lagrangian_structure(data) == [(1, 1)]
    H = [NaN]
    MOI.eval_hessian_lagrangian(data, H, x_val, 0.0, [1.5])
    @test H == 1.5 .* [2.0]
    return
end

function test_constraint_quadratic_multivariate()
    x = MOI.VariableIndex(1)
    y = MOI.VariableIndex(2)
    data = Nonlinear.NonlinearData()
    Nonlinear.add_constraint(data, :($x^2 + $x * $y + $y^2 <= 2.0))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [x, y],
    )
    MOI.initialize(data, [:Grad, :Jac, :Hess])
    g = [NaN]
    x_val = [1.2, 2.3]
    MOI.eval_constraint(data, g, x_val)
    @test g == [x_val[1]^2 + x_val[1] * x_val[2] + x_val[2]^2] .- 2
    @test MOI.jacobian_structure(data) == [(1, 1), (1, 2)]
    J = [NaN, NaN]
    MOI.eval_constraint_jacobian(data, J, x_val)
    @test J == [2 * x_val[1] + x_val[2], x_val[1] + 2 * x_val[2]]
    @test MOI.hessian_lagrangian_structure(data) == [(1, 1), (2, 2), (2, 1)]
    H = [NaN, NaN, NaN]
    MOI.eval_hessian_lagrangian(data, H, x_val, 0.0, [1.5])
    @test H == 1.5 .* [2.0, 2.0, 1.0]
    return
end

function test_hessian_sparsity_registered_function()
    x = MOI.VariableIndex(1)
    y = MOI.VariableIndex(2)
    z = MOI.VariableIndex(3)
    f(x, y) = x^2 + y^2
    data = Nonlinear.NonlinearData()
    Nonlinear.register_operator(data, :f, 2, f)
    Nonlinear.set_objective(data, :(f($x, $z) + $y^2))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [x, y, z],
    )
    # TODO(odow): re-enable when user-defined hessians are supported
    @test_broken :Hess in MOI.features_available(data)
    # MOI.initialize(data, [:Grad, :Jac, :Hess])
    # @test MOI.hessian_lagrangian_structure(data) ==
    #       [(1, 1), (2, 2), (3, 3), (3, 1)]
    # H = fill(NaN, 4)
    # MOI.eval_hessian_lagrangian(data, H, rand(3), 1.5, Float64[])
    # @test H == 1.5 .* [2.0, 2.0, 2.0, 0.0]
    return
end

struct _ColoringGraph
    num_vertices::Int
    edges::Vector{Tuple{Int,Int}}
end

function to_adjlist(graph::_ColoringGraph)
    I = [i for (i, _) in graph.edges]
    J = [j for (_, j) in graph.edges]
    return Coloring.gen_adjlist(I, J, graph.num_vertices)
end

function test_coloring_edge_free_graph()
    graph = _ColoringGraph(10, [])
    _, numcolors = Coloring.acyclic_coloring(to_adjlist(graph))
    @test numcolors == 1
    return
end

function test_coloring_one_edge_graph()
    graph = _ColoringGraph(10, [(2, 4)])
    color, numcolors = Coloring.acyclic_coloring(to_adjlist(graph))
    @test numcolors == 2
    @test color[2] != color[4]
    return
end

function test_coloring_two_edge_graph()
    graph = _ColoringGraph(10, [(2, 4), (2, 3)])
    color, numcolors = Coloring.acyclic_coloring(to_adjlist(graph))
    @test numcolors == 2
    @test color[3] == color[4]
    return
end

function test_coloring_three_edge_graph()
    graph = _ColoringGraph(10, [(2, 4), (2, 3), (3, 4)])
    color, numcolors = Coloring.acyclic_coloring(to_adjlist(graph))
    @test numcolors == 3
    # TODO: What is this testing?
    Coloring.recovery_preprocess(to_adjlist(graph), color, numcolors, Int[])
    return
end

function test_coloring_two_edge_three_vertex_graph()
    graph = _ColoringGraph(3, [(1, 3), (2, 3)])
    _, numcolors = Coloring.acyclic_coloring(to_adjlist(graph))
    @test numcolors == 2
    return
end

function test_coloring_four_edge_four_vertex_graph()
    graph = _ColoringGraph(4, [(1, 2), (2, 3), (3, 4), (4, 1)])
    _, numcolors = Coloring.acyclic_coloring(to_adjlist(graph))
    @test numcolors == 3
    return
end

function test_coloring_topological_sort()
    # graph = _ColoringGraph(6, [(1, 2), (1, 3), (1, 6), (2, 4), (2, 5)])
    vec = [3, 6, 2, 1, 4, 5, 1, 2, 2, 1]
    offset = [1, 4, 7, 8, 9, 10, 11]
    v = Coloring.reverse_topological_sort_by_dfs(vec, offset, 6, zeros(Int, 6))
    @test v[1] == [3, 6, 4, 5, 2, 1]
    @test v[2] == [0, 1, 1, 2, 2, 1]
    return
end

function test_coloring_end_to_end_hessian_coloring_and_recovery()
    I, J, rinfo = Coloring.hessian_color_preprocess(Set([(1, 2)]), 2)
    R = Coloring.seed_matrix(rinfo)
    Coloring.prepare_seed_matrix!(R, rinfo)
    @test I == [1, 2, 2]
    @test J == [1, 2, 1]
    @test R == [1.0 0.0; 0.0 1.0]
    hess = [3.4 2.1; 2.1 1.3]
    matmat = hess * R
    V = zeros(3)
    Coloring.recover_from_matmat!(V, matmat, rinfo, zeros(3))
    @test V == [3.4, 1.3, 2.1]
    return
end

function test_derivatives()
    a = MOI.VariableIndex(1)
    b = MOI.VariableIndex(2)
    data = Nonlinear.NonlinearData()
    Nonlinear.set_objective(data, :(sin($a^2) + cos($b * 4) / 5 - 2.0))
    Nonlinear.set_differentiation_backend(
        data,
        Nonlinear.SparseReverseMode(),
        [a, b],
    )
    MOI.initialize(data, [:Grad, :Jac, :Hess])
    x = [2.0, 3.0]
    @test MOI.eval_objective(data, x) == sin(x[1]^2) + cos(x[2] * 4) / 5 - 2.0
    g = [NaN, NaN]
    MOI.eval_objective_gradient(data, g, x)
    @test g ≈ [2 * x[1] * cos(x[1]^2), -4 * sin(x[2] * 4) / 5]
    @test MOI.hessian_lagrangian_structure(data) == [(1, 1), (2, 2)]
    H = [NaN, NaN]
    MOI.eval_hessian_lagrangian(data, H, x, 1.5, Float64[])
    H_exact = [
        -4 * x[1]^2 * sin(x[1]^2) + 2 * cos(x[1]^2),
        -4 / 5 * 4 * cos(x[2] * 4),
    ]
    @test H == 1.5 .* H_exact
    return
end

function test_NLPBlockData()
    data = Nonlinear.NonlinearData()
    x = MOI.VariableIndex(1)
    Nonlinear.add_constraint(data, :($x <= 1))
    Nonlinear.add_constraint(data, :($x >= 2))
    Nonlinear.add_constraint(data, :($x == 3))
    Nonlinear.add_constraint(data, :(4 <= $x <= 5))
    block = MOI.NLPBlockData(data, [x])
    @test block.constraint_bounds == [
        MOI.NLPBoundsPair(-Inf, 0.0),
        MOI.NLPBoundsPair(0.0, Inf),
        MOI.NLPBoundsPair(0.0, 0.0),
        MOI.NLPBoundsPair(4.0, 5.0),
    ]
    @test block.has_objective == false
    return
end

end  # module

TestReverseAD.runtests()
