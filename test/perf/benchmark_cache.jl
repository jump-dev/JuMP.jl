# Benchmark for PR #4032: caching strategies for moi_function on
# GenericNonlinearExpr.
#
# Four modes are compared:
#   :none           — master behaviour: NO cache. `add_constraint` runs
#                      `check_belongs_to_model` and then `moi_function`, both
#                      of which traverse the expression independently.
#   :per_call_oid   — odow's actual PR snippet: a fresh `Dict{UInt64, ...}`
#                      allocated INSIDE every call to `moi_function`, keyed by
#                      `objectid(arg)` (so cache lookups are O(1)).
#                      `check_belongs_to_model` still runs as a separate walk.
#   :model_struct   — this branch as-is: cache stored on the model, keyed by
#                      the JuMP `GenericNonlinearExpr` itself. Default `==`/
#                      `hash` on that struct is structural and walks the whole
#                      sub-tree, so a cache hit still costs O(size) — which
#                      can defeat the cache for deeply aliased DAGs.
#   :model_oid      — model-level cache, but keyed by `objectid` like odow's
#                      snippet. Shows what the branch could evolve to once
#                      blegat's "use hash as keys" suggestion is applied.
#
# Run from the repo root, in an environment that has BenchmarkTools added on
# top of dev'd JuMP (do NOT add BenchmarkTools to JuMP's own Project.toml):
#
#   julia> using Pkg
#   julia> Pkg.activate(temp = true)
#   julia> Pkg.develop(path = ".")
#   julia> Pkg.add("BenchmarkTools")
#   julia> include("benchmark_cache.jl"); run_all()
#
# The script monkey-patches `JuMP.moi_function(::GenericNonlinearExpr,
# ::GenericModel)` so the same JuMP build can be exercised under each mode
# without rebuilding the package.
#
# Sample numbers from this branch (Julia 1.12, K=14, etc.):
#
#   A: aliased tree (one big DAG, K=14)        ~46 ms — all four modes within
#                                               1% (MOI-tree alloc cost is
#                                               small vs. the rest of
#                                               add_constraint here).
#   B: many independent constraints (N=5000)   none 84 / per_call_oid 86 /
#                                               model_struct 91 / model_oid 86
#                                               → the *current branch* is the
#                                               slowest; struct keys hurt the
#                                               common case. Switching to
#                                               objectid keys fixes it.
#   C: shared big subexpr (N=200, M=200)       model_struct 171 / model_oid 170
#                                               vs none 191 / per_call_oid 195
#                                               → 12–14% win for the model-
#                                               level cache. Per-call cannot
#                                               see the cross-constraint
#                                               sharing.
#   D: many aliased trees (M=200, K=8)         all ~140 ms (per-constraint
#                                               sharing is small and JuMP's
#                                               `+` would flatten without the
#                                               explicit GenericNonlinearExpr
#                                               construction we use here).
#
# Takeaways:
#  * The model-level scope is the right call when subexpressions are shared
#    across constraints (scenario C).
#  * Using the JuMP struct as the cache key (current branch) is what costs
#    in scenario B — the deep `==`/`hash` adds overhead per node. Switching
#    to `objectid` keys (blegat's "use hash as keys" suggestion) recovers
#    most of it.
#  * Scenario A's exponential blow-up is more visible in MEMORY than in
#    TIME at K=14; bump K higher (or use the bench.jl example with
#    |R|*|S| ≈ 3600) to make the time gap obvious.

using Printf
using JuMP
using BenchmarkTools
import MathOptInterface as MOI

const _G = JuMP.GenericNonlinearExpr

# ---------------------------------------------------------------------------
# Three implementations of the GenericNonlinearExpr → ScalarNonlinearFunction
# walk. They share the same skeleton; only the cache differs.
# ---------------------------------------------------------------------------

function _moi_no_cache(f::_G{V}) where {V}
    ret = MOI.ScalarNonlinearFunction(f.head, similar(f.args))
    stack = Tuple{MOI.ScalarNonlinearFunction,Int,_G{V}}[]
    for i in length(f.args):-1:1
        if f.args[i] isa _G{V}
            push!(stack, (ret, i, f.args[i]))
        else
            ret.args[i] = JuMP.moi_function(f.args[i])
        end
    end
    while !isempty(stack)
        parent, i, arg = pop!(stack)
        child = MOI.ScalarNonlinearFunction(arg.head, similar(arg.args))
        parent.args[i] = child
        for j in length(arg.args):-1:1
            if arg.args[j] isa _G{V}
                push!(stack, (child, j, arg.args[j]))
            else
                child.args[j] = JuMP.moi_function(arg.args[j])
            end
        end
    end
    return ret
end

# Cache walk parameterised by the `keyfn` that maps a JuMP NL expression to
# the dict key. `keyfn = identity` matches the branch (struct keys, deep
# hash). `keyfn = objectid` matches odow's snippet (UInt64 keys, O(1) hash).
function _moi_with_cache(f::_G{V}, cache, keyfn::F) where {V,F}
    fk = keyfn(f)
    if haskey(cache, fk)
        return cache[fk]
    end
    ret = MOI.ScalarNonlinearFunction(f.head, similar(f.args))
    stack = Tuple{MOI.ScalarNonlinearFunction,Int,_G{V}}[]
    for i in length(f.args):-1:1
        if f.args[i] isa _G{V}
            push!(stack, (ret, i, f.args[i]))
        else
            ret.args[i] = JuMP.moi_function(f.args[i])
        end
    end
    while !isempty(stack)
        parent, i, arg = pop!(stack)
        argk = keyfn(arg)
        if haskey(cache, argk)
            parent.args[i] = cache[argk]
            continue
        end
        child = MOI.ScalarNonlinearFunction(arg.head, similar(arg.args))
        parent.args[i] = child
        for j in length(arg.args):-1:1
            if arg.args[j] isa _G{V}
                push!(stack, (child, j, arg.args[j]))
            else
                child.args[j] = JuMP.moi_function(arg.args[j])
            end
        end
        cache[argk] = child
    end
    cache[fk] = ret
    return ret
end

# ---------------------------------------------------------------------------
# Switching mechanism. We re-define JuMP's `moi_function(::GenericNonlinearExpr,
# ::GenericModel)` so that `add_constraint` dispatches differently per mode.
#
# Modes :none and :per_call mirror what those PR options would actually
# look like in production: a separate `check_belongs_to_model` walk, then
# the (cached or uncached) build walk.
# Mode :model_level matches the branch: a single integrated walk.
# ---------------------------------------------------------------------------

const BENCH_MODE = Ref(:model_struct)

# Per-mode caches that need a UInt64 (objectid) key type rather than the
# branch's `Dict{Any,...}`. We store them on the model as a side table so
# they survive across constraints inside one build but are cleared between
# builds (because each scenario builds a fresh `Model()`).
const _OID_CACHE = IdDict{JuMP.GenericModel,Dict{UInt64,MOI.ScalarNonlinearFunction}}()

function _oid_cache(model)
    get!(_OID_CACHE, model) do
        Dict{UInt64,MOI.ScalarNonlinearFunction}()
    end
end

@eval JuMP function moi_function(
    f::JuMP.GenericNonlinearExpr{V},
    model::JuMP.GenericModel,
) where {V}
    mode = Main.BENCH_MODE[]
    if mode === :none
        JuMP.check_belongs_to_model(f, model)
        return Main._moi_no_cache(f)
    elseif mode === :per_call_oid
        JuMP.check_belongs_to_model(f, model)
        cache = Dict{UInt64,$(MOI.ScalarNonlinearFunction)}()
        return Main._moi_with_cache(f, cache, objectid)
    elseif mode === :model_struct
        # branch behaviour: belongs-check is integrated; struct keys (deep hash)
        return Main._moi_with_cache(f, model.subexpressions, identity)
    else  # :model_oid
        return Main._moi_with_cache(f, Main._oid_cache(model), objectid)
    end
end

# ---------------------------------------------------------------------------
# Scenarios. Each `scenario_*` returns a 0-arg closure that builds a fresh
# model from scratch, so each BenchmarkTools sample sees an empty
# `model.subexpressions` cache.
# ---------------------------------------------------------------------------

# Scenario A — odow's case: ONE constraint with a deeply aliased binary tree.
# `e_k = e_{k-1} + e_{k-1}`. With no cache this expands to 2^K leaves; with
# either cache it stays linear in K.
#
# We build the `+` node directly via the GenericNonlinearExpr constructor so
# JuMP's `+` operator overload does NOT flatten the n-ary sum. Without this,
# `e + e` returns a flat `:+(args...)` of 2^K identical leaf references and
# the master walk only iterates 2^K times instead of doing 2^K *recursive*
# descents through K levels — which would mask the worst case the PR is
# meant to fix.
function scenario_aliased_tree(; K::Int)
    V = JuMP.VariableRef
    return function ()
        model = Model()
        @variable(model, x)
        e = sin(x)
        for _ in 1:K
            e = JuMP.GenericNonlinearExpr{V}(:+, e, e)
        end
        @constraint(model, e <= 1)
        return model
    end
end

# Scenario B — many small INDEPENDENT NL constraints. No subexpression
# sharing within or across constraints. The model-level cache pays an
# insertion + hash cost for every new node; the per-call cache pays an
# extra Dict allocation per constraint; master pays nothing.
#
# This is the case blegat warned about in the PR thread: "we start
# creating many small dictionaries if there are a lot of constraints,
# but that's probably negligible". The benchmark tells us whether it
# really is negligible.
function scenario_many_independent(; N::Int = 5_000)
    return function ()
        model = Model()
        @variable(model, x[1:N])
        for i in 1:N
            @constraint(model, sin(x[i]) + cos(x[i]) * exp(x[i]) <= 1)
        end
        return model
    end
end

# Scenario C — many constraints sharing one big subexpression.
# `big = sum(sin(x[i]) for i in 1:N)` appears once in each of M constraints.
# The per-call cache cannot help here (each call sees `big` only once); only
# the model-level cache deduplicates `big` across constraints. This is the
# case that genuinely motivates the model-level cache.
function scenario_shared_big(; N::Int = 200, M::Int = 200)
    return function ()
        model = Model()
        @variable(model, x[1:N])
        big = sum(sin(x[i]) for i in 1:N)
        for j in 1:M
            @constraint(model, big * x[j] <= 1)
        end
        return model
    end
end

# Scenario D — many constraints, each with internal aliasing (a smaller
# binary tree per constraint, distinct leaf variable per constraint).
# Both per-call and model-level caches help WITHIN each constraint;
# across constraints there is nothing to share. This is the natural
# territory of odow's per-call dict.
function scenario_many_aliased(; M::Int = 200, K::Int = 8)
    V = JuMP.VariableRef
    return function ()
        model = Model()
        @variable(model, x[1:M])
        for i in 1:M
            e = sin(x[i])
            for _ in 1:K
                e = JuMP.GenericNonlinearExpr{V}(:+, e, e)
            end
            @constraint(model, e <= 1)
        end
        return model
    end
end

# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


const MODES = [:none, :per_call_oid, :model_struct, :model_oid]

function run_bench(label, build)
    println("="^72)
    println(label)
    println("="^72)
    # Warm up each mode (compile its method specializations) before
    # timing, so the first sample isn't biased by JIT.
    for m in MODES
        BENCH_MODE[] = m
        build()
    end
    results = Dict{Symbol,Any}()
    for m in MODES
        BENCH_MODE[] = m
        b = @benchmark $build() samples = 5 evals = 1 seconds = 30
        results[m] = b
        t_ms = minimum(b).time / 1e6
        mem_mb = minimum(b).memory / 1024^2
        allocs = minimum(b).allocs
        @printf("  %-13s  %10.2f ms   %9.2f MiB   %12d allocs\n",
                string(m), t_ms, mem_mb, allocs)
    end
    # Ratios versus the branch's model_struct baseline so we can see
    # where each strategy wins, ties, or loses.
    base = minimum(results[:model_struct]).time
    println()
    for m in MODES
        r = minimum(results[m]).time / base
        @printf("  ratio(%-13s / model_struct) = %.2f\n", string(m), r)
    end
    println()
end

function one_aliased_tree(; K = 16)
    run_bench("A: aliased tree (one big DAG, K=$K)",      scenario_aliased_tree(; K))
end

function independent()
    run_bench("B: many independent constraints (N=5000)", scenario_many_independent(N = 5_000))
end

function shared_subexpr()
    run_bench("C: shared big subexpr (N=200, M=200)",     scenario_shared_big(N = 200, M = 200))
end

function many_aliased_tree()
    run_bench("D: many aliased trees (M=200, K=8)",       scenario_many_aliased(M = 200, K = 8))
end

function run_all()
    one_aliased_tree()
    independent()
    shared_subexpr()
    many_aliased_tree()
end

#run_all()
