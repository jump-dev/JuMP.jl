module JuMPEnzymeExt

using Enzyme
using JuMP

function jump_operator(f::Function)
    @inline function f!(y, x...)
        y[1] = f(x...)
    end
    function gradient!(g::AbstractVector{T}, x::Vararg{T,N}) where {T,N}
        y = zeros(T,1)
        ry = ones(T,1)
        rx = ntuple(N) do i
            Active(x[i])
        end
        g .= autodiff(Reverse, f!, Const, Duplicated(y,ry), rx...)[1][2:end]
        return nothing
    end

    function gradient_deferred!(g, y, ry, rx...)
        g .= autodiff_deferred(Reverse, f!, Const, Duplicated(y,ry), rx...)[1][2:end]
        return nothing
    end

    function hessian!(H::AbstractMatrix{T}, x::Vararg{T,N}) where {T,N}
        y = zeros(T,1)
        dy = ntuple(N) do i
            ones(1)
        end
        g = zeros(T,N)
        dg = ntuple(N) do i
            zeros(T,N)
        end
        ry = ones(1)
        dry = ntuple(N) do i
            zeros(T,1)
        end
        rx = ntuple(N) do i
            Active(x[i])
        end

        args = ntuple(N) do i
            drx = ntuple(N) do j
                if i == j
                    Active(one(T))
                else
                    Active(zero(T))
                end
            end
            BatchDuplicated(rx[i], drx)
        end
        autodiff(Forward, gradient_deferred!, Const, BatchDuplicated(g,dg), BatchDuplicated(y,dy), BatchDuplicated(ry, dry), args...)
        for i in 1:N
            for j in 1:N
                if i <= j
                    H[j,i] = dg[j][i]
                end
            end
        end
        return nothing
    end

    return gradient!, hessian!
end

function JuMP.add_nonlinear_operator(
    model::GenericModel,
    dim::Int,
    f::Function;
    name::Symbol = Symbol(f),
)
    gradient, hessian = jump_operator(f)
    MOI.set(model, MOI.UserDefinedFunction(name, dim), tuple(f, gradient, hessian))
    return NonlinearOperator(f, name)
end
end
