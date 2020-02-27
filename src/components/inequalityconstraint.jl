mutable struct InequalityConstraint{T,N,Cs} <: AbstractConstraint{T,N}
    id::Int64

    constraints::Cs
    pid::Int64
    # bodyid::Int64

    s0::SVector{N,T}
    s1::SVector{N,T}
    γ0::SVector{N,T}
    γ1::SVector{N,T}

    μ::T
    α::T
    τ::T

    function InequalityConstraint(data...)
        bounddata = Tuple{Bound,Int64}[]
        for info in data
            if typeof(info[1]) <: Bound
                push!(bounddata, info)
            else
                for subinfo in info
                    push!(bounddata, subinfo)
                end
            end
        end

        T = getT(bounddata[1][1])

        pid = bounddata[1][2]
        bodyids = Int64[]
        constraints = Bound{T}[]
        N = 0
        for set in bounddata
            push!(constraints, set[1])
            @assert set[2] == pid
            N += 1 # getNc(set[1])
        end
        constraints = Tuple(constraints)
        # Nc = length(constraints)

        s0 = ones(T, N)
        s1 = ones(T, N)
        γ0 = ones(T, N)
        γ1 = ones(T, N)

        α = 1
        μ = 1
        τ = 0.995

        new{T,N,typeof(constraints)}(getGlobalID(), constraints, pid, s0, s1, γ0, γ1, α, μ, τ)
    end
end


Base.length(::InequalityConstraint{T,N}) where {T,N} = N

function resetVars!(ineqc::InequalityConstraint{T,N}) where {T,N}
    ineqc.s0 = @SVector ones(T, N)
    ineqc.s1 = @SVector ones(T, N)
    ineqc.γ0 = @SVector ones(T, N)
    ineqc.γ1 = @SVector ones(T, N)

    ineqc.α = 1
    ineqc.μ = 1
    return 
end

function feasibilityStepLenth(ineqc::InequalityConstraint{T,N}, ineqentry) where {T,N}
    s1 = ineqc.s1
    γ1 = ineqc.γ1
    Δs = ineqentry.Δs
    Δγ = ineqentry.Δγ
    τ = ineqc.τ

    for i = 1:N
        αmax = τ * s1[i] / Δs[i]
        (αmax > 0) && (αmax < ineqc.α) && (ineqc.α = αmax)
        αmax = τ * γ1[i] / Δγ[i]
        (αmax > 0) && (αmax < ineqc.α) && (ineqc.α = αmax)
    end

    return ineqc.α
end


function g(ineqc::InequalityConstraint{T,1}, mechanism) where {T}
    g(ineqc.constraints[1], getbody(mechanism, ineqc.pid), mechanism.dt, mechanism.No)
end

@generated function g(ineqc::InequalityConstraint{T,N}, mechanism) where {T,N}
    vec = [:(g(ineqc.constraints[$i], getbody(mechanism, ineqc.pid), mechanism.dt, mechanism.No)) for i = 1:N]
    :(SVector{N,T}($(vec...)))
end

function gs(ineqc::InequalityConstraint{T,1}, mechanism) where {T}
    g(ineqc.constraints[1], getbody(mechanism, ineqc.pid), mechanism.dt, mechanism.No) - ineqc.s1[1]
end

@generated function gs(ineqc::InequalityConstraint{T,N}, mechanism) where {T,N}
    vec = [:(g(ineqc.constraints[$i], getbody(mechanism, ineqc.pid), mechanism.dt, mechanism.No) - ineqc.s1[$i]) for i = 1:N]
    :(SVector{N,T}($(vec...)))
end

function h(ineqc::InequalityConstraint)
    ineqc.s1 .* ineqc.γ1
end

function hμ(ineqc::InequalityConstraint{T}) where T
    ineqc.s1 .* ineqc.γ1 .- ineqc.μ
end


function schurf(ineqc::InequalityConstraint{T,N}, body, mechanism) where {T,N}
    val = @SVector zeros(T, 6)
    for i = 1:N
        val += schurf(ineqc, ineqc.constraints[i], i, body, mechanism.dt, mechanism.No)
    end
    return val
end

function schurD(ineqc::InequalityConstraint{T,N}, body, dt) where {T,N}
    val = @SMatrix zeros(T, 6, 6)
    for i = 1:N
        val += schurD(ineqc, ineqc.constraints[i], i, body, dt)
    end
    return val
end

@generated function ∂g∂pos(ineqc::InequalityConstraint{T,N}, body, mechanism) where {T,N}
    vec = [:(∂g∂pos(ineqc.constraints[$i], mechanism.No)) for i = 1:N]
    :(vcat($(vec...)))
end

@generated function ∂g∂vel(ineqc::InequalityConstraint{T,N}, body, mechanism) where {T,N}
    vec = [:(∂g∂vel(ineqc.constraints[$i], mechanism.dt, mechanism.No)) for i = 1:N]
    :(vcat($(vec...)))
end
