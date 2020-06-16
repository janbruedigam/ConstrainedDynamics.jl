abstract type Joint{T,N} end

Base.show(io::IO, joint::Joint) = summary(io, joint)

getT(joint::Joint{T}) where T = T
getN(joint::Joint{T,N}) where {T,N} = N

@inline setForce!(joint::Joint) = return
@inline setForce!(joint::Joint, body1::AbstractBody, body2::AbstractBody, Fτ::SVector{0,T}) where T = setForce!(joint)

@inline minimalCoordinates(joint::Joint{T,N}) where {T,N} = szeros(T, 3 - N)
@inline g(joint::Joint{T,N}) where {T,N} = szeros(T, N)

@inline ∂g∂ᵣposa(joint::Joint{T,N}) where {T,N} = szeros(T, N, 6)
@inline ∂g∂ᵣposb(joint::Joint{T,N}) where {T,N} = szeros(T, N, 6)
@inline ∂g∂ᵣvela(joint::Joint{T,N}) where {T,N} = szeros(T, N, 6)
@inline ∂g∂ᵣvelb(joint::Joint{T,N}) where {T,N} = szeros(T, N, 6)

@inline ∂Fτ∂ub(joint::Joint{T,N}) where {T,N} = szeros(T, 6, 3 - N)


# Derivatives accounting for quaternion specialness

@inline function ∂g∂ᵣposa(joint::Joint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posa(joint, xa, qa, xb, qb)
    Q = Q * LVᵀmat(qa)

    return [X Q]
end
@inline function ∂g∂ᵣposb(joint::Joint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posb(joint, xa, qa, xb, qb)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end
@inline function ∂g∂ᵣposb(joint::Joint, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posb(joint, xb, qb)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end

@inline function ∂g∂ᵣvela(joint::Joint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
        x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector, Δt
    )

    X, Q = ∂g∂posa(joint, x2a, q2a, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1a) * derivωbar(ω1a, Δt)

    return [V Ω]
end
@inline function ∂g∂ᵣvelb(joint::Joint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
    x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt
    )

    X, Q = ∂g∂posb(joint, x2a, q2a, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt)

    return [V Ω]
end
@inline function ∂g∂ᵣvelb(joint::Joint, x2b::AbstractVector, q2b::UnitQuaternion,
        x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt
    )

    X, Q = ∂g∂posb(joint, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt)

    return [V Ω]
end


# Convenience calls for dynamics calculation solver
g(joint::Joint, statea::State, stateb::State, Δt) = g(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...)
g(joint::Joint, stateb::State, Δt) = g(joint, posargsnext(stateb, Δt)...)

∂g∂ᵣposa(joint::Joint, statea::State, stateb::State) = ∂g∂ᵣposa(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ᵣposb(joint::Joint, statea::State, stateb::State) = ∂g∂ᵣposb(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ᵣposb(joint::Joint, stateb::State) = ∂g∂ᵣposb(joint, posargsk(stateb)...)

∂g∂ᵣvela(joint::Joint, statea::State, stateb::State, Δt) = ∂g∂ᵣvela(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(statea)..., Δt)
∂g∂ᵣvelb(joint::Joint, statea::State, stateb::State, Δt) = ∂g∂ᵣvelb(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(stateb)..., Δt)
∂g∂ᵣvelb(joint::Joint, stateb::State, Δt) = ∂g∂ᵣvelb(joint, posargsnext(stateb, Δt)..., fullargssol(stateb)..., Δt)
