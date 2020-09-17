abstract type Joint{T,N} end

Joint0 = Joint{T,0} where T
Joint1 = Joint{T,1} where T
Joint2 = Joint{T,2} where T
Joint3 = Joint{T,3} where T

Base.show(io::IO, joint::Joint) = summary(io, joint)

getT(joint::Joint{T}) where T = T
getN(joint::Joint{T,N}) where {T,N} = N

@inline setForce!(joint::Joint) = return
@inline setForce!(joint::Joint, body1::AbstractBody, body2::AbstractBody, Fτ::SVector{0}) = setForce!(joint)

@inline minimalCoordinates(joint::Joint{T,N}) where {T,N} = szeros(T, 3 - N)
@inline g(joint::Joint{T,N}) where {T,N} = szeros(T, N)

@inline ∂g∂ʳposa(joint::Joint{T,N}) where {T,N} = szeros(T, N, 6)
@inline ∂g∂ʳposb(joint::Joint{T,N}) where {T,N} = szeros(T, N, 6)
@inline ∂g∂ʳvela(joint::Joint{T,N}) where {T,N} = szeros(T, N, 6)
@inline ∂g∂ʳvelb(joint::Joint{T,N}) where {T,N} = szeros(T, N, 6)

@inline ∂Fτ∂ub(joint::Joint{T,N}) where {T,N} = szeros(T, 6, 3 - N)


# Derivatives accounting for quaternion specialness

@inline function ∂g∂ʳposa(joint::Joint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posa(joint, xa, qa, xb, qb)
    Q = Q * LVᵀmat(qa)

    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Joint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posb(joint, xa, qa, xb, qb)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Joint, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posb(joint, xb, qb)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end

@inline function ∂g∂ʳvela(joint::Joint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
        x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector, Δt
    )

    X, Q = ∂g∂posa(joint, x2a, q2a, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1a) * derivωbar(ω1a, Δt)

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::Joint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
    x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt
    )

    X, Q = ∂g∂posb(joint, x2a, q2a, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt)

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::Joint, x2b::AbstractVector, q2b::UnitQuaternion,
        x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt
    )

    X, Q = ∂g∂posb(joint, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt)

    return [V Ω]
end


# Calls for dynamics
g(joint::Joint, statea::State, stateb::State, Δt) = g(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...)
g(joint::Joint, stateb::State, Δt) = g(joint, posargsnext(stateb, Δt)...)

∂g∂ʳposa(joint::Joint, statea::State, stateb::State) = ∂g∂ʳposa(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::Joint, statea::State, stateb::State) = ∂g∂ʳposb(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::Joint, stateb::State) = ∂g∂ʳposb(joint, posargsk(stateb)...)

∂g∂ʳvela(joint::Joint, statea::State, stateb::State, Δt) = ∂g∂ʳvela(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(statea)..., Δt)
∂g∂ʳvelb(joint::Joint, statea::State, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(stateb)..., Δt)
∂g∂ʳvelb(joint::Joint, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsnext(stateb, Δt)..., fullargssol(stateb)..., Δt)
