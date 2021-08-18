abstract type Bound{T,N} end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Bound{T,N}) where {T,N}
    summary(io, constraint)
end

### General functions
getT(bound::Bound{T}) where T = T
Base.length(bound::Bound{T,N}) where {T,N} = N


### Constraints and derivatives
## Position level constraint wrappers
g(bound::Bound, body::Body, Δt) = g(bound, body.state, Δt)

## Discrete-time position wrappers (for dynamics)
@inline g(bound::Bound, state::State, Δt) = g(bound, posargsnext(state, Δt)...)

## Discrete-time position derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳpos(bound::Bound, body::Body)
    return ∂g∂ʳpos(bound, body.state)
end

# Wrappers 2
∂g∂ʳpos(bound::Bound, state::State) = ∂g∂ʳpos(bound, posargsk(state)...)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳpos(bound::Bound, x::AbstractVector, q::UnitQuaternion)
    X, Q = ∂g∂pos(bound, x, q)
    Q = Q * LVᵀmat(q)

    return [X Q]
end

## Discrete-time velocity derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳvel(bound::Bound, body::Body, Δt)
    return ∂g∂ʳvel(bound, body.state, Δt)
end

# Wrappers 2
∂g∂ʳvel(bound::Bound, state::State, Δt) = ∂g∂ʳvel(bound, posargsnext(state, Δt)..., fullargssol(state)..., Δt)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳvel(bound::Bound, x2::AbstractVector, q2::UnitQuaternion,
        x1::AbstractVector, v1::AbstractVector, q1::UnitQuaternion, ω1::AbstractVector, Δt
    )

    X, Q = ∂g∂pos(bound, x2, q2)
    V = X * Δt
    Ω = Q * Lmat(q1) * derivωbar(ω1, Δt) * Δt / 2

    return [V Ω]
end