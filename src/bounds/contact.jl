abstract type Contact{T} <: Bound{T} end

### Constraints and derivatives
## Discrete-time position wrappers (for dynamics)
@inline g(contact::Contact, state::State, Δt) = g(contact, posargsnext(state, Δt)...)

## Position level constraints (for dynamics)
# @inline g(contact::Contact{T}) where {T} = szeros(T,6)
@inline g(contact::Contact, x::AbstractVector, q::UnitQuaternion) = contact.ainv3 * (x + vrotate(contact.p,q) - contact.offset)


## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂pos(contact::Contact, x::AbstractVector, q::UnitQuaternion)
    p = contact.p
    X = contact.ainv3
    Q = contact.ainv3 * (VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p)))
    return X, Q
end

## Discrete-time position derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳpos(contact::Contact, body::Body)
    return ∂g∂ʳpos(contact, body.state)
end

#Wrappers 2
@inline ∂g∂ʳpos(contact::Contact{T}) where T = szeros(T, 6)
∂g∂ʳpos(contact::Contact, state::State) = ∂g∂ʳpos(contact, posargsk(state)...)

#Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳpos(contact::Contact, x::AbstractVector, q::UnitQuaternion)
    X, Q = ∂g∂pos(contact, x, q)
    Q = Q * LVᵀmat(q)

    return [X Q]
end

## Discrete-time velocity derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳvel(contact::Contact, body::Body, Δt)
    return ∂g∂ʳvel(contact, body.state, Δt)
end

#Wrappers 2
@inline ∂g∂ʳvel(contact::Contact{T}) where {T} = szeros(T, 6)
∂g∂ʳvel(contact::Contact, state::State, Δt) = ∂g∂ʳvel(contact, posargsnext(state, Δt)..., fullargssol(state)..., Δt)

#Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳvel(contact::Contact,x2::AbstractVector, q2::UnitQuaternion, x1::AbstractVector,v1::AbstractVector, q1::UnitQuaternion,w1::AbstractVector, Δt)
    X, Q = ∂g∂pos(contact, x2, q2)
    V = X * Δt
    Ω = Q * Lmat(q1) * derivωbar(w1, Δt) * Δt / 2
    return [V Ω]
end

##Schurf und SchurD
#Schurf
@inline function schurf(ineqc, contact::Contact, i, body::Body, μ, Δt)
    return schurf(contact, body.state, ineqc.γsol[2][i], ineqc.ssol[2][i], μ, Δt)[SA[1; 2; 3; 4; 5; 6]]
end
@inline function schurf(contact::Contact, state::State, γ, s, μ, Δt)
    φ = g(contact, posargsnext(state, Δt)...)
    return ∂g∂ʳpos(contact, state)' * (γ / s * φ - μ / s)
end

#SchurD
@inline function schurD(ineqc, contact::Contact, i, body::Body, Δt)
    return schurD(contact, body.state, ineqc.γsol[2][i], ineqc.ssol[2][i], Δt)
end
@inline function schurD(contact::Contact,state::State, γ, s, Δt)
    return ∂g∂ʳpos(contact, state)' * γ/s * ∂g∂ʳvel(contact, state, Δt)
end


## Additional force for friction
@inline additionalforce(bound::Contact, state::State) = additionalforce(bound, posargsk(state)...)
@inline additionalforce(bound::Contact, body::Body) = additionalforce(bound, body.state)