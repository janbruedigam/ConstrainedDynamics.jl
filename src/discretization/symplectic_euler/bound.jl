@inline g(contact::Contact, x, q) = contact.Nx * (x + vrotate(contact.p,q) - contact.offset)
@inline g(contact::Contact, state::State, Δt) = g(contact, getx3(state, Δt), getq3(state, Δt))

##Discrete Time Position Derivatives
@inline function ∂g∂pos(contact::Contact, x::AbstractVector, q::UnitQuaternion)
    p = contact.p
    X = contact.Nx
    Q = contact.Nx * (VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p)))
    return X, Q
end


##Schurf und SchurD
#Schurf
@inline function schurf(contact::Contact, state::State, γ, s, μ, Δt)
    φ = g(contact, getx3(state, Δt), getq3(state, Δt))
    return ∂g∂ʳpos(contact, state)' * (γ / s * φ - μ / s)
end

#SchurD
@inline function schurD(contact::Contact,state::State, γ, s, Δt)
    return ∂g∂ʳpos(contact, state)' * γ/s * ∂g∂ʳvel(contact, state, Δt)
end


##Friction
@inline function calcFrictionForce!(mechanism::Mechanism{T}, friction::Contact, body::Body, γ) where T
    cf = friction.cf
    D = friction.D
    state = body.state

    d = state.d
    v = state.vsol[2]
    ω = state.ωsol[2]
    state.vsol[2] = szeros(T, 3)
    state.ωsol[2] = szeros(T, 3)
    dyn = dynamics(mechanism, body)
    state.vsol[2] = v
    ω = state.ωsol[2] = ω
    state.d = d

    b0 = D*dyn + friction.b # remove old friction force

    if norm(b0) > cf*γ
        friction.b = b0/norm(b0)*cf*γ
    else
        friction.b = b0
    end    
    return
end