@inline g(contact::Contact, x) = contact.Nx[SA[1; 2; 3]]' * (x - contact.offset[SA[1; 2; 3]])
@inline g(contact::Contact, state::State, Δt) = g(contact, getx3(state, Δt))

@inline ∂g∂pos(contact::Contact) = contact.Nx
@inline ∂g∂pos(contact::Contact, state::State) = ∂g∂pos(contact)

@inline ∂g∂vel(contact::Contact, Δt) = contact.Nx * Δt
@inline ∂g∂vel(contact::Contact, state::State, Δt) = ∂g∂vel(contact, Δt)

@inline schurf(contact::Contact, φ, γ, s, μ) = ∂g∂pos(contact)' * (γ / s * φ - μ / s)
@inline function schurf(contact::Contact, state::State, γ, s, μ, Δt)
    φ = g(contact, getx3(state, Δt))
    return schurf(contact, φ, γ, s, μ)
end

@inline function schurD(contact::Contact, γ, s, Δt)
    Nx = ∂g∂pos(contact)
    Nv = ∂g∂vel(contact, Δt)

    return Nx' * γ / s * Nv
end
@inline schurD(contact::Contact, state::State, γ, s, Δt) = schurD(contact, γ, s, Δt)


@inline function calcFrictionForce!(mechanism, friction::Contact, body::Body, γ)
    cf = friction.cf
    D = friction.D
    state = body.state

    d = state.d
    v = state.vsol[2]
    ω = state.ωsol[2]
    state.vsol[2] = @SVector zeros(3)
    state.ωsol[2] = @SVector zeros(3)
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