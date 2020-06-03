abstract type Contact{T} <: Bound{T} end

@inline g(contact::Contact, body::Body, Δt) = g(contact, body.state, Δt)

@inline ∂g∂pos(contact::Contact, body::Body) = ∂g∂pos(contact, body.state)
@inline ∂g∂vel(contact::Contact, body::Body, Δt) = ∂g∂vel(contact, body.state, Δt)

@inline function schurf(ineqc, contact::Contact, i, body::Body, μ, Δt)
    schurf(contact, body.state, ineqc.γsol[2][i], ineqc.ssol[2][i], μ, Δt)
end

@inline function schurD(ineqc, contact::Contact, i, body::Body, Δt)
    schurD(contact, body.state, ineqc.γsol[2][i], ineqc.ssol[2][i], Δt)
end