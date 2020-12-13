abstract type Contact{T} <: Bound{T} end

###Constraints and Derivatives
##Position level constraint wrappers
g(contact::Contact, body::Body, Δt) = g(contact, body.state, Δt)
@inline g(contact::Contact{T}) where {T} = szeros(T,6)


##Discrete Time Position Derivatives (for dynamics)
#Wrappers 1
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

##Discrete Time Velocity Derivatives 
#Wrappers 1
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
    Ω = Q * Lmat(q1) * derivωbar(w1, Δt)
    return [V Ω]
end


#Schurf
@inline function schurf(ineqc, contact::Contact, i, body::Body, μ, Δt)
    return schurf(contact, body.state, ineqc.γsol[2][i], ineqc.ssol[2][i], μ, Δt)[1:6]
end


#SchurD
@inline function schurD(ineqc, contact::Contact, i, body::Body, Δt)
    return schurD(contact, body.state, ineqc.γsol[2][i], ineqc.ssol[2][i], Δt)
end

