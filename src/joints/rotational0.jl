@inline function getPositionDelta(joint::Rotational0, body1::AbstractBody, body2::Body{T}, θ::SVector{3,T}) where T
    # axis angle representation
    if norm(θ) == 0
        q = one(UnitQuaternion{T})
    else
        q = UnitQuaternion(cos(norm(θ)/2),(θ/norm(θ)*sin(norm(θ)/2))..., false)
    end
    
    Δq = joint.qoff * q # in body1 frame
    return Δq
end
@inline function getVelocityDelta(joint::Rotational0, body1::Body, body2::Body{T}, ω::SVector{3,T}) where T
    Δω = vrotate(ω, inv(body2.state.qc)*body1.state.qc*joint.qoff) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational0, body1::Origin, body2::Body{T}, ω::SVector{3,T}) where T
    Δω = vrotate(ω, inv(body2.state.qc)*joint.qoff) # in body2 frame
    return Δω
end

@inline function setForce!(joint::Rotational0, body1::Body, body2::Body{T}, τ::SVector{3,T}) where T
    setForce!(joint, body1.state, body2.state, τ)
    return
end
@inline function setForce!(joint::Rotational0, body1::Origin, body2::Body{T}, τ::SVector{3,T}) where T
    setForce!(joint, body2.state, τ)
    return
end

@inline function minimalCoordinates(joint::Rotational0, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    q = g(joint, statea.qc, stateb.qc)
    return rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational0, body1::Origin, body2::Body)
    stateb = body2.state
    q = g(joint, stateb.qc)
    return rotation_vector(q)
end

@inline g(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt) = g(joint)

@inline ∂g∂posa(joint::Rotational0, body1::AbstractBody, body2::AbstractBody) = ∂g∂posa(joint)
@inline ∂g∂posb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody) = ∂g∂posb(joint)
@inline ∂g∂vela(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt) = ∂g∂vela(joint)
@inline ∂g∂velb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt) = ∂g∂velb(joint)
