@inline function getPositionDelta(joint::Rotational2, body1::AbstractBody, body2::Body{T}, θ::SVector{1,T}) where T
    q = Quaternion(cos(θ[1]/2),(joint.V3*sin(θ[1]/2))...)
    Δq = joint.qoff * q # in body1 frame
    return Δq
end
@inline function getVelocityDelta(joint::Rotational2, body1::Body, body2::Body{T}, ω::SVector{1,T}) where T
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qc)*body1.state.qc*joint.qoff) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational2, body1::Origin, body2::Body{T}, ω::SVector{1,T}) where T
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qc)*joint.qoff) # in body2 frame
    return Δω
end

@inline function setForce!(joint::Rotational2, body1::Body, body2::Body{T}, τ::SVector{1,T}) where T
    setForce!(joint, body1.state, body2.state, joint.V3' * τ)
    return
end
@inline function setForce!(joint::Rotational2, body1::Origin, body2::Body{T}, τ::SVector{1,T}) where T
    setForce!(joint, body2.state, joint.V3' * τ)
    return
end

@inline function minimalCoordinates(joint::Rotational2, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    q = g(joint, statea.qc, stateb.qc)
    joint.V3 * axis(q) * angle(q) 
end
@inline function minimalCoordinates(joint::Rotational2, body1::Origin, body2::Body)
    stateb = body2.state
    q = g(joint, stateb.qc)
    joint.V3 * axis(q) * angle(q)
end

@inline g(joint::Rotational2, body1::Body, body2::Body, Δt) = joint.V12 * g(joint, body1.state, body2.state, Δt)
@inline g(joint::Rotational2, body1::Origin, body2::Body, Δt) = joint.V12 * g(joint, body2.state, Δt)

@inline function ∂g∂posa(joint::Rotational2, body1::Body, body2::Body)
    if body2.id == joint.cid
        return joint.V12 * ∂g∂posa(joint, body1.state, body2.state)
    else
        return ∂g∂posa(joint)
    end
end
@inline function ∂g∂posb(joint::Rotational2, body1::Body, body2::Body)
    if body2.id == joint.cid
        return joint.V12 * ∂g∂posb(joint, body1.state, body2.state)
    else
        return ∂g∂posb(joint)
    end
end
@inline function ∂g∂posb(joint::Rotational2, body1::Origin, body2::Body)
    if body2.id == joint.cid
        return joint.V12 * ∂g∂posb(joint, body2.state)
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational2, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V12 * ∂g∂vela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂vela(joint)
    end
end
@inline function ∂g∂velb(joint::Rotational2, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V12 * ∂g∂velb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end
@inline function ∂g∂velb(joint::Rotational2, body1::Origin, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V12 * ∂g∂velb(joint, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end