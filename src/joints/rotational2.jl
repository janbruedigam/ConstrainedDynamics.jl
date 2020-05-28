@inline function getPositionDelta(joint::Rotational2, body1::AbstractBody, body2::Body{T}, θ::SVector{1,T}) where T
    q = Quaternion(cos(θ[1]/2),(joint.V3*sin(θ[1]/2))...)
    Δq = joint.qoff * q # in body1 frame
    return Δq
end

@inline function getVelocityDelta(joint::Rotational2, body1::Body, body2::Body{T}, ω::SVector{1,T}) where T
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qk[2])*body1.state.qk[2]*joint.qoff) # in body2 frame
    return Δω
end

@inline function getVelocityDelta(joint::Rotational2, body1::Origin, body2::Body{T}, ω::SVector{1,T}) where T
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qk[2])*joint.qoff) # in body2 frame
    return Δω
end

@inline function setForce!(joint::Rotational2, body1::Body, body2::Body{T}, τ::SVector{1,T}, No) where T
    clearForce!(joint, body1, body2, No)

    q1 = body1.state.qk[No]
    q2 = body2.state.qk[No]

    τ1 = vrotate(joint.V3' * -τ, q1*joint.qoff) # in world coordinates
    τ2 = -τ1 # in world coordinates

    τ1 = vrotate(τ1,inv(q1)) # in local coordinates
    τ2 = vrotate(τ2,inv(q2)) # in local coordinates

    F1 =  @SVector zeros(T,3)
    F2 =  @SVector zeros(T,3)

    updateForce!(joint, body1, body2, F1, τ1, F2, τ2, No)
    return
end

@inline function setForce!(joint::Rotational2, body1::Origin, body2::Body{T}, τ::SVector{1,T}, No) where T
    clearForce!(joint, body2, No)

    q2 = body2.state.qk[No]

    τ1 = vrotate(joint.V3' * -τ, joint.qoff) # in world coordinates
    τ2 = -τ1 # in world coordinates

    τ2 = vrotate(τ2,inv(q2)) # in local coordinates

    F2 = @SVector zeros(T,3)

    updateForce!(joint, body2, F2, τ2, No)
    return
end


@inline function minimalCoordinates(joint::Rotational2, body1::Origin, body2::Body, No)
    q2 = joint.qoff \ body2.state.qk[No]
    joint.V3 * axis(q2) * angle(q2) 
end

@inline function minimalCoordinates(joint::Rotational2, body1::Body, body2::Body, No)
    q2 = joint.qoff \ (body1.state.qk[No] \ body2.state.qk[No])
    joint.V3 * axis(q2) * angle(q2)
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