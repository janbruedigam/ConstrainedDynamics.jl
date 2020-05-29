# No idea what kind of joint this actually is...
@inline function getPositionDelta(joint::Rotational1, body1::AbstractBody, body2::Body{T}, θ) where T
    #TODO define this function
    @error("not defined for rot2")
end

@inline function getVelocityDelta(joint::Rotational1, body1::AbstractBody, body2::Body{T}, ω) where T
    #TODO define this function
    @error("not defined for rot2")
end

@inline function setForce!(joint::Rotational1, body1::Body, body2::Body{T}, τ::SVector{2,T}, No) where T
    clearForce!(joint, body1, body2, No)

    q1 = body1.state.qk[No]
    q2 = body2.state.qk[No]

    τ1 = vrotate(joint.V12' * -τ, q1*joint.qoff) # in world coordinates
    τ2 = -τ1 # in world coordinates

    τ1 = vrotate(τ1,inv(q1)) # in local coordinates
    τ2 = vrotate(τ2,inv(q2)) # in local coordinates

    F1 =  @SVector zeros(T,3)
    F2 =  @SVector zeros(T,3)

    updateForce!(joint, body1, body2, F1, τ1, F2, τ2, No)
    return
end

@inline function setForce!(joint::Rotational1, body1::Origin, body2::Body{T}, τ::SVector{2,T}, No) where T
    clearForce!(joint, body2, No)

    q2 = body2.state.qk[No]

    τ1 = vrotate(joint.V12' * -τ, joint.qoff) # in world coordinates
    τ2 = -τ  # in world coordinates
    
    τ2 = vrotate(τ2,inv(q2)) # in local coordinates

    F2 = @SVector zeros(T,3)

    updateForce!(joint, body2, F2, τ2, No)
    return
end

@inline function minimalCoordinates(joint::Rotational1, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    q = g(joint, statea.qc, stateb.qc)
    joint.V12 * axis(q) * angle(q) 
end
@inline function minimalCoordinates(joint::Rotational1, body1::Origin, body2::Body)
    stateb = body2.state
    q = g(joint, stateb.qc)
    joint.V12 * axis(q) * angle(q)
end

@inline g(joint::Rotational1, body1::Body, body2::Body, Δt) = joint.V3 * g(joint, body1.state, body2.state, Δt)
@inline g(joint::Rotational1, body1::Origin, body2::Body, Δt) = joint.V3 * g(joint, body2.state, Δt)


@inline function ∂g∂posa(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂posa(joint, body1.state, body2.state)
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂posb(joint, body1.state, body2.state)
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational1, body1::Origin, body2::Body)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂posb(joint, body2.state)
    else
        return ∂g∂posb(joint)
    end
end


@inline function ∂g∂vela(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂vela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂velb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1, body1::Origin, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂velb(joint, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end
