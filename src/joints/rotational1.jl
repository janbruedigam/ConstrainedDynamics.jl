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

    q1 = body1.state.qd[No]
    q2 = body2.state.qd[No]

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

    q2 = body2.state.qd[No]

    τ1 = vrotate(joint.V12' * -τ, joint.qoff) # in world coordinates
    τ2 = -τ  # in world coordinates
    
    τ2 = vrotate(τ2,inv(q2)) # in local coordinates

    F2 = @SVector zeros(T,3)

    updateForce!(joint, body2, F2, τ2, No)
    return
end


@inline function minimalCoordinates(joint::Rotational1, body1::Body, body2::Body, No)
    joint.V12 * (VLᵀmat(joint.qoff) * Lᵀmat(body1.state.qd[No]) * body2.state.qd[No])
end

@inline function minimalCoordinates(joint::Rotational1, body1::Origin, body2::Body, No)
    joint.V12 * (VLᵀmat(joint.qoff) * body2.state.qd[No])
end


@inline function g(joint::Rotational1, body1::Body, body2::Body, Δt)
    joint.V3 * g(joint, getx2(body1, Δt), getq2(body1, Δt), getx2(body2, Δt), getq2(body2, Δt))
end

@inline function g(joint::Rotational1, body1::Origin, body2::Body, Δt)
    joint.V3 * g(joint, getx2(body2, Δt), getq2(body2, Δt))
end


@inline function ∂g∂posa(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂posa(joint, getxd2(body1), getqd2(body1), getxd2(body2), getqd2(body2))
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂posb(joint, getxd2(body1), getqd2(body1), getxd2(body2), getqd2(body2))
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational1, body1::Origin, body2::Body)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂posb(joint, getxd2(body2), getqd2(body2))
    else
        return ∂g∂posb(joint)
    end
end


@inline function ∂g∂vela(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂vela(joint, getxd2(body1), getx2(body1, Δt), getqd2(body1), getq2(body1, Δt), getv2(body1), getω2(body1), getx2(body2, Δt), getq2(body2, Δt), Δt)
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂velb(joint, getx2(body1, Δt), getq2(body1, Δt), getxd2(body2), getx2(body2, Δt), getqd2(body2), getq2(body2, Δt), getv2(body2), getω2(body2), Δt)
    else
        return ∂g∂velb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1, body1::Origin, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂velb(joint, getxd2(body2), getx2(body2, Δt), getqd2(body2), getq2(body2, Δt), getv2(body2), getω2(body2), Δt)
    else
        return ∂g∂velb(joint)
    end
end
