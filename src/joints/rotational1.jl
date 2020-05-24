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


@inline function g(joint::Rotational1, body1::Body, body2::Body, Δt, No)
    joint.V3 * (VLᵀmat(joint.qoff) * Lᵀmat(getq2(body1, Δt)) * getq2(body2, Δt))
end

@inline function g(joint::Rotational1, body1::Origin, body2::Body, Δt, No)
    joint.V3 * VLᵀmat(joint.qoff) * getq2(body2, Δt)
end


@inline function ∂g∂posa(joint::Rotational1{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 1, 3)
        R = -joint.V3 * VLᵀmat(joint.qoff) * Rmat(body2.state.qd[No]) * RᵀVᵀmat(body1.state.qd[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational1{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 1, 3)
        R = joint.V3 * VLᵀmat(joint.qoff) * Lᵀmat(body1.state.qd[No]) * LVᵀmat(body2.state.qd[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational1{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 1, 3)
        Ω = joint.V3 * VLᵀmat(joint.qoff) * Rmat(ωbar(getω2(body2), Δt)) * Rmat(body2.state.qd[No]) * Rᵀmat(body1.state.qd[No]) * Tmat(T) * derivωbar(getω2(body1), Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 1, 3)
        Ω = joint.V3 * VLᵀmat(joint.qoff) * Lᵀmat(ωbar(getω2(body1), Δt)) * Lᵀmat(body1.state.qd[No]) * Lmat(body2.state.qd[No]) * derivωbar(getω2(body2), Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function ∂g∂posb(joint::Rotational1{T}, body1::Origin, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 1, 3)
        R = joint.V3 * VLᵀmat(joint.qoff) * LVᵀmat(body2.state.qd[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1{T}, body1::Origin, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 1, 3)
        Ω = joint.V3 * VLᵀmat(joint.qoff) * Lmat(body2.state.qd[No]) * derivωbar(getω2(body2), Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
