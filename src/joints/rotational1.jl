# No idea what kind of joint this actually is...
@inline function getVelocityDelta(joint::Rotational1, body1::AbstractBody, body2::Body{T}, ω) where T
    #TODO define this function
    @error("not defined for rot2")
end

@inline function getPositionDelta(joint::Rotational1, body1::AbstractBody, body2::Body{T}, θ) where T
    #TODO define this function
    @error("not defined for rot2")
end

@inline function setForce!(joint::Rotational1, body1::Body, body2::Body{T}, τ::SVector{2,T}, No) where T
    τ1 = vrotate(joint.V12' * -τ, body1.q[No] * joint.qoff)
    τ2 = -τ1

    body1.τ[No] = τ1
    body2.τ[No] = τ2
    return
end

@inline function setForce!(joint::Rotational1, body1::Origin, body2::Body{T}, τ::SVector{2,T}, No) where T
    body2.τ[No] = vrotate(joint.V12' * τ, joint.qoff)
    return
end


@inline function minimalCoordinates(joint::Rotational1, body1::Body, body2::Body, No)
    joint.V12 * (VLᵀmat(joint.qoff) * Lᵀmat(body1.q[No]) * body2.q[No])
end

@inline function minimalCoordinates(joint::Rotational1, body1::Origin, body2::Body, No)
    joint.V12 * (VLᵀmat(joint.qoff) * body2.q[No])
end


@inline function g(joint::Rotational1, body1::Body, body2::Body, Δt, No)
    joint.V3 * (VLᵀmat(joint.qoff) * Lᵀmat(getq3(body1, Δt)) * getq3(body2, Δt))
end

@inline function g(joint::Rotational1, body1::Origin, body2::Body, Δt, No)
    joint.V3 * VLᵀmat(joint.qoff) * getq3(body2, Δt)
end


@inline function ∂g∂posa(joint::Rotational1{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 1, 3)
        R = -joint.V3 * VLᵀmat(joint.qoff) * Rmat(body2.q[No]) * RᵀVᵀmat(body1.q[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational1{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 1, 3)
        R = joint.V3 * VLᵀmat(joint.qoff) * Lᵀmat(body1.q[No]) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational1{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 1, 3)
        Ω = joint.V3 * VLᵀmat(joint.qoff) * Rmat(ωbar(body2, Δt)) * Rmat(body2.q[No]) * Rᵀmat(body1.q[No]) * Tmat(T) * derivωbar(body1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 1, 3)
        Ω = joint.V3 * VLᵀmat(joint.qoff) * Lᵀmat(ωbar(body1, Δt)) * Lᵀmat(body1.q[No]) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function ∂g∂posb(joint::Rotational1{T}, body1::Origin, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 1, 3)
        R = joint.V3 * VLᵀmat(joint.qoff) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1{T}, body1::Origin, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 1, 3)
        Ω = joint.V3 * VLᵀmat(joint.qoff) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
