@inline function getVelocityDelta(joint::Rotational2, body1::AbstractBody, body2::Body{T}, ω::Union{T,SVector{1,T}}) where T
    ω = joint.V3' * ω
    Δω = vrotate(ω, joint.qoff) # in body1 frame
    return Δω
end

@inline function getPositionDelta(joint::Rotational2, body1::AbstractBody, body2::Body{T}, θ::Union{T,SVector{1,T}}) where T
    q = Quaternion(cos(θ/2),(joint.V3*sin(θ/2))...)
    Δq = joint.qoff * q # in body1 frame
    return Δq
end

@inline function setForce!(joint::Rotational2, body1::Body, body2::Body{T}, τ::Union{T,SVector{1,T}}, No) where T
    τ1 = vrotate(joint.V3' * -τ, body1.q[No] * joint.qoff)
    τ2 = -τ1

    body1.τ[No] = τ1
    body2.τ[No] = τ2
    return
end

@inline function setForce!(joint::Rotational2, body1::Origin, body2::Body{T}, τ::Union{T,SVector{1,T}}, No) where T
    body2.τ[No] = vrotate(joint.V3' * τ, joint.qoff)
    return
end


@inline function minimalCoordinates(joint::Rotational2, body1::Origin, body2::Body, No)
    q2 = joint.qoff \ body2.q[No]
    joint.V3 * axis(q2) * angle(q2) 
end

@inline function minimalCoordinates(joint::Rotational2, body1::Body, body2::Body, No)
    q2 = joint.qoff \ (body1.q[No] \ body2.q[No])
    joint.V3 * axis(q2) * angle(q2)
end


@inline function g(joint::Rotational2, body1::Origin, body2::Body, Δt, No)
    joint.V12 * VLᵀmat(joint.qoff) * getq3(body2, Δt)
end

@inline function g(joint::Rotational2, body1::Body, body2::Body, Δt, No)
    joint.V12 * (VLᵀmat(joint.qoff) * Lᵀmat(getq3(body1, Δt)) * getq3(body2, Δt))
end


@inline function ∂g∂posa(joint::Rotational2{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 2, 3)
        R = -joint.V12 * VLᵀmat(joint.qoff) * Rmat(body2.q[No]) * RᵀVᵀmat(body1.q[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational2{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 2, 3)
        R = joint.V12 * VLᵀmat(joint.qoff) * Lᵀmat(body1.q[No]) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational2{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 2, 3)
        Ω = joint.V12 * VLᵀmat(joint.qoff) * Rmat(ωbar(body2, Δt)) * Rmat(body2.q[No]) * Rᵀmat(body1.q[No]) * Tmat(T) * derivωbar(body1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational2{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 2, 3)
        Ω = joint.V12 * VLᵀmat(joint.qoff) * Lᵀmat(ωbar(body1, Δt)) * Lᵀmat(body1.q[No]) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function ∂g∂posb(joint::Rotational2{T}, body1::Origin, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 2, 3)
        R = joint.V12 * VLᵀmat(joint.qoff) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational2{T}, body1::Origin, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 2, 3)
        Ω = joint.V12 * VLᵀmat(joint.qoff) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
