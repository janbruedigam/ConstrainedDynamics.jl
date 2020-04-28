@inline function getVelocityDelta(joint::Translational2, body1::AbstractBody, body2::Body{T}, v::Union{T,SVector{1,T}}) where T
    Δv = joint.V3' * v # in body1 frame
    return Δv
end

@inline function getPositionDelta(joint::Translational2, body1::AbstractBody, body2::Body{T}, x::Union{T,SVector{1,T}}) where T
    Δx = joint.V3' * x # in body1 frame
    return Δx
end

@inline function setForce!(joint::Translational2, body1::Body, body2::Body{T}, F::Union{T,SVector{1,T}}, No) where T
    q1 = body1.q[No]
    q2 = body2.q[No]

    F1 = vrotate(joint.V3' * -F, q1)
    F2 = -F1

    τ1 = torqueFromForce(F1, vrotate(joint.vertices[1], q1))
    τ2 = torqueFromForce(F2, vrotate(joint.vertices[2], q2))

    body1.F[No] = F1
    body1.τ[No] = τ1

    body2.F[No] = F2
    body2.τ[No] = τ2
    return
end

@inline function setForce!(joint::Translational2, body1::Origin, body2::Body{T}, F::Union{T,SVector{1,T}}, No) where T
    F2 = joint.V3' * F
    τ2 = torqueFromForce(F2, vrotate(joint.vertices[2], body2.q[No]))

    body2.F[No] = F2
    body2.τ[No] = τ2
    return
end


@inline function minimalCoordinates(joint::Translational2, body1::Body, body2::Body, No)
    vertices = joint.vertices
    q1 = body1.q[No]
    joint.V3 * vrotate(body2.x[No] + vrotate(vertices[2], body2.q[No]) - (body1.x[No] + vrotate(vertices[1], q1)), inv(q1))
end

@inline function minimalCoordinates(joint::Translational2, body1::Origin, body2::Body, No)
    vertices = joint.vertices
    joint.V3 * (body2.x[No] + vrotate(vertices[2], body2.q[No]) - vertices[1])
end


@inline function g(joint::Translational2, body1::Body, body2::Body, Δt, No)
    vertices = joint.vertices
    q1 = getq3(body1, Δt)
    joint.V12 * vrotate(getx3(body2, Δt) + vrotate(vertices[2], getq3(body2, Δt)) - (getx3(body1, Δt) + vrotate(vertices[1], q1)), inv(q1))
end

@inline function g(joint::Translational2, body1::Origin, body2::Body, Δt, No)
    vertices = joint.vertices
    joint.V12 * (getx3(body2, Δt) + vrotate(vertices[2], getq3(body2, Δt)) - vertices[1])
end


@inline function ∂g∂posa(joint::Translational2{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        point2 = body2.x[No] + vrotate(joint.vertices[2], body2.q[No])

        X = -joint.V12 * VLᵀmat(q1) * RVᵀmat(q1)
        R = joint.V12 * 2 * VLᵀmat(q1) * (Lmat(Quaternion(point2)) - Lmat(Quaternion(body1.x[No]))) * LVᵀmat(q1)

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Translational2{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        q2 = body2.q[No]

        X = joint.V12 * VLᵀmat(q1)RVᵀmat(q1)
        R = joint.V12 * 2 * VLᵀmat(q1) * Rmat(q1) * Rᵀmat(q2) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(q2)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Translational2{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        ωbar1 = ωbar(body1, Δt)
        point2 = body2.x[No] + Δt * getvnew(body2) + vrotate(vrotate(joint.vertices[2], ωbar(body2, Δt)), body2.q[No])

        V = -Δt * joint.V12 * VLᵀmat(ωbar1) * Lᵀmat(q1) * Rmat(ωbar1) * RVᵀmat(q1)
        Ω = 2 * joint.V12 * VLᵀmat(ωbar1) * Lᵀmat(q1) * (Lmat(Quaternion(point2)) - Lmat(Quaternion(body1.x[No] + Δt * getvnew(body1)))) * Lmat(q1) * derivωbar(body1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Translational2{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        q2 = body2.q[No]
        ωbar1 = ωbar(body1, Δt)

        V = Δt * joint.V12 * VLᵀmat(ωbar1)Lᵀmat(q1)Rmat(ωbar1)RVᵀmat(q1)
        Ω = 2 * joint.V12 * VLᵀmat(ωbar1) * Lᵀmat(q1) * Lmat(q2) * Rmat(ωbar1) * Rmat(q1) * Rᵀmat(q2) * Rᵀmat(ωbar(body2, Δt)) * Rmat(Quaternion(joint.vertices[2])) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function ∂g∂posb(joint::Translational2{T}, body1::Origin, body2::Body, No) where T
    if body2.id == joint.cid
        q2 = body2.q[No]

        X = joint.V12
        R = joint.V12 * 2 * VRᵀmat(q2) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(q2)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Translational2{T}, body1::Origin, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q2 = body2.q[No]

        V = Δt * joint.V12
        Ω = 2 * joint.V12 * VLmat(q2) * Rᵀmat(q2) * Rᵀmat(ωbar(body2, Δt)) * Rmat(Quaternion(joint.vertices[2])) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
