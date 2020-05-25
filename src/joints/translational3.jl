@inline function getPositionDelta(joint::Translational3, body1::AbstractBody, body2::Body{T}, x::SVector{0,T}) where T
    Δx = @SVector zeros(T,3)
    return Δx
end

@inline function getVelocityDelta(joint::Translational3, body1::AbstractBody, body2::Body{T}, v::SVector{0,T}) where T
    Δv = @SVector zeros(T,3)
    return Δv
end

@inline function minimalCoordinates(joint::Translational3, body1::AbstractBody{T}, body2::Body, No) where T
    SVector{0,T}()
end


@inline function g(joint::Translational3, body1::Body, body2::Body, Δt)
    vertices = joint.vertices
    q1 = getq2(body1, Δt)
    # getx2(body2, Δt) + vrotate(vertices[2], getq2(body2, Δt)) - (getx2(body1, Δt) + vrotate(vertices[1], getq2(body1, Δt)))
    vrotate(getx2(body2, Δt) + vrotate(vertices[2], getq2(body2, Δt)) - (getx2(body1, Δt) + vrotate(vertices[1], q1)), inv(q1))
end

@inline function g(joint::Translational3, body1::Origin, body2::Body, Δt)
    vertices = joint.vertices
    # getx2(body2, Δt) + vrotate(vertices[2], getq2(body2, Δt)) - vertices[1]
    (getx2(body2, Δt) + vrotate(vertices[2], getq2(body2, Δt)) - vertices[1])
end


@inline function ∂g∂posa(joint::Translational3{T}, body1::Body, body2::Body) where T
    if body2.id == joint.cid
        q1 = body1.state.qd[2]
        point2 = body2.state.xd[2] + vrotate(joint.vertices[2], body2.state.qd[2])

        # X = SMatrix{3,3,T,9}(-I)
        # R = -2 * VRᵀmat(q1) * Rmat(Quaternion(joint.vertices[1])) * LVᵀmat(q1)
        X = -VLᵀmat(q1) * RVᵀmat(q1)
        R = 2 * VLᵀmat(q1) * (Lmat(Quaternion(point2)) - Lmat(Quaternion(body1.state.xd[2]))) * LVᵀmat(q1)

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Translational3{T}, body1::Body, body2::Body) where T
    if body2.id == joint.cid
        q1 = body1.state.qd[2]
        q2 = body2.state.qd[2]

        # X = SMatrix{3,3,T,9}(I)
        # R = 2 * VRᵀmat(q2) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(q2)
        X = VLᵀmat(q1)RVᵀmat(q1)
        R = 2 * VLᵀmat(q1) * Rmat(q1) * Rᵀmat(q2) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(q2)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Translational3{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q1 = body1.state.qd[No]
        ω1 = getω2(body1)
        ω2 = getω2(body2)
        ωbar1 = ωbar(ω1, Δt)
        point2 = body2.state.xd[No] + Δt * getv2(body2) + vrotate(vrotate(joint.vertices[2], ωbar(ω2, Δt)), body2.state.qd[No])

        # V = SMatrix{3,3,T,9}(-Δt * I)        
        # Ω = -2 * VRᵀmat(q1) * Lmat(q1) * Rᵀmat(ωbar(ω1, Δt)) * Rmat(Quaternion(joint.vertices[1])) * derivωbar(ω1, Δt)
        V = -Δt * VLᵀmat(ωbar1) * Lᵀmat(q1) * Rmat(ωbar1) * RVᵀmat(q1)
        Ω = 2 * VLᵀmat(ωbar1) * Lᵀmat(q1) * (Lmat(Quaternion(point2)) - Lmat(Quaternion(body1.state.xd[No] + Δt * getv2(body1)))) * Lmat(q1) * derivωbar(ω1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Translational3{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q1 = body1.state.qd[No]
        q2 = body2.state.qd[No]
        ω1 = getω2(body1)
        ω2 = getω2(body2)
        ωbar1 = ωbar(ω1, Δt)

        # V = SMatrix{3,3,T,9}(Δt * I)
        # Ω = 2 * VRᵀmat(q2) * Lmat(q2) * Rᵀmat(ωbar(ω2, Δt)) * Rmat(Quaternion(joint.vertices[2])) * derivωbar(ω2, Δt)
        V = Δt * VLᵀmat(ωbar1)Lᵀmat(q1)Rmat(ωbar1)RVᵀmat(q1)
        Ω = 2 * VLᵀmat(ωbar1) * Lᵀmat(q1) * Lmat(q2) * Rmat(ωbar1) * Rmat(q1) * Rᵀmat(q2) * Rᵀmat(ωbar(ω2, Δt)) * Rmat(Quaternion(joint.vertices[2])) * derivωbar(ω2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function ∂g∂posb(joint::Translational3{T}, body1::Origin, body2::Body) where T
    if body2.id == joint.cid
        q2 = body2.state.qd[2]

        X = SMatrix{3,3,T,9}(I)
        R = 2 * VRᵀmat(q2) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(q2)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Translational3{T}, body1::Origin, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q2 = body2.state.qd[No]
        ω2 = getω2(body2)

        V = SMatrix{3,3,T,9}(Δt * I)
        Ω = 2 * VLmat(q2) * Rᵀmat(q2) * Rᵀmat(ωbar(ω2, Δt)) * Rmat(Quaternion(joint.vertices[2])) * derivωbar(ω2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end