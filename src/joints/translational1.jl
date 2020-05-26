@inline function getPositionDelta(joint::Translational1, body1::AbstractBody, body2::Body{T}, x::SVector{2,T}) where T
    Δx = joint.V12' * x # in body1 frame
    return Δx
end

@inline function getVelocityDelta(joint::Translational1, body1::AbstractBody, body2::Body{T}, v::SVector{2,T}) where T
    Δv = joint.V12' * v # in body1 frame
    return Δv
end

@inline function setForce!(joint::Translational1, body1::Body, body2::Body{T}, F::SVector{2,T}, No) where T
    clearForce!(joint, body1, body2, No)

    q1 = body1.state.qd[No]
    q2 = body2.state.qd[No]

    F1 = vrotate(joint.V12' * -F, q1)
    F2 = -F1

    τ1 = vrotate(torqueFromForce(F1, vrotate(joint.vertices[1], q1)),inv(q1)) # in local coordinates
    τ2 = vrotate(torqueFromForce(F2, vrotate(joint.vertices[2], q2)),inv(q2)) # in local coordinates

    updateForce!(joint, body1, body2, F1, τ1, F2, τ2, No)
    return
end

@inline function setForce!(joint::Translational1, body1::Origin, body2::Body{T}, F::SVector{2,T}, No) where T
    clearForce!(joint, body2, No)

    q2 = body2.state.qd[No]

    F2 = joint.V12' * F
    τ2 = vrotate(torqueFromForce(F2, vrotate(joint.vertices[2], q2)),inv(q2)) # in local coordinates

    updateForce!(joint, body2, F2, τ2, No)
    return
end


@inline function minimalCoordinates(joint::Translational1, body1::Body, body2::Body, No)
    vertices = joint.vertices
    q1 = body1.state.qd[No]
    joint.V12 * vrotate(body2.state.xd[No] + vrotate(vertices[2], body2.state.qd[No]) - (body1.state.xd[No] + vrotate(vertices[1], q1)), inv(q1))
end

@inline function minimalCoordinates(joint::Translational1, body1::Origin, body2::Body, No)
    vertices = joint.vertices
    joint.V12 * (body2.state.xd[No] + vrotate(vertices[2], body2.state.qd[No]) - vertices[1])
end


@inline function g(joint::Translational1, body1::Body, body2::Body, Δt)
    joint.V3 * g(joint, getx2(body1, Δt), getq2(body1, Δt), getx2(body2, Δt), getq2(body2, Δt))
end

@inline function g(joint::Translational1, body1::Origin, body2::Body, Δt)
    joint.V3 * g(joint, getx2(body2, Δt), getq2(body2, Δt))
end


@inline function ∂g∂posa(joint::Translational1, body1::Body, body2::Body)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂posa(joint, getxd2(body1), getqd2(body1), getxd2(body2), getqd2(body2))
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Translational1, body1::Body, body2::Body)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂posb(joint, getxd2(body1), getqd2(body1), getxd2(body2), getqd2(body2))
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂posb(joint::Translational1, body1::Origin, body2::Body)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂posb(joint, getxd2(body2), getqd2(body2))
    else
        return ∂g∂posb(joint)
    end
end


@inline function ∂g∂vela(joint::Translational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂vela(joint, getx1(body1), getx2(body1, Δt), getq1(body1), getq2(body1, Δt), getvupdate(body1), getωupdate(body1), getx2(body2, Δt), getq2(body2, Δt), Δt)
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Translational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂velb(joint, getx2(body1, Δt), getq2(body1, Δt), getx1(body2), getx2(body2, Δt), getq1(body2), getq2(body2, Δt), getvupdate(body2), getωupdate(body2), Δt)
    else
        return ∂g∂velb(joint)
    end
end

@inline function ∂g∂velb(joint::Translational1, body1::Origin, body2::Body, Δt)
    if body2.id == joint.cid
        return joint.V3 * ∂g∂velb(joint, getx1(body2), getx2(body2, Δt), getq1(body2), getq2(body2, Δt), getvupdate(body2), getωupdate(body2), Δt)
    else
        return ∂g∂velb(joint)
    end
end
