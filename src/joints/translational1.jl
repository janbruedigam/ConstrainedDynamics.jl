@inline function getPositionDelta(joint::Translational1, body1::AbstractBody, body2::Body, x::SVector{2,T}) where T
    Δx = joint.V12' * x # in body1 frame
    return Δx
end
@inline function getVelocityDelta(joint::Translational1, body1::AbstractBody, body2::Body, v::SVector{2,T}) where T
    Δv = joint.V12' * v # in body1 frame
    return Δv
end

@inline function setForce!(joint::Translational1, body1::Body, body2::Body, F::SVector{2,T}) where T
    setForce!(joint, body1.state, body2.state, joint.V12' * F)
    return
end
@inline function setForce!(joint::Translational1, body1::Origin, body2::Body, F::SVector{2,T}) where T
    setForce!(joint, body2.state, joint.V12' * F)
    return
end

@inline function ∂Fτa∂u(joint::Translational1, body1::Body)
    return ∂Fτa∂u(joint, body1.state) * joint.V12'
end
@inline function ∂Fτb∂u(joint::Translational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂Fτb∂u(joint, body1.state, body2.state) * joint.V12'
    else
        return ∂Fτb∂u(joint)
    end
end
@inline function ∂Fτb∂u(joint::Translational1, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return return ∂Fτb∂u(joint, body2.state) * joint.V12'
    else
        return ∂Fτb∂u(joint)
    end
end

@inline function minimalCoordinates(joint::Translational1, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return joint.V12 * g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
end
@inline function minimalCoordinates(joint::Translational1, body1::Origin, body2::Body)
    stateb = body2.state
    return joint.V12 * g(joint, stateb.xc, stateb.qc)
end

@inline g(joint::Translational1, body1::Body, body2::Body, Δt) = joint.V3 * g(joint, body1.state, body2.state, Δt)
@inline g(joint::Translational1, body1::Origin, body2::Body, Δt) = joint.V3 * g(joint, body2.state, Δt)

@inline function ∂g∂posa(joint::Translational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂posa(joint, body1.state, body2.state)
    else
        return ∂g∂posa(joint)
    end
end
@inline function ∂g∂posb(joint::Translational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂posb(joint, body1.state, body2.state)
    else
        return ∂g∂posb(joint)
    end
end
@inline function ∂g∂posb(joint::Translational1, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂posb(joint, body2.state)
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Translational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂vela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂vela(joint)
    end
end
@inline function ∂g∂velb(joint::Translational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂velb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end
@inline function ∂g∂velb(joint::Translational1, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂velb(joint, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end
