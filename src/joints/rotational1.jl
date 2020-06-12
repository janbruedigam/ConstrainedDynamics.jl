# No idea what kind of joint this actually is...
@inline function getPositionDelta(joint::Rotational1, body1::AbstractBody, body2::Body, θ::SVector{2,T}) where T
    #TODO define this function
    @error("not defined for rot2")
end
@inline function getVelocityDelta(joint::Rotational1, body1::AbstractBody, body2::Body, ω::SVector{2,T}) where T
    #TODO define this function
    @error("not defined for rot2")
end

@inline function setForce!(joint::Rotational1, body1::Body, body2::Body, τ::SVector{2,T}) where T
    setForce!(joint, body1.state, body2.state, joint.V12' * τ)
    return
end
@inline function setForce!(joint::Rotational1, body1::Origin, body2::Body, τ::SVector{2,T}) where T
    setForce!(joint, body2.state, joint.V12' * τ)
    return
end

@inline function ∂Fτa∂u(joint::Rotational1, body1::Body)
    return ∂Fτa∂u(joint, body1.state) * joint.V12'
end
@inline function ∂Fτb∂u(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂Fτb∂u(joint, body1.state, body2.state) * joint.V12'
    else
        return ∂Fτb∂u(joint)
    end
end
@inline function ∂Fτb∂u(joint::Rotational1, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return return ∂Fτb∂u(joint, body2.state) * joint.V12'
    else
        return ∂Fτb∂u(joint)
    end
end

@inline function minimalCoordinates(joint::Rotational1, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    q = g(joint, statea.qc, stateb.qc)
    return joint.V12 * rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational1, body1::Origin, body2::Body)
    stateb = body2.state
    q = g(joint, stateb.qc)
    return joint.V12 * rotation_vector(q)
end

@inline g(joint::Rotational1, body1::Body, body2::Body, Δt) = joint.V3 * g(joint, body1.state, body2.state, Δt)
@inline g(joint::Rotational1, body1::Origin, body2::Body, Δt) = joint.V3 * g(joint, body2.state, Δt)

@inline function ∂g∂posa(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂posa(joint, body1.state, body2.state)
    else
        return ∂g∂posa(joint)
    end
end
@inline function ∂g∂posb(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂posb(joint, body1.state, body2.state)
    else
        return ∂g∂posb(joint)
    end
end
@inline function ∂g∂posb(joint::Rotational1, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂posb(joint, body2.state)
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂vela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂vela(joint)
    end
end
@inline function ∂g∂velb(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂velb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end
@inline function ∂g∂velb(joint::Rotational1, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂velb(joint, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end
