# No idea what kind of joint this actually is...
@inline function getPositionDelta(joint::Rotational1, body1::AbstractBody, body2::Body, θ::SVector{2})
    #TODO define this function
    @error("not defined for rot2")
end
@inline function getVelocityDelta(joint::Rotational1, body1::AbstractBody, body2::Body, ω::SVector{2})
    #TODO define this function
    @error("not defined for rot2")
end

@inline function setForce!(joint::Rotational1, body1::Body, body2::Body, τ::SVector{2})
    setForce!(joint, body1.state, body2.state, joint.V12' * τ)
    return
end
@inline function setForce!(joint::Rotational1, body1::Origin, body2::Body, τ::SVector{2})
    setForce!(joint, body2.state, joint.V12' * τ)
    return
end

@inline function ∂Fτ∂ua(joint::Rotational1, body1::Body)
    return ∂Fτ∂ua(joint, body1.state) * joint.V12'
end
@inline function ∂Fτ∂ub(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂Fτ∂ub(joint, body1.state, body2.state) * joint.V12'
    else
        return ∂Fτ∂ub(joint)
    end
end
@inline function ∂Fτ∂ub(joint::Rotational1, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return return ∂Fτ∂ub(joint, body2.state) * joint.V12'
    else
        return ∂Fτ∂ub(joint)
    end
end

@inline function minimalCoordinates(joint::Rotational1, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    # q = g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
    q = statea.qc \ stateb.qc / joint.qoffset
    return joint.V12 * rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational1, body1::Origin, body2::Body)
    stateb = body2.state
    # q = g(joint, stateb.xc, stateb.qc)
    q = stateb.qc / joint.qoffset
    return joint.V12 * rotation_vector(q)
end

@inline g(joint::Rotational1, body1::Body, body2::Body, Δt) = joint.V3 * g(joint, body1.state, body2.state, Δt)
@inline g(joint::Rotational1, body1::Origin, body2::Body, Δt) = joint.V3 * g(joint, body2.state, Δt)

@inline function ∂g∂ʳposa(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ʳposa(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposa(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ʳposb(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Rotational1, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ʳposb(joint, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end

@inline function ∂g∂ʳvela(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ʳvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvela(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ʳvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Rotational1, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ʳvelb(joint, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end

@inline reductionmat(joint::Rotational1{T}) where T = joint.V3