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
    q = joint.qoff \ (statea.qc \ stateb.qc)
    return joint.V12 * rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational1, body1::Origin, body2::Body)
    stateb = body2.state
    # q = g(joint, stateb.xc, stateb.qc)
    q = joint.qoff \ stateb.qc
    return joint.V12 * rotation_vector(q)
end

@inline g(joint::Rotational1, body1::Body, body2::Body, Δt) = joint.V3 * g(joint, body1.state, body2.state, Δt)
@inline g(joint::Rotational1, body1::Origin, body2::Body, Δt) = joint.V3 * g(joint, body2.state, Δt)

@inline function ∂g∂ᵣposa(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣposa(joint, body1.state, body2.state)
    else
        return ∂g∂ᵣposa(joint)
    end
end
@inline function ∂g∂ᵣposb(joint::Rotational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣposb(joint, body1.state, body2.state)
    else
        return ∂g∂ᵣposb(joint)
    end
end
@inline function ∂g∂ᵣposb(joint::Rotational1, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣposb(joint, body2.state)
    else
        return ∂g∂ᵣposb(joint)
    end
end

@inline function ∂g∂ᵣvela(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ᵣvela(joint)
    end
end
@inline function ∂g∂ᵣvelb(joint::Rotational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ᵣvelb(joint)
    end
end
@inline function ∂g∂ᵣvelb(joint::Rotational1, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣvelb(joint, body2.state, Δt)
    else
        return ∂g∂ᵣvelb(joint)
    end
end
