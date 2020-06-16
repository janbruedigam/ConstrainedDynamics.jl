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

@inline function ∂Fτ∂ua(joint::Translational1, body1::Body)
    return ∂Fτ∂ua(joint, body1.state) * joint.V12'
end
@inline function ∂Fτ∂ub(joint::Translational1, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂Fτ∂ub(joint, body1.state, body2.state) * joint.V12'
    else
        return ∂Fτ∂ub(joint)
    end
end
@inline function ∂Fτ∂ub(joint::Translational1, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return return ∂Fτ∂ub(joint, body2.state) * joint.V12'
    else
        return ∂Fτ∂ub(joint)
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

@inline function ∂g∂ᵣposa(joint::Translational1, body1::Body, body2::Body, args...)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣposa(joint, body1.state, body2.state, args...)
    else
        return ∂g∂ᵣposa(joint)
    end
end
@inline function ∂g∂ᵣposb(joint::Translational1, body1::Body, body2::Body, args...)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣposb(joint, body1.state, body2.state, args...)
    else
        return ∂g∂ᵣposb(joint)
    end
end
@inline function ∂g∂ᵣposb(joint::Translational1, body1::Origin, body2::Body, args...)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣposb(joint, body2.state, args...)
    else
        return ∂g∂ᵣposb(joint)
    end
end

@inline function ∂g∂ᵣvela(joint::Translational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ᵣvela(joint)
    end
end
@inline function ∂g∂ᵣvelb(joint::Translational1, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ᵣvelb(joint)
    end
end
@inline function ∂g∂ᵣvelb(joint::Translational1, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V3 * ∂g∂ᵣvelb(joint, body2.state, Δt)
    else
        return ∂g∂ᵣvelb(joint)
    end
end
