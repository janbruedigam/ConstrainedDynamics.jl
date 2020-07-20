@inline function getPositionDelta(joint::Translational2, body1::AbstractBody, body2::Body, x::SVector{1})
    Δx = joint.V3' * x # in body1 frame
    return Δx
end
@inline function getVelocityDelta(joint::Translational2, body1::AbstractBody, body2::Body, v::SVector{1})
    Δv = joint.V3' * v # in body1 frame
    return Δv
end

@inline function setForce!(joint::Translational2, body1::Body, body2::Body, F::SVector{1})
    setForce!(joint, body1.state, body2.state, joint.V3' * F)
    return
end
@inline function setForce!(joint::Translational2, body1::Origin, body2::Body, F::SVector{1})
    setForce!(joint, body2.state, joint.V3' * F)
    return
end

@inline function ∂Fτ∂ua(joint::Translational2, body1::Body)
    return ∂Fτ∂ua(joint, body1.state) * joint.V3'
end
@inline function ∂Fτ∂ub(joint::Translational2, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂Fτ∂ub(joint, body1.state, body2.state) * joint.V3'
    else
        return ∂Fτ∂ub(joint)
    end
end
@inline function ∂Fτ∂ub(joint::Translational2, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return return ∂Fτ∂ub(joint, body2.state) * joint.V3'
    else
        return ∂Fτ∂ub(joint)
    end
end

@inline function minimalCoordinates(joint::Translational2, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return joint.V3 * g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
end
@inline function minimalCoordinates(joint::Translational2, body1::Origin, body2::Body)
    stateb = body2.state
    return joint.V3 * g(joint, stateb.xc, stateb.qc)
end

@inline g(joint::Translational2, body1::Body, body2::Body, Δt) = joint.V12 * g(joint, body1.state, body2.state, Δt)
@inline g(joint::Translational2, body1::Origin, body2::Body, Δt) = joint.V12 * g(joint, body2.state, Δt)

@inline function ∂g∂ʳposa(joint::Translational2, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳposa(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposa(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Translational2, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳposb(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Translational2, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳposb(joint, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end

@inline function ∂g∂ʳvela(joint::Translational2, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvela(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Translational2, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Translational2, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳvelb(joint, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end

@inline reductionmat(joint::Translational2{T}) where T = joint.V12