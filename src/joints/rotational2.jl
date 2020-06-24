@inline function getPositionDelta(joint::Rotational2, body1::AbstractBody, body2::Body, θ::SVector{1})
    q = UnitQuaternion(cos(θ[1]/2), (joint.V3*sin(θ[1]/2))..., false)
    Δq = joint.qoff * q # in body1 frame
    return Δq
end
@inline function getVelocityDelta(joint::Rotational2, body1::Body, body2::Body, ω::SVector{1})
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qc)*body1.state.qc*joint.qoff) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational2, body1::Origin, body2::Body, ω::SVector{1})
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qc)*joint.qoff) # in body2 frame
    return Δω
end

@inline function setForce!(joint::Rotational2, body1::Body, body2::Body, τ::SVector{1})
    setForce!(joint, body1.state, body2.state, joint.V3' * τ)
    return
end
@inline function setForce!(joint::Rotational2, body1::Origin, body2::Body, τ::SVector{1})
    setForce!(joint, body2.state, joint.V3' * τ)
    return
end

@inline function ∂Fτ∂ua(joint::Rotational2, body1::Body)
    return ∂Fτ∂ua(joint, body1.state) * joint.V3'
end
@inline function ∂Fτ∂ub(joint::Rotational2, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂Fτ∂ub(joint, body1.state, body2.state) * joint.V3'
    else
        return ∂Fτ∂ub(joint)
    end
end
@inline function ∂Fτ∂ub(joint::Rotational2, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return return ∂Fτ∂ub(joint, body2.state) * joint.V3'
    else
        return ∂Fτ∂ub(joint)
    end
end

@inline function minimalCoordinates(joint::Rotational2, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    # q = g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
    q = joint.qoff \ (statea.qc \ stateb.qc)
    return joint.V3 * rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational2, body1::Origin, body2::Body)
    stateb = body2.state
    # q = g(joint, stateb.xc, stateb.qc)
    q = joint.qoff \ stateb.qc
    return joint.V3 * rotation_vector(q)
end

@inline g(joint::Rotational2, body1::Body, body2::Body, Δt) = joint.V12 * g(joint, body1.state, body2.state, Δt)
@inline g(joint::Rotational2, body1::Origin, body2::Body, Δt) = joint.V12 * g(joint, body2.state, Δt)

@inline function ∂g∂ʳposa(joint::Rotational2, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳposa(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposa(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Rotational2, body1::Body, body2::Body)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳposb(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Rotational2, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳposb(joint, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end

@inline function ∂g∂ʳvela(joint::Rotational2, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvela(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Rotational2, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Rotational2, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ʳvelb(joint, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end