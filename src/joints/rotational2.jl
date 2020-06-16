@inline function getPositionDelta(joint::Rotational2, body1::AbstractBody, body2::Body, θ::SVector{1,T}) where T
    q = UnitQuaternion(cos(θ[1]/2), (joint.V3*sin(θ[1]/2))..., false)
    Δq = joint.qoff * q # in body1 frame
    return Δq
end
@inline function getVelocityDelta(joint::Rotational2, body1::Body, body2::Body, ω::SVector{1,T}) where T
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qc)*body1.state.qc*joint.qoff) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational2, body1::Origin, body2::Body, ω::SVector{1,T}) where T
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qc)*joint.qoff) # in body2 frame
    return Δω
end

@inline function setForce!(joint::Rotational2, body1::Body, body2::Body, τ::SVector{1,T}) where T
    setForce!(joint, body1.state, body2.state, joint.V3' * τ)
    return
end
@inline function setForce!(joint::Rotational2, body1::Origin, body2::Body, τ::SVector{1,T}) where T
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

@inline function ∂g∂ᵣposa(joint::Rotational2, body1::Body, body2::Body, args...)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ᵣposa(joint, body1.state, body2.state, args...)
    else
        return ∂g∂ᵣposa(joint)
    end
end
@inline function ∂g∂ᵣposb(joint::Rotational2, body1::Body, body2::Body, args...)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ᵣposb(joint, body1.state, body2.state, args...)
    else
        return ∂g∂ᵣposb(joint)
    end
end
@inline function ∂g∂ᵣposb(joint::Rotational2, body1::Origin, body2::Body, args...)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ᵣposb(joint, body2.state, args...)
    else
        return ∂g∂ᵣposb(joint)
    end
end

@inline function ∂g∂ᵣvela(joint::Rotational2, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ᵣvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ᵣvela(joint)
    end
end
@inline function ∂g∂ᵣvelb(joint::Rotational2, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ᵣvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ᵣvelb(joint)
    end
end
@inline function ∂g∂ᵣvelb(joint::Rotational2, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return joint.V12 * ∂g∂ᵣvelb(joint, body2.state, Δt)
    else
        return ∂g∂ᵣvelb(joint)
    end
end