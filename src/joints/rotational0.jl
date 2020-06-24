@inline function getPositionDelta(joint::Rotational0, body1::AbstractBody, body2::Body, θ::SVector{3,T}) where T
    # axis angle representation
    if norm(θ) == 0
        q = one(UnitQuaternion{T})
    else
        q = UnitQuaternion(cos(norm(θ)/2),(θ/norm(θ)*sin(norm(θ)/2))..., false)
    end
    
    Δq = joint.qoff * q # in body1 frame
    return Δq
end
@inline function getVelocityDelta(joint::Rotational0, body1::Body, body2::Body, ω::SVector{3})
    Δω = vrotate(ω, inv(body2.state.qc)*body1.state.qc*joint.qoff) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational0, body1::Origin, body2::Body, ω::SVector{3})
    Δω = vrotate(ω, inv(body2.state.qc)*joint.qoff) # in body2 frame
    return Δω
end

@inline function setForce!(joint::Rotational0, body1::Body, body2::Body, τ::SVector{3})
    setForce!(joint, body1.state, body2.state, τ)
    return
end
@inline function setForce!(joint::Rotational0, body1::Origin, body2::Body, τ::SVector{3})
    setForce!(joint, body2.state, τ)
    return
end

@inline function ∂Fτ∂ua(joint::Rotational0, body1::Body)
    return ∂Fτ∂ua(joint, body1.state)
end
@inline function ∂Fτ∂ub(joint::Rotational0, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂Fτ∂ub(joint, body1.state, body2.state)
    else
        return ∂Fτ∂ub(joint)
    end
end
@inline function ∂Fτ∂ub(joint::Rotational0, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return return ∂Fτ∂ub(joint, body2.state)
    else
        return ∂Fτ∂ub(joint)
    end
end

@inline function minimalCoordinates(joint::Rotational0, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    # q = g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
    q = joint.qoff \ (statea.qc \ stateb.qc)
    return rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational0, body1::Origin, body2::Body)
    stateb = body2.state
    # q = g(joint, stateb.xc, stateb.qc)
    q = joint.qoff \ stateb.qc
    return rotation_vector(q)
end

@inline g(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt) = g(joint)

@inline ∂g∂ʳposa(joint::Rotational0, body1::AbstractBody, body2::AbstractBody) = ∂g∂ʳposa(joint)
@inline ∂g∂ʳposb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody) = ∂g∂ʳposb(joint)
@inline ∂g∂ʳvela(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt) = ∂g∂ʳvela(joint)
@inline ∂g∂ʳvelb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt) = ∂g∂ʳvelb(joint)
