@inline function getPositionDelta(joint::Translational0, body1::AbstractBody, body2::Body, x::SVector{3})
    Δx = x # in body1 frame
    return Δx
end
@inline function getVelocityDelta(joint::Translational0, body1::AbstractBody, body2::Body, v::SVector{3})
    Δv = v # in body1 frame
    return Δv
end
@inline function getPositionDelta(joint::Rotational0, body1::AbstractBody, body2::Body, θ::SVector{3,T}) where T
    # axis angle representation
    nθ = norm(θ)
    if nθ == 0
        q = one(UnitQuaternion{T})
    else
        q = UnitQuaternion(cos(nθ/2),(θ/nθ*sin(nθ/2))..., false)
    end
    
    Δq = q * joint.qoffset # in body1 frame
    return Δq
end
@inline function getVelocityDelta(joint::Rotational0, body1::Body, body2::Body, ω::SVector{3})
    Δω = vrotate(ω, inv(body2.state.qc)*body1.state.qc) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational0, body1::Origin, body2::Body, ω::SVector{3})
    Δω = vrotate(ω, inv(body2.state.qc)) # in body2 frame
    return Δω
end

@inline function setForce!(joint::Joint0, body1::Body, body2::Body, Fτ::SVector{3})
    setForce!(joint, body1.state, body2.state, Fτ)
    return
end
@inline function setForce!(joint::Joint0, body1::Origin, body2::Body, Fτ::SVector{3})
    setForce!(joint, body2.state, Fτ)
    return
end

@inline function ∂Fτ∂ua(joint::Joint0, body1::Body)
    return ∂Fτ∂ua(joint, body1.state)
end
@inline function ∂Fτ∂ub(joint::Joint0, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂Fτ∂ub(joint, body1.state, body2.state)
    else
        return ∂Fτ∂ub(joint)
    end
end
@inline function ∂Fτ∂ub(joint::Joint0, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return return ∂Fτ∂ub(joint, body2.state)
    else
        return ∂Fτ∂ub(joint)
    end
end

@inline function minimalCoordinates(joint::Translational0, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
end
@inline function minimalCoordinates(joint::Translational0, body1::Origin, body2::Body)
    stateb = body2.state
    return g(joint, stateb.xc, stateb.qc)
end
@inline function minimalCoordinates(joint::Rotational0, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    # q = g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
    q = statea.qc \ stateb.qc / joint.qoffset
    return rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational0, body1::Origin, body2::Body)
    stateb = body2.state
    # q = g(joint, stateb.xc, stateb.qc)
    q = stateb.qc / joint.qoffset
    return rotation_vector(q)
end

@inline g(joint::Joint0, body1::AbstractBody, body2::AbstractBody, Δt) = g(joint)

@inline ∂g∂ʳposa(joint::Joint0, body1::AbstractBody, body2::AbstractBody) = ∂g∂ʳposa(joint)
@inline ∂g∂ʳposb(joint::Joint0, body1::AbstractBody, body2::AbstractBody) = ∂g∂ʳposb(joint)
@inline ∂g∂ʳvela(joint::Joint0, body1::AbstractBody, body2::AbstractBody, Δt) = ∂g∂ʳvela(joint)
@inline ∂g∂ʳvelb(joint::Joint0, body1::AbstractBody, body2::AbstractBody, Δt) = ∂g∂ʳvelb(joint)

@inline reductionmat(joint::Joint0{T}) where T = szeros(T,0,3)