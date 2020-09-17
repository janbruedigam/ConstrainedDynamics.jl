@inline function getPositionDelta(joint::Translational2, body1::AbstractBody, body2::Body, x::SVector{1})
    Δx = joint.V3' * x # in body1 frame
    return Δx
end
@inline function getVelocityDelta(joint::Translational2, body1::AbstractBody, body2::Body, v::SVector{1})
    Δv = joint.V3' * v # in body1 frame
    return Δv
end
@inline function getPositionDelta(joint::Rotational2, body1::AbstractBody, body2::Body, θ::SVector{1})
    q = UnitQuaternion(cos(θ[1]/2), (joint.V3*sin(θ[1]/2))..., false)
    Δq = q * joint.qoffset # in body1 frame
    return Δq
end
@inline function getVelocityDelta(joint::Rotational2, body1::Body, body2::Body, ω::SVector{1})
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qc)*body1.state.qc) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational2, body1::Origin, body2::Body, ω::SVector{1})
    ω = joint.V3' * ω
    Δω = vrotate(ω, inv(body2.state.qc)) # in body2 frame
    return Δω
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
@inline function minimalCoordinates(joint::Rotational2, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    # q = g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
    q = statea.qc \ stateb.qc / joint.qoffset
    return joint.V3 * rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational2, body1::Origin, body2::Body)
    stateb = body2.state
    # q = g(joint, stateb.xc, stateb.qc)
    q = stateb.qc / joint.qoffset
    return joint.V3 * rotation_vector(q)
end

@inline constraintmat(joint::Joint2) = joint.V12
@inline nullspacemat(joint::Joint2) = joint.V3