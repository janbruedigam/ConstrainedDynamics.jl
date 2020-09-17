@inline function getPositionDelta(joint::Translational1, body1::AbstractBody, body2::Body, x::SVector{2})
    Δx = joint.V12' * x # in body1 frame
    return Δx
end
@inline function getVelocityDelta(joint::Translational1, body1::AbstractBody, body2::Body, v::SVector{2})
    Δv = joint.V12' * v # in body1 frame
    return Δv
end
# No idea what kind of joint this actually is...
@inline function getPositionDelta(joint::Rotational1, body1::AbstractBody, body2::Body, θ::SVector{2})
    #TODO define this function
    @error("not defined for rot2")
end
@inline function getVelocityDelta(joint::Rotational1, body1::AbstractBody, body2::Body, ω::SVector{2})
    #TODO define this function
    @error("not defined for rot2")
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

@inline constraintmat(joint::Joint1) = joint.V3
@inline nullspacemat(joint::Joint1) = joint.V12