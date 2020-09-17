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

@inline constraintmat(::Joint0{T}) where T = szeros(T,0,3)
@inline nullspacemat(::Joint0{T}) where T = SMatrix{3,3,T,9}(I)