@inline function getPositionDelta(joint::Translational3, body1::AbstractBody, body2::Body, x::SVector{0,T}) where T
    Δx = szeros(T,3)
    return Δx
end
@inline function getVelocityDelta(joint::Translational3, body1::AbstractBody, body2::Body, v::SVector{0,T}) where T
    Δv = szeros(T,3)
    return Δv
end
@inline function getPositionDelta(joint::Rotational3, body1::AbstractBody, body2::Body, θ::SVector{0})
    Δq = joint.qoffset
    return Δq
end
@inline function getVelocityDelta(joint::Rotational3, body1::AbstractBody, body2::Body, ω::SVector{0,T}) where T
    Δω = szeros(T,3)
    return Δω
end

@inline minimalCoordinates(joint::Joint3{T}, body1::AbstractBody, body2::Body) where T = SA{T}[]

@inline constraintmat(::Joint3{T}) where T = SMatrix{3,3,T,9}(I)
@inline nullspacemat(::Joint3{T}) where T = szeros(T,0,3)