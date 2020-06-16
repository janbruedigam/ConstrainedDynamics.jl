@inline function getPositionDelta(joint::Translational3, body1::AbstractBody, body2::Body, x::SVector{0,T}) where T
    Δx = szeros(T,3)
    return Δx
end
@inline function getVelocityDelta(joint::Translational3, body1::AbstractBody, body2::Body, v::SVector{0,T}) where T
    Δv = szeros(T,3)
    return Δv
end

@inline function ∂Fτ∂ua(joint::Translational3{T}, body1::Body) where T
    return szeros(T, 6, 0)
end
@inline function ∂Fτ∂ub(joint::Translational3{T}, body1::AbstractBody, body2::Body) where T
    return szeros(T, 6, 0)
end

@inline minimalCoordinates(joint::Translational3{T}, body1::AbstractBody, body2::Body) where T = SA_F64[]

@inline g(joint::Translational3, body1::Body, body2::Body, Δt) = g(joint, body1.state, body2.state, Δt)
@inline g(joint::Translational3, body1::Origin, body2::Body, Δt) = g(joint, body2.state, Δt)

@inline function ∂g∂ᵣposa(joint::Translational3, body1::Body, body2::Body, args...)
    if body2.id == joint.childid
        return ∂g∂ᵣposa(joint, body1.state, body2.state, args...)
    else
        return ∂g∂ᵣposa(joint)
    end
end
@inline function ∂g∂ᵣposb(joint::Translational3, body1::Body, body2::Body, args...)
    if body2.id == joint.childid
        return ∂g∂ᵣposb(joint, body1.state, body2.state, args...)
    else
        return ∂g∂ᵣposb(joint)
    end
end
@inline function ∂g∂ᵣposb(joint::Translational3, body1::Origin, body2::Body, args...)
    if body2.id == joint.childid
        return ∂g∂ᵣposb(joint, body2.state, args...)
    else
        return ∂g∂ᵣposb(joint)
    end
end

@inline function ∂g∂ᵣvela(joint::Translational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ᵣvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ᵣvela(joint)
    end
end
@inline function ∂g∂ᵣvelb(joint::Translational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ᵣvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ᵣvelb(joint)
    end
end
@inline function ∂g∂ᵣvelb(joint::Translational3, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ᵣvelb(joint, body2.state, Δt)
    else
        return ∂g∂ᵣvelb(joint)
    end
end