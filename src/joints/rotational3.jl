@inline function getPositionDelta(joint::Rotational3, body1::AbstractBody, body2::Body, θ::SVector{0,T}) where T
    Δq = joint.qoff
    return Δq
end
@inline function getVelocityDelta(joint::Rotational3, body1::AbstractBody, body2::Body, ω::SVector{0,T}) where T
    Δω = szeros(T,3)
    return Δω
end

@inline function ∂Fτ∂ua(joint::Rotational3{T}, body1::Body) where T
    return szeros(T, 6, 0)
end
@inline function ∂Fτ∂ub(joint::Rotational3{T}, body1::AbstractBody, body2::Body) where T
    return szeros(T, 6, 0)
end

@inline minimalCoordinates(joint::Rotational3{T}, body1::AbstractBody, body2::Body) where T = SA_F64[]

@inline g(joint::Rotational3, body1::Body, body2::Body, Δt) = g(joint, body1.state, body2.state, Δt)
@inline g(joint::Rotational3, body1::Origin, body2::Body, Δt) = g(joint, body2.state, Δt)

@inline function ∂g∂ʳposa(joint::Rotational3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂ʳposa(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposa(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Rotational3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂ʳposb(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Rotational3, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return ∂g∂ʳposb(joint, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end

@inline function ∂g∂ʳvela(joint::Rotational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ʳvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvela(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Rotational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ʳvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Rotational3, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ʳvelb(joint, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end