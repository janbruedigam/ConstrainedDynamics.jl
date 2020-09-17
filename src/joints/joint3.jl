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

@inline function ∂Fτ∂ua(joint::Joint3{T}, body1::Body) where T
    return szeros(T, 6, 0)
end
@inline function ∂Fτ∂ub(joint::Joint3{T}, body1::AbstractBody, body2::Body) where T
    return szeros(T, 6, 0)
end

@inline minimalCoordinates(joint::Joint3{T}, body1::AbstractBody, body2::Body) where T = SA{T}[]

@inline g(joint::Joint3, body1::Body, body2::Body, Δt) = g(joint, body1.state, body2.state, Δt)
@inline g(joint::Joint3, body1::Origin, body2::Body, Δt) = g(joint, body2.state, Δt)
@inline g(joint::Joint3, body1::Body, body2::Body) = g(joint, body1.state, body2.state)
@inline g(joint::Joint3, body1::Origin, body2::Body) = g(joint, body2.state)

@inline function ∂g∂ʳposa(joint::Joint3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂ʳposa(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposa(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Joint3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂ʳposb(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Joint3, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return ∂g∂ʳposb(joint, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end

@inline function ∂g∂posac(joint::Joint3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posac(joint, body1.state, body2.state)
    else
        return ∂g∂posac(joint)
    end
end
@inline function ∂g∂posbc(joint::Joint3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posbc(joint, body1.state, body2.state)
    else
        return ∂g∂posbc(joint)
    end
end
@inline function ∂g∂posbc(joint::Joint3, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posbc(joint, body2.state)
    else
        return ∂g∂posbc(joint)
    end
end

@inline function ∂g∂ʳvela(joint::Joint3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ʳvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvela(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Joint3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ʳvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Joint3, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂ʳvelb(joint, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end

@inline reductionmat(joint::Joint3{T}) where T = SMatrix{3,3,T,9}(I)