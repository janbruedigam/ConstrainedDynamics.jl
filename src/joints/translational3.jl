@inline function getVelocityDelta(joint::Translational3, body1::AbstractBody, body2::Body{T}, v::Union{Nothing,SVector{0,T}}) where T
    Δv = @SVector zeros(T,3)
    return Δv
end

@inline function getPositionDelta(joint::Translational3, body1::AbstractBody, body2::Body{T}, x::Union{Nothing,SVector{0,T}}) where T
    Δx = @SVector zeros(T,3)
    return Δx
end

@inline function minimalCoordinates(joint::Translational3, body1::AbstractBody{T}, body2::Body, No) where T
    SVector{0,T}()
end


@inline function g(joint::Translational3, body1::Body, body2::Body, Δt, No)
    vertices = joint.vertices
    getx3(body2, Δt) + vrotate(vertices[2], getq3(body2, Δt)) - (getx3(body1, Δt) + vrotate(vertices[1], getq3(body1, Δt)))
end

@inline function g(joint::Translational3, body1::Origin, body2::Body, Δt, No)
    vertices = joint.vertices
    getx3(body2, Δt) + vrotate(vertices[2], getq3(body2, Δt)) - vertices[1]
end


@inline function ∂g∂posa(joint::Translational3{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]

        X = SMatrix{3,3,T,9}(-I)
        R = -2 * VRᵀmat(q1) * Rmat(Quaternion(joint.vertices[1])) * LVᵀmat(q1)

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Translational3{T}, body1::AbstractBody, body2::Body, No) where T
    if body2.id == joint.cid
        q2 = body2.q[No]

        X = SMatrix{3,3,T,9}(I)
        R = 2 * VRᵀmat(q2) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(q2)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Translational3{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]

        V = SMatrix{3,3,T,9}(-Δt * I)        
        Ω = -2 * VRᵀmat(q1) * Lmat(q1) * Rᵀmat(ωbar(body1, Δt)) * Rmat(Quaternion(joint.vertices[1])) * derivωbar(body1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Translational3{T}, body1::AbstractBody, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q2 = body2.q[No]

        V = SMatrix{3,3,T,9}(Δt * I)
        Ω = 2 * VRᵀmat(q2) * Lmat(q2) * Rᵀmat(ωbar(body2, Δt)) * Rmat(Quaternion(joint.vertices[2])) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
