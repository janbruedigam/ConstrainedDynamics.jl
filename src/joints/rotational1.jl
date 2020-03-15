# No idea what kind of joint this actually is...

mutable struct Rotational1{T,Nc} <: Joint{T,Nc}
    V12::SMatrix{2,3,T,6}
    V3::Adjoint{T,SVector{3,T}}
    qoff::Quaternion{T}
    cid::Int64

    function Rotational1(body1::AbstractBody{T}, body2::AbstractBody{T}, axis::AbstractVector{T};offset::Quaternion{T} = Quaternion{T}()) where T
        Nc = 1

        axis = vrotate(SVector(axis...),inv(offset))
        axis = axis / norm(axis)
        A = svd(skew(axis)).Vt # in frame of body1
        V12 = A[1:2,:]
        V3 = axis' # instead of A[3,:] for correct sign: abs(axis) = abs(A[3,:])
        cid = body2.id

        new{T,Nc}(V12, V3, offset, cid), body1.id, body2.id
    end
end

function setForce!(joint::Rotational1, body1::AbstractBody, body2::Body, τ::SVector{0,T}, No) where T
    return
end

function setForce!(joint::Rotational1, body1::AbstractBody, body2::Body, τ::Nothing, No)
    return
end

function setForce!(joint::Rotational1, body1::Body, body2::Body{T}, τ::SVector{2,T}, No) where T
    τ1 = vrotate(joint.V12' * -τ,body1.q[No])
    τ2 = -τ1

    body1.τ[No] = τ1
    body2.τ[No] = τ2
    return
end

@inline minimalCoordinates(joint::Rotational1, body1::Body, body2::Body, No) = joint.V12 * (VLᵀmat(body1.q[No]) * body2.q[No])

@inline g(joint::Rotational1, body1::Body, body2::Body, Δt, No) = joint.V3 * (VLᵀmat(joint.qoff) * Lᵀmat(getq3(body1, Δt)) * getq3(body2, Δt))

@inline function ∂g∂posa(joint::Rotational1{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 1, 3)

        R = -joint.V3 * VLᵀmat(joint.qoff) * Rmat(body2.q[No]) * RᵀVᵀmat(body1.q[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational1{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 1, 3)

        R = joint.V3 * VLᵀmat(joint.qoff) * Lᵀmat(body1.q[No]) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational1{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 1, 3)

        Ω = joint.V3 * VLᵀmat(joint.qoff) * Rmat(ωbar(body2, Δt)) * Rmat(body2.q[No]) * Rᵀmat(body1.q[No]) * Tmat(T) * derivωbar(body1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 1, 3)

        Ω = joint.V3 * VLᵀmat(joint.qoff) * Lᵀmat(ωbar(body1, Δt)) * Lᵀmat(body1.q[No]) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end

function setForce!(joint::Rotational1, body1::Origin, body2::Body{T}, τ::SVector{2,T}, No) where T
    body2.τ[No] = joint.V12' * τ
    return
end

@inline minimalCoordinates(joint::Rotational1, body1::Origin, body2::Body, No) = joint.V12 * Vmat(body2.q[No])

@inline g(joint::Rotational1, body1::Origin, body2::Body, Δt, No) = joint.V3 * VLᵀmat(joint.qoff) * getq3(body2, Δt)

@inline function ∂g∂posb(joint::Rotational1{T}, body1::Origin, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 1, 3)

        R = joint.V3 * VLᵀmat(joint.qoff) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational1{T}, body1::Origin, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 1, 3)

        Ω = joint.V3 * VLᵀmat(joint.qoff) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
