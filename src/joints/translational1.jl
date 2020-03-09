mutable struct Translational1{T,Nc} <: Joint{T,Nc}
    vertices::NTuple{2,SVector{3,T}}
    V12::SMatrix{2,3,T,6}
    V3::Adjoint{T,SVector{3,T}}
    cid::Int64

    function Translational1(body1::AbstractBody{T}, body2::AbstractBody{T}, p1::AbstractVector{T}, p2::AbstractVector{T}, axis::AbstractVector{T}) where T
        Nc = 1
        vertices = (p1, p2)

        axis = axis / norm(axis)
        A = svd(skew(axis)).Vt
        V12 = A[1:2,:]
        V3 = axis' # instead of A[3,:] for correct sign: abs(axis) = abs(A[3,:])
        cid = body2.id

        new{T,Nc}(vertices, V12, V3, cid), body1.id, body2.id
    end
end

function setForce!(joint::Translational1, body1::AbstractBody, body2::Body, F::SVector{0,T}, No) where T
    return
end

function setForce!(joint::Translational1, body1::AbstractBody, body2::Body, F::Nothing, No)
    return
end

function setForce!(joint::Translational1, body1::Body, body2::Body{T}, F::SVector{2,T}, No) where T
    F1 = joint.V12' * -F
    F2 = -F1
    τ1 = torqueFromForce(F1, joint.vertices[1])
    τ2 = torqueFromForce(F2, joint.vertices[2])

    body1.F[No] = F1
    body1.τ[No] = τ1

    body2.F[No] = F2
    body2.τ[No] = τ2
    return
end

function minimalCoordinates(joint::Translational1, body1::Body, body2::Body, No)
    vertices = joint.vertices
    q1 = body1.q[No]
    joint.V12 * vrotate(body2.x[No] + vrotate(vertices[2], body2.q[No]) - (body1.x[No] + vrotate(vertices[1], q1)), inv(q1))
end

@inline function g(joint::Translational1, body1::Body, body2::Body, Δt, No)
    vertices = joint.vertices
    q1 = getq3(body1, Δt)
    joint.V3 * vrotate(getx3(body2, Δt) + vrotate(vertices[2], getq3(body2, Δt)) - (getx3(body1, Δt) + vrotate(vertices[1], q1)), inv(q1))
end

@inline function ∂g∂posa(joint::Translational1{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        point2 = body2.x[No] + vrotate(joint.vertices[2], body2.q[No])

        X = -joint.V3 * VLᵀmat(q1) * RVᵀmat(q1)

        R = joint.V3 * 2 * VLᵀmat(q1) * (Lmat(Quaternion(point2)) - Lmat(Quaternion(body1.x[No]))) * LVᵀmat(q1)

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Translational1{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        q2 = body2.q[No]

        X = joint.V3 * VLᵀmat(q1)RVᵀmat(q1)

        R = joint.V3 * 2 * VLᵀmat(q1) * Rmat(q1) * Rᵀmat(q2) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(q2)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Translational1{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        ωbar1 = ωbar(body1, Δt)
        point2 = body2.x[No] + Δt * getvnew(body2) + vrotate(vrotate(joint.vertices[2], ωbar(body2, Δt)), body2.q[No])

        V = -Δt * joint.V3 * VLᵀmat(ωbar1)Lᵀmat(q1)Rmat(ωbar1)RVᵀmat(q1)

        Ω = 2 * joint.V3 * VLᵀmat(ωbar1) * Lᵀmat(q1) * (Lmat(Quaternion(point2)) - Lmat(Quaternion(body1.x[No] + Δt * getvnew(body1)))) * Lmat(q1) * derivωbar(body1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Translational1{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        q2 = body2.q[No]
        ωbar1 = ωbar(body1, Δt)
        V = Δt * joint.V3 * VLᵀmat(ωbar1)Lᵀmat(q1)Rmat(ωbar1)RVᵀmat(q1)

        Ω = 2 * joint.V3 * VLᵀmat(ωbar1) * Lᵀmat(q1) * Lmat(q2) * Rmat(ωbar1) * Rmat(q1) * Rᵀmat(q2) * Rᵀmat(ωbar(body2, Δt)) * Rmat(Quaternion(joint.vertices[2])) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


function setForce!(joint::Translational1, body1::Origin, body2::Body{T}, F::SVector{2,T}, No) where T
    F2 = joint.V12' * F
    τ2 = torqueFromForce(F2, joint.vertices[2])

    body2.F[No] = F2
    body2.τ[No] = τ2
    return
end

@inline function minimalCoordinates(joint::Translational1, body1::Origin, body2::Body, No)
    vertices = joint.vertices
    joint.V12 * (body2.x[No] + vrotate(vertices[2], body2.q[No]) - vertices[1])
end

@inline function g(joint::Translational1, body1::Origin, body2::Body, Δt, No)
    vertices = joint.vertices
    joint.V3 * (getx3(body2, Δt) + vrotate(vertices[2], getq3(body2, Δt)) - vertices[1])
end

@inline function ∂g∂posb(joint::Translational1{T}, body1::Origin, body2::Body, No) where T
    if body2.id == joint.cid
        q2 = body2.q[No]

        X = joint.V3

        R = joint.V3 * 2 * VRᵀmat(q2) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(q2)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Translational1{T}, body1::Origin, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        q2 = body2.q[No]

        V = Δt * joint.V3

        Ω = 2 * joint.V3 * VLmat(q2) * Rᵀmat(q2) * Rᵀmat(ωbar(body2, Δt)) * Rmat(Quaternion(joint.vertices[2])) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
