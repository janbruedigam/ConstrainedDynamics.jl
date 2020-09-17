mutable struct Translational{T,N} <: Joint{T,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    vertices::NTuple{2,SVector{3,T}} # in body1's & body2's frames
    
    F::SVector{3,T}
    τ::SVector{3,T}

    childid::Int64

    function Translational{T,N}(body1::AbstractBody, body2::AbstractBody;
            p1::AbstractVector = zeros(3), p2::AbstractVector = zeros(3), axis::AbstractVector = zeros(3)
        ) where {T,N}
        
        vertices = (p1, p2)
        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]

        F = zeros(T,3)
        τ = zeros(T,3)

        childid = body2.id

        new{T,N}(V3, V12, vertices, F, τ, childid), body1.id, body2.id
    end
end

Translational0 = Translational{T,0} where T
Translational1 = Translational{T,1} where T
Translational2 = Translational{T,2} where T
Translational3 = Translational{T,3} where T


@inline function getPositionDelta(joint::Translational, body1::AbstractBody, body2::Body, x::SVector)
    Δx = zerodimstaticadjoint(nullspacemat(joint)) * x # in body1 frame
    return Δx
end
@inline function getVelocityDelta(joint::Translational, body1::AbstractBody, body2::Body, v::SVector)
    Δv = zerodimstaticadjoint(nullspacemat(joint)) * v # in body1 frame
    return Δv
end

@inline function minimalCoordinates(joint::Translational, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return nullspacemat(joint) * g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
end
@inline function minimalCoordinates(joint::Translational, body1::Origin, body2::Body)
    stateb = body2.state
    return nullspacemat(joint) * g(joint, stateb.xc, stateb.qc)
end

# Position level constraints

@inline function g(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
end
@inline function g(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return xb + vrotate(vertices[2], qb) - vertices[1]
end


# Derivatives NOT accounting for quaternion specialness

@inline function ∂g∂posa(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    point2 = xb + vrotate(joint.vertices[2], qb)

    X = -VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * (Lmat(UnitQuaternion(point2)) - Lmat(UnitQuaternion(xa)))

    return X, Q
end
@inline function ∂g∂posb(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X = VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * Rmat(qa) * Rᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))

    return X, Q
end
@inline function ∂g∂posb(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
    X = I
    Q = 2 * VRᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))

    return X, Q
end


# vec(G) Jacobian (also NOT accounting for quaternion specialness)

@inline function ∂2g∂posaa(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    Lpos = Lmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))
    Ltpos = Lᵀmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))

    XX = szeros(T, 9, 3)
    XQ = -kron(Vmat(T),VLᵀmat(qa))*∂R∂qsplit(T) - kron(VRᵀmat(qa),Vmat(T))*∂Lᵀ∂qsplit(T)
    QX = -kron(VLᵀmat(qa),2*VLᵀmat(qa))*∂L∂qsplit(T)[:,SA[2; 3; 4]]
    QQ = kron(Vmat(T),2*VLᵀmat(qa)*Lpos)*∂L∂qsplit(T) + kron(VLᵀmat(qa)*Ltpos,2*Vmat(T))*∂Lᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posab(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = kron(VLᵀmat(qa),2*VLᵀmat(qa))*∂L∂qsplit(T)[:,SA[2; 3; 4]]
    QQ = kron(VLᵀmat(qa),2*VLᵀmat(qa)*Lmat(qb)*Lmat(UnitQuaternion(joint.vertices[2])))*∂Lᵀ∂qsplit(T) + kron(VLᵀmat(qa)*Lmat(qb)*Lᵀmat(UnitQuaternion(joint.vertices[2])),2*VLᵀmat(qa))*∂L∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posba(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = kron(Vmat(T),VLᵀmat(qa))*∂R∂qsplit(T) + kron(VRᵀmat(qa),Vmat(T))*∂Lᵀ∂qsplit(T)
    QX = szeros(T, 9, 3)
    QQ = kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2]))*Rmat(qb),2*VLᵀmat(qa))*∂R∂qsplit(T) + kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2]))*Rmat(qb)*Rᵀmat(qa),2*Vmat(T))*∂Lᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(Vmat(T),2*VLᵀmat(qa)*Rmat(qa)*Rᵀmat(qb)*Rmat(UnitQuaternion(joint.vertices[2])))*∂L∂qsplit(T) + kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2])),2*VLᵀmat(qa)*Rmat(qa))*∂Rᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Translational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(Vmat(T),2*VRᵀmat(qb)*Rmat(UnitQuaternion(joint.vertices[2])))*∂L∂qsplit(T) + kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2])),2*Vmat(T))*∂Rᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end