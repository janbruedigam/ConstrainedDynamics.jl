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

# vec(G) Jacobian
@inline function ∂2g∂posaa(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    VLtqa = VLᵀmat(qa)
    VRqa = VRmat(qa)
    R1, R2, R3, R4 = ∂R∂qsplit(T)
    Lt1, Lt2, Lt3, Lt4 = ∂Lᵀ∂qsplit(T)
    L1, L2, L3, L4 = ∂L∂qsplit(T)
    L1 = L1[:,2:4]
    L2 = L2[:,2:4]
    L3 = L3[:,2:4]
    L4 = L4[:,2:4]
    Ltpos = Lmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))

    Q1 = -VLtqa * R2
    Q2 = -VLtqa * R3
    Q3 = -VLtqa * R4

    Q4 = -VRqa * Lt2
    Q5 = -VRqa * Lt3
    Q6 = -VRqa * Lt4

    Q7 = -2*VLtqa * L1
    Q8 = -2*VLtqa * L2
    Q9 = -2*VLtqa * L3
    Q10 = -2*VLtqa * L4

    XX = szeros(T, 9, 3)
    XQ = [Q1; Q2; Q3] + [Q4; Q5; Q6]
    QX = [Q7; Q8; Q9; Q10]
    QQ = kron(Ltpos,Vmat()) * [Lt1; Lt2; Lt3; Lt4]

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posab(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    VLtqa = VLᵀmat(qa)
    Lr1, Lr2, Lr3, Lr4 = ∂L∂qsplit(T)
    Lr1 = Lr1[:,2:4]
    Lr2 = Lr2[:,2:4]
    Lr3 = Lr3[:,2:4]
    Lr4 = Lr4[:,2:4]
    VLtqaLqbLpb = VLᵀmat(qa)*Lmat(qb)*Lmat(UnitQuaternion(joint.vertices[2]))
    VLtqaRtqbRpb = VLᵀmat(qa)*Rᵀmat(qb)*Rmat(UnitQuaternion(joint.vertices[2]))
    Lt1, Lt2, Lt3, Lt4 = ∂Lᵀ∂qsplit(T)
    L1, L2, L3, L4 = ∂L∂qsplit(T)

    Q1 = 2*VLtqa * Lr1
    Q2 = 2*VLtqa * Lr2
    Q3 = 2*VLtqa * Lr3
    Q4 = 2*VLtqa * Lr4

    Q5 = 2*VLtqaLqbLpb * Lt1
    Q6 = 2*VLtqaLqbLpb * Lt2
    Q7 = 2*VLtqaLqbLpb * Lt3
    Q8 = 2*VLtqaLqbLpb * Lt4

    Q9 = 2*VLtqaRtqbRpb * L1
    Q10 = 2*VLtqaRtqbRpb * L2
    Q11 = 2*VLtqaRtqbRpb * L3
    Q12 = 2*VLtqaRtqbRpb * L4

    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = [Q1; Q2; Q3; Q4]
    QQ = [Q5; Q6; Q7; Q8] + [Q9; Q10; Q11; Q12]

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posba(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    VLtqa = VLᵀmat(qa)
    VRqa = VRmat(qa)
    R1, R2, R3, R4 = ∂R∂qsplit(T)
    Lt1, Lt2, Lt3, Lt4 = ∂Lᵀ∂qsplit(T)
    RtpbRqb = Rᵀmat(UnitQuaternion(joint.vertices[2]))*Rmat(qb)
    RtpbRqbRtqa = Rᵀmat(UnitQuaternion(joint.vertices[2]))*Rmat(qb)*Rᵀmat(qa)
    
    Q1 = VLtqa * R2
    Q2 = VLtqa * R3
    Q3 = VLtqa * R4

    Q4 = VRqa * Lt2
    Q5 = VRqa * Lt3
    Q6 = VRqa * Lt4

    XX = szeros(T, 9, 3)
    XQ = [Q1; Q2; Q3] + [Q4; Q5; Q6]
    QX = szeros(T, 12, 3)
    QQ = kron(RtpbRqb,2*VLᵀmat(qa))*[R1;R2;R3;R4] + kron(RtpbRqbRtqa,2*Vmat())*[Lt1;Lt2;Lt3;Lt4]

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    Rt1, Rt2, Rt3, Rt4 = ∂Rᵀ∂qsplit(T)
    VLtqaRqa = VLᵀmat(qa)*Rmat(qa)

    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 12, 3)
    QQ = kron(Rᵀmat(UnitQuaternion(joint.vertices[2])),2*VLtqaRqa)*[Rt1;Rt2;Rt3;Rt4]


    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Translational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
    Rt1, Rt2, Rt3, Rt4 = ∂Rᵀ∂qsplit(T)

    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 12, 3)
    QQ = kron(Rᵀmat(UnitQuaternion(joint.vertices[2])),2*Vmat())*[Rt1;Rt2;Rt3;Rt4]

    return XX, XQ, QX, QQ
end