mutable struct Rotational{T,N} <: Joint{T,N}
    V12::SMatrix{2,3,T,6}
    V3::Adjoint{T,SVector{3,T}}
    qoff::UnitQuaternion{T}

    F::SVector{3,T}
    τ::SVector{3,T}

    childid::Int64

    function Rotational{T,N}(body1::AbstractBody{T}, body2::AbstractBody{T}; 
            axis::AbstractVector{T} = zeros(3), qoffset::UnitQuaternion{T} = one(UnitQuaternion{T})
        ) where {T,N}
        
        axis = vrotate(SVector(axis...), inv(qoffset))
        if norm(axis) != 0
            axis = axis / norm(axis)
        end
        A = svd(skew(axis)).Vt # in frame of body1
        V12 = A[1:2,:]
        V3 = axis' # instead of A[3,:] for correct sign: abs(axis) = abs(A[3,:])

        F = zeros(T,3)
        τ = zeros(T,3)

        childid = body2.id

        new{T,N}(V12, V3, qoffset, F, τ, childid), body1.id, body2.id
    end
end

Rotational0 = Rotational{T,0} where T
Rotational1 = Rotational{T,1} where T
Rotational2 = Rotational{T,2} where T
Rotational3 = Rotational{T,3} where T


# Position level constraints

@inline function g(joint::Rotational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return Vmat(joint.qoff \ (qa \ qb))
end
@inline function g(joint::Rotational, xb::AbstractVector, qb::UnitQuaternion)
    return Vmat(joint.qoff \ qb)
end


# Naive derivatives without accounting for quaternion specialness

@inline function ∂g∂posa(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    X = @SMatrix zeros(T, 3, 3)
    Q = VLᵀmat(joint.qoff) * Rmat(qb) * Tmat()

    return X, Q
end
@inline function ∂g∂posb(joint::Rotational{T}, qa::UnitQuaternion, qb::UnitQuaternion) where T
    X = @SMatrix zeros(T, 3, 3)
    Q = VLᵀmat(joint.qoff) * Lᵀmat(qa)

    return X, Q
end
@inline function ∂g∂posb(joint::Rotational{T}, qb::UnitQuaternion) where T
    X = @SMatrix zeros(T, 3, 3)
    Q = VLᵀmat(joint.qoff)

    return X, Q
end