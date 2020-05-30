mutable struct Rotational{T,N} <: Joint{T,N}
    V12::SMatrix{2,3,T,6}
    V3::Adjoint{T,SVector{3,T}}
    qoff::Quaternion{T}

    F::SVector{3,T}
    τ::SVector{3,T}

    cid::Int64

    function Rotational{T,N}(body1::AbstractBody{T}, body2::AbstractBody{T}; axis::AbstractVector{T} = zeros(3), qoffset::Quaternion{T} = Quaternion{T}()) where {T,N}
        axis = vrotate(SVector(axis...), inv(qoffset))
        if norm(axis) != 0
            axis = axis / norm(axis)
        end
        A = svd(skew(axis)).Vt # in frame of body1
        V12 = A[1:2,:]
        V3 = axis' # instead of A[3,:] for correct sign: abs(axis) = abs(A[3,:])

        F = zeros(T,3)
        τ = zeros(T,3)

        cid = body2.id

        new{T,N}(V12, V3, qoffset, F, τ, cid), body1.id, body2.id
    end
end

Rotational0 = Rotational{T,0} where T
Rotational1 = Rotational{T,1} where T
Rotational2 = Rotational{T,2} where T
Rotational3 = Rotational{T,3} where T