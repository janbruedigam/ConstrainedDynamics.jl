mutable struct Rotational{T,Nc} <: Joint{T,Nc}
    V12::SMatrix{2,3,T,6}
    V3::Adjoint{T,SVector{3,T}}
    qoff::Quaternion{T}
    cid::Int64

    function Rotational{T,Nc}(body1::AbstractBody{T}, body2::AbstractBody{T}; axis::AbstractVector{T}=zeros(3), offset::Quaternion{T} = Quaternion{T}()) where {T,Nc}
        axis = vrotate(SVector(axis...),inv(offset))
        if norm(axis) != 0
            axis = axis / norm(axis)
        end
        A = svd(skew(axis)).Vt # in frame of body1
        V12 = A[1:2,:]
        V3 = axis' # instead of A[3,:] for correct sign: abs(axis) = abs(A[3,:])
        cid = body2.id

        new{T,Nc}(V12, V3, offset, cid), body1.id, body2.id
    end
end

Rotational0 = Rotational{T,0} where T
Rotational1 = Rotational{T,1} where T
Rotational2 = Rotational{T,2} where T
Rotational3 = Rotational{T,3} where T