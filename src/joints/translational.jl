mutable struct Translational{T,N} <: Joint{T,N}
    vertices::NTuple{2,SVector{3,T}}
    V12::SMatrix{2,3,T,6}
    V3::Adjoint{T,SVector{3,T}}

    F::SVector{3,T}
    τ::SVector{3,T}

    cid::Int64

    function Translational{T,N}(body1::AbstractBody{T}, body2::AbstractBody{T};
        p1::AbstractVector{T} = zeros(3), p2::AbstractVector{T} = zeros(3), axis::AbstractVector{T} = zeros(3)) where {T,N}
        
        vertices = (p1, p2)
        if norm(axis) != 0
            axis = axis / norm(axis)
        end
        A = svd(skew(axis)).Vt # in frame of body1
        V12 = A[1:2,:]
        V3 = axis' # instead of A[3,:] for correct sign: abs(axis) = abs(A[3,:])

        F = zeros(T,3)
        τ = zeros(T,3)

        cid = body2.id

        new{T,N}(vertices, V12, V3, F, τ, cid), body1.id, body2.id
    end
end

Translational0 = Translational{T,0} where T
Translational1 = Translational{T,1} where T
Translational2 = Translational{T,2} where T
Translational3 = Translational{T,3} where T
