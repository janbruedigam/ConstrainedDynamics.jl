mutable struct Friction{T} <: Contact{T}
    Nx::Adjoint{T,SVector{6,T}}
    D::SMatrix{2,6,T,12}
    cf::T
    b::SVector{2,T}
    offset::SVector{6,T}


    function Friction(body::Body{T}, normal::AbstractVector{T}, cf::T;offset::AbstractVector{T} = zeros(3)) where T
        @assert cf>0
        normal = normal / norm(normal)

        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        A = Array(svd(skew(normal)).V) # gives two plane vectors
        A[:,3] = normal # to ensure correct sign
        Ainv = inv(A)
        ainv3 = Ainv[3,:]
        Nx = [ainv3;0;0;0]'
        D = [A[:,1:2];zeros(3,2)]'
        offset = [offset;0;0;0]

        new{T}(Nx, D, cf, zeros(2),offset), body.id
    end
end


@inline additionalforce(friction::Friction) = friction.D'*friction.b
@inline function calcFrictionForce!(mechanism, ineqc, friction::Friction, i, body::Body)
    calcFrictionForce!(mechanism, friction, body, ineqc.γsol[2][i])
    return
end