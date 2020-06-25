mutable struct Friction{T} <: Contact{T}
    Nx::Adjoint{T,SVector{6,T}}
    D::SMatrix{2,6,T,12}
    cf::T
    b::SVector{2,T}
    offset::SVector{6,T}


    function Friction(body::Body{T}, normal::AbstractVector, cf::Real; offset::AbstractVector = zeros(3)) where T
        @assert cf>0

        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]
        Nx = [ainv3;0.0;0.0;0.0]'
        D = [[V1 V2];szeros(3,2)]'
        offset = [offset;0.0;0.0;0.0]

        new{T}(Nx, D, cf, zeros(2),offset), body.id
    end
end


@inline additionalforce(friction::Friction) = friction.D'*friction.b
@inline function calcFrictionForce!(mechanism, ineqc, friction::Friction, i, body::Body)
    calcFrictionForce!(mechanism, friction, body, ineqc.γsol[2][i])
    return
end