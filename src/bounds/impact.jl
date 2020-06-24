mutable struct Impact{T} <: Contact{T}
    Nx::Adjoint{T,SVector{6,T}}
    offset::SVector{6,T}


    function Impact(body::Body{T}, normal::AbstractVector; offset::AbstractVector = zeros(3)) where T
        normal = normal / norm(normal)

        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        A = Array(svd(skew(normal)).V) # gives two plane vectors
        A[:,3] = normal # to ensure correct sign
        Ainv = inv(A)
        ainv3 = Ainv[3,:]
        Nx = [ainv3;0;0;0]'
        offset = [offset;0;0;0]

        new{T}(Nx, offset), body.id
    end
end
