
mutable struct Impact{T} <: Contact{T}
    Nx::Adjoint{T,SVector{3,T}}
    p::SVector{3,T}
    offset::SVector{6,T}
    

    function Impact(body::Body{T}, normal::AbstractVector; p::AbstractVector = zeros(3),  offset::AbstractVector = zeros(3)) where T
        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]
        Nx = ainv3'
        offset = [offset;0.0;0.0;0.0]

        new{T}(Nx, p, offset), body.id
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Impact{T}) where {T}
    summary(io, constraint)
    println(io,"")
    println(io, " Nx:     "*string(constraint.Nx))
    println(io, " offset: "*string(constraint.offset))
end