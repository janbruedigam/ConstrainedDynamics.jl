
mutable struct Impact{T} <: Contact{T}
    ainv3::Adjoint{T,SVector{3,T}}
    p::SVector{3,T}
    offset::SVector{3,T}
    

    function Impact(body::Body{T}, normal::AbstractVector; p::AbstractVector = zeros(3),  offset::AbstractVector = zeros(3)) where T
        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'

        new{T}(ainv3, p, offset), body.id
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Impact{T}) where {T}
    summary(io, constraint)
    println(io,"")
    println(io, " ainv3:  "*string(constraint.ainv3))
    println(io, " offset: "*string(constraint.offset))
end


## Additional force for friction default
@inline additionalforce(bound::Impact{T}, x::AbstractVector, q::UnitQuaternion) where T = szeros(T, 6)