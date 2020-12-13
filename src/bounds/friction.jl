mutable struct Friction{T} <: Contact{T}
    Nx::Adjoint{T,SVector{3,T}}
    p::SVector{3,T}
    D::SMatrix{2,6,T,12}
    cf::T
    offset::SVector{3,T}
    b::SVector{2,T}
    

    function Friction(body::Body{T}, normal::AbstractVector, cf::Real; p::AbstractVector = zeros(3), offset::AbstractVector = zeros(3)) where T
        @assert cf>0

        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]
        Nx = ainv3'
        D = [[V1 V2];szeros(3,2)]'

        new{T}(Nx, p, D, cf, offset, zeros(2)), body.id
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Friction{T}) where {T}
    summary(io, constraint)
    println(io,"")
    println(io, " Nx:     "*string(constraint.Nx))
    println(io, " D:      "*string(constraint.D))
    println(io, " cf:     "*string(constraint.cf))
    println(io, " offset: "*string(constraint.offset))
end


@inline additionalforce(friction::Friction) = friction.D'*friction.b
@inline function calcFrictionForce!(mechanism, ineqc, friction::Friction, i, body::Body)
    calcFrictionForce!(mechanism, friction, body, ineqc.γsol[2][i])
    return
end