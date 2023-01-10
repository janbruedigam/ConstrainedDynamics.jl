mutable struct Impact{T,N} <: Bound{T,N}
    ainv3::Adjoint{T,SVector{3,T}}
    p::SVector{3,T}
    offset::SVector{3,T}
    

    function Impact(body::Body{T}, normal::AbstractVector; p::AbstractVector = zeros(3),  offset::AbstractVector = zeros(3)) where T
        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'

        new{T,2}(ainv3, p, offset), body.id, nothing
    end
end


### Constraints and derivatives
## Position level constraints (for dynamics)
@inline g(impact::Impact{T}, x::AbstractVector, q::Quaternion) where T = SVector{1,T}(impact.ainv3 * (x + vrotate(impact.p,q) - impact.offset))

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂pos(impact::Impact, x::AbstractVector, q::Quaternion)
    p = impact.p
    X = impact.ainv3
    Q = impact.ainv3 * (VLmat(q) * Lmat(Quaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(Quaternion(p)))
    return X, Q
end
