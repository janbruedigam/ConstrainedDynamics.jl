Rotations.UnitQuaternion(w::T, v::StaticVector{3,T}, normalize::Bool = true) where T = UnitQuaternion{T}(w, v[1], v[2], v[3], normalize)
Rotations.UnitQuaternion(w::T, v::Vector{T}, normalize::Bool = true) where T = (@assert length(v)==3; UnitQuaternion{T}(w, v[1], v[2], v[3], normalize))
Rotations.UnitQuaternion(v::StaticVector{3,T}) where T = pure_quaternion(v)
Rotations.UnitQuaternion(v::Vector) = (@assert length(v)==3; pure_quaternion(v))

@inline imag(q::UnitQuaternion) = Rotations.vector(q)

qrotate(q1::UnitQuaternion,q2::UnitQuaternion) = q2 * q1 / q2
vrotate(v::AbstractVector,q::UnitQuaternion) = imag(qrotate(pure_quaternion(v), q))

Vmat(::Type{T}) where T = SMatrix{3,4,T,12}(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
Vmat() = Vmat(Float64)
Vmat(q::SVector) = q[SUnitRange(2, 4)]
Vmat(A::SMatrix) = A[SUnitRange(2, 4),:]
Vmat(q::UnitQuaternion) = imag(q)
Vᵀmat(::Type{T}) where T = SMatrix{4,3,T,12}(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
Vᵀmat() = Vᵀmat(Float64)
Tmat(::Type{T}) where T = SMatrix{4,4,T,16}(1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1)
Tmat() = Tmat(Float64)

Lmat(q::UnitQuaternion{T}) where T = SMatrix{4,4,T,16}(q.w, q.x, q.y, q.z, -q.x, q.w, q.z, -q.y, -q.y, -q.z, q.w, q.x, -q.z, q.y, -q.x, q.w)
Lᵀmat(q) = Lmat(q)'
Rmat(q::UnitQuaternion{T})  where T = SMatrix{4,4,T,16}(q.w, q.x, q.y, q.z, -q.x, q.w, -q.z, q.y, -q.y, q.z, q.w, -q.x, -q.z, -q.y, q.x, q.w)
Rᵀmat(q) = Rmat(q)'

VLmat(q::UnitQuaternion{T}) where T = SMatrix{3,4,T,12}(q.x, q.y, q.z, q.w, q.z, -q.y, -q.z, q.w, q.x, q.y, -q.x, q.w)
VLᵀmat(q::UnitQuaternion{T}) where T = SMatrix{3,4,T,12}(-q.x, -q.y, -q.z, q.w, -q.z, q.y, q.z, q.w, -q.x, -q.y, q.x, q.w)
VRmat(q::UnitQuaternion{T})  where T = SMatrix{3,4,T,12}(q.x, q.y, q.z, q.w, -q.z, q.y, q.z, q.w, -q.x, -q.y, q.x, q.w)
VRᵀmat(q::UnitQuaternion{T})  where T = SMatrix{3,4,T,12}(-q.x, -q.y, -q.z, q.w, q.z, -q.y, -q.z, q.w, q.x, q.y, -q.x, q.w)

LVᵀmat(q::UnitQuaternion{T}) where T = SMatrix{4,3,T,12}(-q.x, q.w, q.z, -q.y, -q.y, -q.z, q.w, q.x, -q.z, q.y, -q.x, q.w)
LᵀVᵀmat(q::UnitQuaternion{T}) where T = SMatrix{4,3,T,12}(q.x, q.w, -q.z, q.y, q.y, q.z, q.w, -q.x, q.z, -q.y, q.x, q.w)
RVᵀmat(q::UnitQuaternion{T}) where T = SMatrix{4,3,T,12}(-q.x, q.w, -q.z, q.y, -q.y, q.z, q.w, -q.x, -q.z, -q.y, q.x, q.w)
RᵀVᵀmat(q::UnitQuaternion{T}) where T = SMatrix{4,3,T,12}(q.x, q.w, q.z, -q.y, q.y, -q.z, q.w, q.x, q.z, q.y, -q.x, q.w)

angle(q) = rotation_angle(AngleAxis(q))
axis(q) = rotation_axis(AngleAxis(q))

# Rotations.:*(A::AbstractMatrix,q::UnitQuaternion) = A*params(q)
# Rotations.:*(A::Adjoint{T,<:AbstractArray{T,1}}, q::UnitQuaternion{T}) where T = A*params(q)



# struct UnitQuaternion{T <: Real} <: FieldVector{4,T}
#     s::T
#     v1::T
#     v2::T
#     v3::T
# end

# # Constructors
# UnitQuaternion(s::Real,v1::Real,v2::Real,v3::Real) = UnitQuaternion(promote(s, v1, v2, v3)...)
# UnitQuaternion(s::Real) = UnitQuaternion(s, 0, 0, 0)
# UnitQuaternion(v::Vector) = (@assert length(v)==3; UnitQuaternion(0, v[1], v[2], v[3]))
# UnitQuaternion(s::T,v::Vector{T}) where T = (@assert length(v)==3; UnitQuaternion(s, v[1], v[2], v[3]))
# UnitQuaternion(v::SVector{3,T}) where T = UnitQuaternion(0, v[1], v[2], v[3])
# UnitQuaternion(s::T,v::SVector{3,T}) where T = UnitQuaternion(s, v[1], v[2], v[3])
# UnitQuaternion(R::Rotation) = (q = UnitQuaternion(R); UnitQuaternion(q.w, q.x, q.y, q.z))
# one(UnitQuaternion{T}) where T = UnitQuaternion{T}(1, 0, 0, 0)

# # Basic quaternion operations
# LinearAlgebra.real(q::UnitQuaternion) = q[1]
# LinearAlgebra.imag(q::UnitQuaternion) = q[SUnitRange(2, 4)]

# LinearAlgebra.conj(q::UnitQuaternion) = UnitQuaternion(q.s, -q.v1, -q.v2, -q.v3)
# Base.abs(q::UnitQuaternion) = sqrt(q.s * q.s + q.v1 * q.v1 + q.v2 * q.v2 + q.v3 * q.v3)
# Base.abs2(q::UnitQuaternion) = q.s * q.s + q.v1 * q.v1 + q.v2 * q.v2 + q.v3 * q.v3
# Base.inv(q::UnitQuaternion) = conj(q)

# Base.:*(q1::UnitQuaternion, q2::UnitQuaternion) = UnitQuaternion(  q1.s * q2.s - q1.v1 * q2.v1 - q1.v2 * q2.v2 - q1.v3 * q2.v3,
#                                                        q1.s * q2.v1 + q1.v1 * q2.s + q1.v2 * q2.v3 - q1.v3 * q2.v2,
#                                                        q1.s * q2.v2 - q1.v1 * q2.v3 + q1.v2 * q2.s + q1.v3 * q2.v1,
#                                                        q1.s * q2.v3 + q1.v1 * q2.v2 - q1.v2 * q2.v1 + q1.v3 * q2.s)
# Base.:*(q::UnitQuaternion, x::Number) = UnitQuaternion(q.s * x, q.v1 * x, q.v2 * x, q.v3 * x)
# Base.:*(x::Number, q::UnitQuaternion) = q * x

# Base.:/(q1::UnitQuaternion, q2::UnitQuaternion) = q1 * inv(q2)
# Base.:\(q1::UnitQuaternion, q2::UnitQuaternion) = inv(q1) * q2

# Base.:-(q::UnitQuaternion) = UnitQuaternion(-q.s,-q.v1,-q.v2,-q.v3)

# angleaxis(q::UnitQuaternion) = angle(q), axis(q)
# angle(q::UnitQuaternion) = 2 * atan(sqrt(q.v1^2 + q.v2^2 + q.v3^2), q.s) # This is the rotation angle φ (φ = 2θ)
# function axis(q::UnitQuaternion{T}) where T
#     qv = SVector(q.v1, q.v2, q.v3)
#     if norm(qv) == 0
#         return SVector{3,T}(0, 0, 0)
#     else 
#         return qv / norm(qv)
#     end
# end

# qrotate(x::UnitQuaternion,q::UnitQuaternion) = q * x / q
# vrotate(x::AbstractVector,q::UnitQuaternion) = imag(qrotate(UnitQuaternion(x), q))

# # Matrix equivalences
# # 𝟙(::Type{T}) where T = UnitQuaternion(one(T))
# # 𝟙() = 𝟙(Float64)
# Vmat(::Type{T}) where T = SMatrix{3,4,T,12}(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
# Vmat() = Vmat(Float64)
# Vmat(q::SVector) = q[SUnitRange(2, 4)]
# Vmat(A::SMatrix) = A[SUnitRange(2, 4),:]
# Vmat(q::UnitQuaternion) = imag(q)
# Vᵀmat(::Type{T}) where T = SMatrix{4,3,T,12}(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
# Vᵀmat() = Vᵀmat(Float64)
# Tmat(::Type{T}) where T = SMatrix{4,4,T,16}(1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1)
# Tmat() = Tmat(Float64)

# Lmat(q::UnitQuaternion{T}) where T = SMatrix{4,4,T,16}(q.s, q.v1, q.v2, q.v3, -q.v1, q.s, q.v3, -q.v2, -q.v2, -q.v3, q.s, q.v1, -q.v3, q.v2, -q.v1, q.s)
# Lᵀmat(q::UnitQuaternion) = Lmat(q)'
# Rmat(q::UnitQuaternion{T})  where T = SMatrix{4,4,T,16}(q.s, q.v1, q.v2, q.v3, -q.v1, q.s, -q.v3, q.v2, -q.v2, q.v3, q.s, -q.v1, -q.v3, -q.v2, q.v1, q.s)
# Rᵀmat(q::UnitQuaternion) = Rmat(q)'

# VLmat(q::UnitQuaternion{T}) where T = SMatrix{3,4,T,12}(q.v1, q.v2, q.v3, q.s, q.v3, -q.v2, -q.v3, q.s, q.v1, q.v2, -q.v1, q.s)
# VLᵀmat(q::UnitQuaternion{T}) where T = SMatrix{3,4,T,12}(-q.v1, -q.v2, -q.v3, q.s, -q.v3, q.v2, q.v3, q.s, -q.v1, -q.v2, q.v1, q.s)
# VRmat(q::UnitQuaternion{T})  where T = SMatrix{3,4,T,12}(q.v1, q.v2, q.v3, q.s, -q.v3, q.v2, q.v3, q.s, -q.v1, -q.v2, q.v1, q.s)
# VRᵀmat(q::UnitQuaternion{T})  where T = SMatrix{3,4,T,12}(-q.v1, -q.v2, -q.v3, q.s, q.v3, -q.v2, -q.v3, q.s, q.v1, q.v2, -q.v1, q.s)

# LVᵀmat(q::UnitQuaternion{T}) where T = SMatrix{4,3,T,12}(-q.v1, q.s, q.v3, -q.v2, -q.v2, -q.v3, q.s, q.v1, -q.v3, q.v2, -q.v1, q.s)
# LᵀVᵀmat(q::UnitQuaternion{T}) where T = SMatrix{4,3,T,12}(q.v1, q.s, -q.v3, q.v2, q.v2, q.v3, q.s, -q.v1, q.v3, -q.v2, q.v1, q.s)
# RVᵀmat(q::UnitQuaternion{T}) where T = SMatrix{4,3,T,12}(-q.v1, q.s, -q.v3, q.v2, -q.v2, q.v3, q.s, -q.v1, -q.v3, -q.v2, q.v1, q.s)
# RᵀVᵀmat(q::UnitQuaternion{T}) where T = SMatrix{4,3,T,12}(q.v1, q.s, q.v3, -q.v2, q.v2, -q.v3, q.s, q.v1, q.v3, q.v2, -q.v1, q.s)

# skewplusdiag(v::AbstractVector{T},w::T) where T = SMatrix{3,3,T,9}(w, v[3], -v[2], -v[3], w, v[1], v[2], -v[1], w)

# function slerp(q1,q2,h)
#     s = q1'*q2
#     if s < 0
#         s = -s
#         q2 = -q2
#     end

#     qdiff = q1\q2
#     φdiff, udiff = angleaxis(qdiff)
#     φint = φdiff*h
#     qint = UnitQuaternion(cos(φint/2),udiff*sin(φint/2))
    
#     return q1*qint
# end