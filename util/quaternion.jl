# struct Quaternion{T <: Real} <: FieldVector{4, T}
#     s::T
#     v1::T
#     v2::T
#     v3::T
# end
#
#
# @inline Base.imag(q::Quaternion) = q[SVector{3}(2:4)]
#
#
# Quaternion(s::T, v1::T, v2::T, v3::T) where T = Quaternion{T}(s,v1,v2,v3)
# Quaternion(s::T, v::AbstractVector{T}) where T = Quaternion{T}(s,v[1],v[2],v[3])
# Quaternion(R::Rotation{3,T}) where T = Quaternion{T}(Quat(R).w,Quat(R).x,Quat(R).y,Quat(R).z)
# Quaternion(v::Vector{T}) where T = Quaternion{T}(0,v[1],v[2],v[3])
# Quaternion(v::SVector{3,T}) where T = Quaternion{T}(0,v[1],v[2],v[3])
# Quaternion(q::SVector{4,T}) where T = Quaternion{T}(q)
# Quaternion{T}() where T = Quaternion{T}(1,0,0,0)
#
#
# @inline skew(v::AbstractVector{T}) where T = SMatrix{3,3,T,9}(0,v[3],-v[2], -v[3],0,v[1], v[2],-v[1],0)
# @inline skewplusdiag(v::AbstractVector{T},w::T) where T = SMatrix{3,3,T,9}(w,v[3],-v[2], -v[3],w,v[1], v[2],-v[1],w)
#
# @inline Vmat(T::Type) = SMatrix{3,4,T,12}(0,0,0, 1,0,0, 0,1,0, 0,0,1)
# @inline Vmat() = Vmat(Float64)
# @inline VTmat(T::Type) = Vmat(T)'
# @inline VTmat() = Vmat(Float64)'
# @inline Vmat(q::SVector) = q[SVector{3}(2,3,4)]
# @inline Vmat(q::Quaternion) = q[SVector{3}(2,3,4)]
# @inline Vmat(A::SMatrix) = A[SVector{3}(2:4),:]
# @inline VTmat(A::SMatrix) = A[:,SVector{3}(2:4)]
#
#
# @inline Tmat(T::Type) = SMatrix{4,4,T,16}(1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,-1)
# @inline Tmat() = Tmat(Float64)
#
# @inline Lmat(q::Quaternion{T}) where T = SMatrix{4,4,T,16}(q.s,q.v1,q.v2,q.v3, -q.v1, q.s,q.v3,-q.v2, -q.v2, -q.v3,q.s,q.v1, -q.v3, q.v2,-q.v1,q.s)
# @inline LTmat(q::Quaternion{T}) where T = Lmat(q)'
# @inline Rmat(q::Quaternion{T})  where T = SMatrix{4,4,T,16}(q.s,q.v1,q.v2,q.v3, -q.v1, q.s,-q.v3,q.v2, -q.v2, q.v3,q.s,-q.v1, -q.v3, -q.v2,q.v1,q.s)
# @inline RTmat(q::Quaternion{T})  where T = Rmat(q)'
#
# @inline rotate(x::SVector, q::Quaternion) = q.s^2*x + 2*q.s*cross(imag(q),x) + 2*(imag(q)'*x)*imag(q) - (imag(q)'*imag(q))*x
#
# @inline function angleAxis(q::Quaternion{T}) where T
#     n = norm(Vmat(T)*q)
#     if n == 0
#         return zero(T), SVector{3,T}(0,0,0)
#     else
#         return 2*atan(n,q.s), Vmat()*q./n
#     end
# end

struct Quaternion{T<:Real} <: FieldVector{4,T}
    s::T
    v1::T
    v2::T
    v3::T
end

# Constructors
Quaternion(s::Real,v1::Real,v2::Real,v3::Real) = Quaternion(promote(s,v1,v2,v3)...)
Quaternion(s::Real) = Quaternion(s,0,0,0)
Quaternion(v::Vector) = Quaternion(0,v[1],v[2],v[3])
Quaternion(v::SVector{3,T}) where T = Quaternion(0,v[1],v[2],v[3])
Quaternion(s::T,v::SVector{3,T}) where T = Quaternion(s,v[1],v[2],v[3])
Quaternion(R::Rotation) = Quaternion(Quat(R).w,Quat(R).x,Quat(R).y,Quat(R).z)
Quaternion{T}() where T = Quaternion{T}(1,0,0,0)

# Basic quaternion operations
real(q::Quaternion) = q[1]
imag(q::Quaternion) = q[SUnitRange(2,4)]

conj(q::Quaternion) = Quaternion(q.s, -q.v1, -q.v2, -q.v3)
abs(q::Quaternion) = sqrt(q.s * q.s + q.v1 * q.v1 + q.v2 * q.v2 + q.v3 * q.v3)
abs2(q::Quaternion) = q.s * q.s + q.v1 * q.v1 + q.v2 * q.v2 + q.v3 * q.v3
Base.inv(q::Quaternion) = conj(q)

Base.:*(q1::Quaternion, q2::Quaternion) = Quaternion(  q1.s * q2.s - q1.v1 * q2.v1 - q1.v2 * q2.v2 - q1.v3 * q2.v3,
                                                       q1.s * q2.v1 + q1.v1 * q2.s + q1.v2 * q2.v3 - q1.v3 * q2.v2,
                                                       q1.s * q2.v2 - q1.v1 * q2.v3 + q1.v2 * q2.s + q1.v3 * q2.v1,
                                                       q1.s * q2.v3 + q1.v1 * q2.v2 - q1.v2 * q2.v1 + q1.v3 * q2.s)

Base.:/(q1::Quaternion, q2::Quaternion) = q1*inv(q2)
Base.:\(q1::Quaternion, q2::Quaternion) = inv(q1)*q2

angleaxis(q::Quaternion) = angle(q), axis(q)
angle(q::Quaternion) = 2*atan(sqrt(q.v1^2 + q.v2^2 + q.v3^2), q.s)
axis(q::Quaternion{T}) where T = q.s==1 ? SVector{3,T}(0,0,0) : SVector(q.v1, q.v2, q.v3)

qrotate(x::Quaternion,q::Quaternion) = q*x/q
vrotate(x::AbstractVector,q::Quaternion) = imag(qrotate(Quaternion(x),q))

# Matrix equivalences
𝟙(::Type{T}) where T = Quaternion(one(T))
𝟙() = 𝟙(Float64)
Vmat(::Type{T}) where T = SMatrix{3,4,T,12}(0,0,0, 1,0,0, 0,1,0, 0,0,1)
Vmat() = Vmat(Float64)
Vmat(q::SVector) = q[SUnitRange(2,4)]
Vmat(q::Quaternion) = imag(q)
Vᵀmat(::Type{T}) where T = SMatrix{4,3,T,12}(0,1,0,0, 0,0,1,0, 0,0,0,1)
Vᵀmat() = Vᵀmat(Float64)
Tmat(::Type{T}) where T = SMatrix{4,4,T,16}(1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,-1)
Tmat() = Tmat(Float64)

Lmat(q::Quaternion{T}) where T = SMatrix{4,4,T,16}(q.s,q.v1,q.v2,q.v3, -q.v1, q.s,q.v3,-q.v2, -q.v2, -q.v3,q.s,q.v1, -q.v3, q.v2,-q.v1,q.s)
Lᵀmat(q::Quaternion) = Lmat(q)'
Rmat(q::Quaternion{T})  where T = SMatrix{4,4,T,16}(q.s,q.v1,q.v2,q.v3, -q.v1, q.s,-q.v3,q.v2, -q.v2, q.v3,q.s,-q.v1, -q.v3, -q.v2,q.v1,q.s)
Rᵀmat(q::Quaternion) = Rmat(q)'

VLmat(q::Quaternion{T}) where T = SMatrix{3,4,T,12}(q.v1,q.v2,q.v3, q.s,q.v3,-q.v2, -q.v3,q.s,q.v1, q.v2,-q.v1,q.s)
VLᵀmat(q::Quaternion{T}) where T = SMatrix{3,4,T,12}(-q.v1,-q.v2,-q.v3, q.s,-q.v3,q.v2, q.v3,q.s,-q.v1, -q.v2,q.v1,q.s)
VRmat(q::Quaternion{T})  where T = SMatrix{3,4,T,12}(q.v1,q.v2,q.v3, q.s,-q.v3,q.v2, q.v3,q.s,-q.v1, -q.v2,q.v1,q.s)
VRᵀmat(q::Quaternion{T})  where T = SMatrix{3,4,T,12}(-q.v1,-q.v2,-q.v3, q.s,q.v3,-q.v2, -q.v3,q.s,q.v1, q.v2,-q.v1,q.s)

LVᵀmat(q::Quaternion{T}) where T = SMatrix{4,3,T,12}(-q.v1, q.s,q.v3,-q.v2, -q.v2, -q.v3,q.s,q.v1, -q.v3, q.v2,-q.v1,q.s)
LᵀVᵀmat(q::Quaternion{T}) where T = SMatrix{4,3,T,12}(q.v1, q.s,-q.v3,q.v2, q.v2, q.v3,q.s,-q.v1, q.v3, -q.v2,q.v1,q.s)
RVᵀmat(q::Quaternion{T}) where T = SMatrix{4,3,T,12}(-q.v1, q.s,-q.v3,q.v2, -q.v2, q.v3,q.s,-q.v1, -q.v3, -q.v2,q.v1,q.s)
RᵀVᵀmat(q::Quaternion{T}) where T = SMatrix{4,3,T,12}(q.v1, q.s,q.v3,-q.v2, q.v2, -q.v3,q.s,q.v1, q.v3, q.v2,-q.v1,q.s)

skewplusdiag(v::AbstractVector{T},w::T) where T = SMatrix{3,3,T,9}(w,v[3],-v[2], -v[3],w,v[1], v[2],-v[1],w)
