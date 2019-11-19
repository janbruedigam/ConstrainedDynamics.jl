using StaticArrays
using LinearAlgebra
using Rotations

struct Quaternion{T} <: FieldVector{4, T}
    w::T
    v1::T
    v2::T
    v3::T
end


Quaternion(w::T, v1::T, v2::T, v3::T) where T = Quaternion{T}(w,v1,v2,v3)
Quaternion(w::T, v::AbstractVector{T}) where T = Quaternion{T}(w,v[1],v[2],v[3])
Quaternion(R::Rotation{3,T}) where T = Quaternion{T}(Quat(R).w,Quat(R).x,Quat(R).y,Quat(R).z)
Quaternion(v::AbstractVector{T}) where T = Quaternion{T}(0,v[1],v[2],v[3])
Quaternion(v::SVector{3,T}) where T = Quaternion{T}(0,v[1],v[2],v[3])
Quaternion(q::SVector{4,T}) where T = Quaternion{T}(q)
Quaternion{T}() where T = Quaternion{T}(1,0,0,0)
Quaternion() = Quaternion{Float64}()


skew(v::AbstractVector{T}) where T = SMatrix{3,3,T}(0,-v[3],v[2], v[3],0,-v[1], -v[2],v[1],0)

Vmat(T::Type) = SMatrix{3,4,T}(0,0,0, 1,0,0, 0,1,0, 0,0,1)
Vmat() = Vmat(Float64)
VTmat(T::Type) = Vmat(T)'
VTmat() = Vmat(Float64)'

Tmat(T::Type) = SMatrix{4,4,T}(1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,-1)
Tmat(::T) where T = Tmat(T)
Tmat() = Tmat(Float64)

Lmat(q::Quaternion{T}) where T = SMatrix{4,4,T}(q[1],q[2],q[3],q[4], -q[2], q[1],q[4],-q[3], -q[3], -q[4],q[1],q[2], -q[4], q[3],-q[2],q[1])
LTmat(q::Quaternion) = Lmat(q)'
Rmat(q::Quaternion{T})  where T = SMatrix{4,4,T}(q[1],q[2],q[3],q[4], -q[2], q[1],-q[4],q[3], -q[3], q[4],q[1],-q[2], -q[4], -q[3],q[2],q[1])
RTmat(q::Quaternion) = Rmat(q)'

rotate(v::Quaternion{T}, q::Quaternion{T}) where T = Vmat(T)*RTmat(q)*Lmat(q)*v

function angleAxis(q::Quaternion{T}) where T
    n = norm(Vmat(T)*q)
    if n == 0
        return zero(T), SVector{3,T}(0,0,0)
    else
        return 2*atan(n,q[1]), Vmat()*q./n
    end
end

srepeat(v::T,n) where T= T(repeat(v,n)...)
mrepeat(v::AbstractVector,n) = MVector(repeat(v,n)...)
sbracket(v::T) where T = SVector{1,T}(v)
mbracket(v::T) where T = MVector{1,T}(v)
