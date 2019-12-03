struct Quaternion{T} <: FieldVector{4, T}
    w::T
    v1::T
    v2::T
    v3::T
end

function Base.getproperty(Q::Quaternion,x::Symbol)
    if x == :v
        return Q[SVector{3}(2:4)]
    else
        return getfield(Q,x)
    end
end


Quaternion(w::T, v1::T, v2::T, v3::T) where T = Quaternion{T}(w,v1,v2,v3)
Quaternion(w::T, v::AbstractVector{T}) where T = Quaternion{T}(w,v[1],v[2],v[3])
Quaternion(R::Rotation{3,T}) where T = Quaternion{T}(Quat(R).w,Quat(R).x,Quat(R).y,Quat(R).z)
Quaternion(v::Vector{T}) where T = Quaternion{T}(0,v[1],v[2],v[3])
Quaternion(v::SVector{3,T}) where T = Quaternion{T}(0,v[1],v[2],v[3])
Quaternion(q::SVector{4,T}) where T = Quaternion{T}(q)
Quaternion{T}() where T = Quaternion{T}(1,0,0,0)
Quaternion() = Quaternion{Float64}()


@inline skew(v::AbstractVector{T}) where T = SMatrix{3,3,T,9}(0,v[3],-v[2], -v[3],0,v[1], v[2],-v[1],0)
@inline skewplusdiag(v::AbstractVector{T},w::T) where T = SMatrix{3,3,T,9}(w,v[3],-v[2], -v[3],w,v[1], v[2],-v[1],w)

@inline Vmat(T::Type) = SMatrix{3,4,T,12}(0,0,0, 1,0,0, 0,1,0, 0,0,1)
@inline Vmat() = Vmat(Float64)
@inline VTmat(T::Type) = Vmat(T)'
@inline VTmat() = Vmat(Float64)'
@inline Vmat(q::SVector) = q[SVector{3}(2,3,4)]
@inline Vmat(q::Quaternion) = q[SVector{3}(2,3,4)]
@inline Vmat(A::SMatrix) = A[SVector{3}(2:4),:]
@inline VTmat(A::SMatrix) = A[:,SVector{3}(2:4)]


@inline Tmat(T::Type) = SMatrix{4,4,T,16}(1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,-1)
@inline Tmat() = Tmat(Float64)

@inline Lmat(q::Quaternion{T}) where T = SMatrix{4,4,T,16}(q[1],q[2],q[3],q[4], -q[2], q[1],q[4],-q[3], -q[3], -q[4],q[1],q[2], -q[4], q[3],-q[2],q[1])
@inline LTmat(q::Quaternion) = Lmat(q)'
@inline Rmat(q::Quaternion{T})  where T = SMatrix{4,4,T,16}(q[1],q[2],q[3],q[4], -q[2], q[1],-q[4],q[3], -q[3], q[4],q[1],-q[2], -q[4], -q[3],q[2],q[1])
@inline RTmat(q::Quaternion) = Rmat(q)'

@inline rotate(x::SVector, q::Quaternion) = q.w^2*x + 2*q.w*cross(q.v,x) + 2*(q.v'*x)*q.v - (q.v'*q.v)*x

@inline function angleAxis(q::Quaternion{T}) where T
    n = norm(Vmat(T)*q)
    if n == 0
        return zero(T), SVector{3,T}(0,0,0)
    else
        return 2*atan(n,q[1]), Vmat()*q./n
    end
end
