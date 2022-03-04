Rotations.params(q::Quaternion) = SVector(q.s,q.v1,q.v2,q.v3)

Base.:*(q1::QuatRotation, q2::Quaternion) = q1.q*q2
Base.:*(q1::Quaternion, q2::QuatRotation) = q1*q2.q
Base.:/(q1::QuatRotation, q2::Quaternion) = q1.q/q2
Base.:/(q1::Quaternion, q2::QuatRotation) = q1/q2.q
Base.:\(q1::QuatRotation, q2::Quaternion) = q1.q\q2
Base.:\(q1::Quaternion, q2::QuatRotation) = q1\q2.q

Base.:*(q1::QuatRotation, q2::Real) = QuatRotation(q1.q*q2)
Base.:*(q1::Real, q2::QuatRotation) = QuatRotation(q1*q2.q)
Base.:/(q1::QuatRotation, q2::Real) = QuatRotation(q1.q/q2)
Base.:/(q1::Real, q2::QuatRotation) = QuatRotation(q1/q2.q)
Base.:\(q1::QuatRotation, q2::Real) = QuatRotation(q1.q\q2)
Base.:\(q1::Real, q2::QuatRotation) = QuatRotation(q1\q2.q)

Rotations.QuatRotation(w::T, v::StaticVector{3,T}, normalize::Bool = true) where T = QuatRotation{T}(w, v[1], v[2], v[3], normalize)
Rotations.QuatRotation(w::T, v::Vector{T}, normalize::Bool = true) where T = (@assert length(v)==3; QuatRotation{T}(w, v[1], v[2], v[3], normalize))
Rotations.QuatRotation(v::StaticVector{3,T}) where T = pure_quaternion(v)
Rotations.QuatRotation(v::Vector) = (@assert length(v)==3; pure_quaternion(v))
Rotations.QuatRotation(q::Quaternion) = QuatRotation(q.s,q.v1,q.v2,q.v3,false)

@inline imag(q::QuatRotation) = Rotations.vector(q)
@inline imag(q::Quaternion) = SVector(q.v1,q.v2,q.v3)

qrotate(q1::QuatRotation,q2::QuatRotation) = q2 * q1 / q2
qrotate(q1::QuatRotation,q2::Quaternion) = q2 * q1.q / q2
qrotate(q1::Quaternion,q2::QuatRotation) = q2.q * q1 / q2.q
qrotate(q1::Quaternion,q2::Quaternion) = q2 * q1 / q2
vrotate(v::Vector,q::QuatRotation) = imag(qrotate(pure_quaternion(v), q))
vrotate(v::StaticVector,q::QuatRotation) = q*v
vrotate(v::Vector,q::Quaternion) = imag(qrotate(pure_quaternion(v), q))
vrotate(v::StaticVector,q::Quaternion) = imag(qrotate(pure_quaternion(v), q))

rotation_vector(q::QuatRotation) = rotation_angle(q) * rotation_axis(q)

Lmat(q) = lmult(q)
Lᵀmat(q) = lmult(q)'
Rmat(q) = rmult(q)
Rᵀmat(q) = rmult(q)'

Tmat(::Type{T}=Float64) where T = tmat(T)
Tᵀmat(::Type{T}=Float64) where T = tmat(T)
Vmat(::Type{T}=Float64) where T = vmat(T)
Vᵀmat(::Type{T}=Float64) where T = hmat(T)
Vmat(q::QuatRotation) = imag(q)

function VLmat(q::QuatRotation)
    SA[
        q.x  q.w -q.z  q.y;
        q.y  q.z  q.w -q.x;
        q.z -q.y  q.x  q.w;
    ]
end
function VLᵀmat(q::QuatRotation)
    SA[
        -q.x  q.w  q.z -q.y;
        -q.y -q.z  q.w  q.x;
        -q.z  q.y -q.x  q.w;
    ]
end
function VRmat(q::QuatRotation)
    SA[
        q.x  q.w  q.z -q.y;
        q.y -q.z  q.w  q.x;
        q.z  q.y -q.x  q.w;
    ]
end
function VRᵀmat(q::QuatRotation)
    SA[
        -q.x  q.w -q.z  q.y;
        -q.y  q.z  q.w -q.x;
        -q.z -q.y  q.x  q.w;
    ]
end

function LVᵀmat(q::QuatRotation)
    SA[
        -q.x -q.y -q.z;
         q.w -q.z  q.y;
         q.z  q.w -q.x;
        -q.y  q.x  q.w;
    ]
end
function LᵀVᵀmat(q::QuatRotation)
    SA[
         q.x  q.y  q.z;
         q.w  q.z -q.y;
        -q.z  q.w  q.x;
         q.y -q.x  q.w;
    ]
end
function RVᵀmat(q::QuatRotation)
    SA[
        -q.x -q.y -q.z;
         q.w  q.z -q.y;
        -q.z  q.w  q.x;
         q.y -q.x  q.w;
    ]
end
function RᵀVᵀmat(q::QuatRotation)
    SA[
         q.x  q.y  q.z;
         q.w -q.z  q.y;
         q.z  q.w -q.x;
        -q.y  q.x  q.w;
    ]
end

function ∂L∂qsplit(::Type{T}) where T
    return SA{T}[
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1

        0 -1 0 0
        1 0 0 0
        0 0 0 1
        0 0 -1 0

        0 0 -1 0
        0 0 0 -1
        1 0 0 0
        0 1 0 0
        
        0 0 0 -1
        0 0 1 0
        0 -1 0 0
        1 0 0 0
    ]
end
function ∂Lᵀ∂qsplit(::Type{T}) where T
    return SA{T}[
        1 0 0 0
        0 -1 0 0
        0 0 -1 0
        0 0 0 -1
        
        0 1 0 0
        1 0 0 0
        0 0 0 -1
        0 0 1 0
        
        0 0 1 0
        0 0 0 1
        1 0 0 0
        0 -1 0 0
        
        0 0 0 1
        0 0 -1 0
        0 1 0 0
        1 0 0 0
    ]
end
function ∂R∂qsplit(::Type{T}) where T
    return SA{T}[
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1
        
        0 -1 0 0
        1 0 0 0
        0 0 0 -1
        0 0 1 0
        
        0 0 -1 0
        0 0 0 1
        1 0 0 0
        0 -1 0 0
        
        0 0 0 -1
        0 0 -1 0
        0 1 0 0
        1 0 0 0
    ]
end
function ∂Rᵀ∂qsplit(::Type{T}) where T
    return SA{T}[
        1 0 0 0
        0 -1 0 0
        0 0 -1 0
        0 0 0 -1
        
        0 1 0 0
        1 0 0 0
        0 0 0 1
        0 0 -1 0
        
        0 0 1 0
        0 0 0 -1
        1 0 0 0
        0 1 0 0
        
        0 0 0 1
        0 0 1 0
        0 -1 0 0
        1 0 0 0
    ]
end

function slerp(q1,q2,h)
    s = params(q1)'*params(q2)
    if s < 0
        s = -s
        q2 = -q2
    end

    qdiff = q1\q2
    φdiff = rotation_angle(qdiff) 
    udiff = rotation_axis(qdiff)
    φint = φdiff*h
    qint = QuatRotation(cos(φint/2),udiff*sin(φint/2),false)
    
    return q1*qint
end