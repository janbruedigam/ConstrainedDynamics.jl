# Quaternions.Quaternion(w::T, v::StaticVector{3,T}) where T = Quaternion{T}(w, v[1], v[2], v[3])
# Quaternions.Quaternion(w::T, v::Vector{T}) where T = (@assert length(v)==3; Quaternion{T}(w, v[1], v[2], v[3]))
# Quaternions.Quaternion(v::StaticVector{3,T}) where T = Quaternion{T}(0, v[1], v[2], v[3])
# Quaternions.Quaternion(v::StaticVector{4,T}) where T = Quaternion{T}(v[1], v[2], v[3], v[4])
Quaternions.Quaternion(v::AbstractVector) = (@assert length(v)==3; Quaternion(0, v[1], v[2], v[3]))

RotX(θ) = Quaternion(cos(θ/2), sin(θ/2), 0, 0)
RotY(θ) = Quaternion(cos(θ/2), 0, sin(θ/2), 0)
RotZ(θ) = Quaternion(cos(θ/2), 0, 0, sin(θ/2))

quateltype(x) = eltype(x) # TODO not super elegant
quateltype(::Quaternion{T}) where T = T

Rotations.params(q::Quaternion{T}) where T = SVector{4,T}(q.s, q.v1, q.v2, q.v3)

@inline imag(q::Quaternion) = SA[q.v1, q.v2, q.v3]

qrotate(q1::Quaternion,q2::Quaternion) = q2 * q1 / q2
vrotate(v::AbstractVector,q::Quaternion) = imag(qrotate(Quaternion(v), q))
# vrotate(v::StaticVector,q::Quaternion) = q*v

rotation_vector(q::Quaternion) = rotation_angle(q) * rotation_axis(q)

rotation_angle(q::Quaternion) = 2*atan(norm(SVector{3,Float64}(q.v1,q.v2,q.v3)),q.s)
function rotation_axis(q::Quaternion)
    if q.s == 1
        return SVector{3,Float64}(1,0,0)
    else
        qv = SVector{3,Float64}(q.v1,q.v2,q.v3)
        return qv/norm(qv)
    end
end 

Lmat(q) = lmult(q)
Lᵀmat(q) = lmult(q)'
Rmat(q) = rmult(q)
Rᵀmat(q) = rmult(q)'

Tmat(::Type{T}=Float64) where T = tmat(T)
Tᵀmat(::Type{T}=Float64) where T = tmat(T)
Vmat(::Type{T}=Float64) where T = vmat(T)
Vᵀmat(::Type{T}=Float64) where T = hmat(T)
Vmat(q::Quaternion) = imag(q)

function VLmat(q::Quaternion)
    SA[
        q.v1  q.s -q.v3  q.v2;
        q.v2  q.v3  q.s -q.v1;
        q.v3 -q.v2  q.v1  q.s;
    ]
end
function VLᵀmat(q::Quaternion)
    SA[
        -q.v1  q.s  q.v3 -q.v2;
        -q.v2 -q.v3  q.s  q.v1;
        -q.v3  q.v2 -q.v1  q.s;
    ]
end
function VRmat(q::Quaternion)
    SA[
        q.v1  q.s  q.v3 -q.v2;
        q.v2 -q.v3  q.s  q.v1;
        q.v3  q.v2 -q.v1  q.s;
    ]
end
function VRᵀmat(q::Quaternion)
    SA[
        -q.v1  q.s -q.v3  q.v2;
        -q.v2  q.v3  q.s -q.v1;
        -q.v3 -q.v2  q.v1  q.s;
    ]
end

function LVᵀmat(q::Quaternion)
    SA[
        -q.v1 -q.v2 -q.v3;
         q.s -q.v3  q.v2;
         q.v3  q.s -q.v1;
        -q.v2  q.v1  q.s;
    ]
end
function LᵀVᵀmat(q::Quaternion)
    SA[
         q.v1  q.v2  q.v3;
         q.s  q.v3 -q.v2;
        -q.v3  q.s  q.v1;
         q.v2 -q.v1  q.s;
    ]
end
function RVᵀmat(q::Quaternion)
    SA[
        -q.v1 -q.v2 -q.v3;
         q.s  q.v3 -q.v2;
        -q.v3  q.s  q.v1;
         q.v2 -q.v1  q.s;
    ]
end
function RᵀVᵀmat(q::Quaternion)
    SA[
         q.v1  q.v2  q.v3;
         q.s -q.v3  q.v2;
         q.v3  q.s -q.v1;
        -q.v2  q.v1  q.s;
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
    qint = Quaternion(cos(φint/2),sin(φint/2)*udiff...)
    
    return q1*qint
end