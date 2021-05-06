mutable struct Rotational{T,N} <: Joint{T,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    qoffset::UnitQuaternion{T} # in body1's frame

    spring::T
    damper::T

    Fτ::SVector{3,T}

    childid::Int64

    function Rotational{T,N}(body1::AbstractBody, body2::AbstractBody; 
            axis::AbstractVector = zeros(3), qoffset::UnitQuaternion = one(UnitQuaternion), spring, damper
        ) where {T,N}
        
        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]

        Fτ = zeros(T,3)

        childid = body2.id

        new{T,N}(V3, V12, qoffset, spring, damper, Fτ, childid), body1.id, body2.id
    end
end

Rotational0 = Rotational{T,0} where T
Rotational1 = Rotational{T,1} where T
Rotational2 = Rotational{T,2} where T
Rotational3 = Rotational{T,3} where T

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Rotational{T,N}) where {T,N}
    summary(io, constraint)
    println(io,"")
    println(io, " V3:      "*string(constraint.V3))
    println(io, " V12:     "*string(constraint.V12))
    println(io, " qoffset: "*string(constraint.qoffset))
end

### Constraints and derivatives
## Position level constraints (for dynamics)
@inline function g(joint::Rotational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return Vmat(qa \ qb / joint.qoffset)
end
@inline function g(joint::Rotational, xb::AbstractVector, qb::UnitQuaternion)
    return Vmat(qb / joint.qoffset)
end

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)

    return X, Q
end
@inline function ∂g∂posb(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)

    return X, Q
end
@inline function ∂g∂posb(joint::Rotational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset)

    return X, Q
end

## vec(G) Jacobian (also NOT accounting for quaternion specialness in the second derivative: ∂(∂ʳg∂posx)∂y)
@inline function ∂2g∂posaa(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(Vmat(T),VRᵀmat(joint.qoffset)*Rmat(qb)*Tmat(T))*∂L∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posab(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(VLᵀmat(qa)*Tmat(T),VRᵀmat(joint.qoffset))*∂R∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posba(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(VLᵀmat(qb),VRᵀmat(joint.qoffset))*∂Lᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(Vmat(T),VRᵀmat(joint.qoffset)*Lᵀmat(qa))*∂L∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Rotational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(Vmat(T),VRᵀmat(joint.qoffset))*∂L∂qsplit(T)

    return XX, XQ, QX, QQ
end


### Forcing
## Application of joint forces (for dynamics)
@inline function applyFτ!(joint::Rotational{T}, statea::State, stateb::State, clear::Bool) where T
    τ = joint.Fτ
    xa, qa = posargsk(statea)
    xb, qb = posargsk(stateb)    

    τa = vrotate(-τ, qa) # in world coordinates
    τb = -τa # in world coordinates

    τa = vrotate(τa,inv(qa)) # in local coordinates
    τb = vrotate(τb,inv(qb)) # in local coordinates

    statea.τk[end] += τa
    stateb.τk[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end
@inline function applyFτ!(joint::Rotational{T}, stateb::State, clear::Bool) where T
    τ = joint.Fτ
    xb, qb = posargsk(stateb)

    τa = -τ # in world coordinates
    τb = -τa # in world coordinates

    τb = vrotate(τb,inv(qb)) # in local coordinates

    stateb.τk[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end

## Forcing derivatives (for linearization)
# Control derivatives
@inline function ∂Fτ∂ua(joint::Rotational{T}, statea::State, stateb::State) where T
    BFa = (szeros(T, 3, 3))
    Bτa = -I

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::Rotational{T}, statea::State, stateb::State) where T
    xa, qa = posargsk(statea)
    xb, qb = posargsk(stateb)
    qbinvqa = qb\qa

    BFb = (szeros(T, 3, 3))
    Bτb = VLmat(qbinvqa) * RᵀVᵀmat(qbinvqa)

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Rotational{T}, stateb::State) where T
    xb, qb = posargsk(stateb)

    BFb = (szeros(T, 3, 3))
    Bτb = VLᵀmat(qb) * RVᵀmat(qb)

    return [BFb; Bτb]
end

# Position derivatives
@inline function ∂Fτ∂posa(joint::Rotational{T}, statea::State, stateb::State) where T
    xa, qa = posargsk(statea)
    xb, qb = posargsk(stateb)
    τ = joint.Fτ

    FaXa = szeros(T,3,3)
    FaQa = szeros(T,3,3)
    τaXa = szeros(T,3,3)
    τaQa = szeros(T,3,3)
    FbXa = szeros(T,3,3)
    FbQa = szeros(T,3,3)
    τbXa = szeros(T,3,3)
    τbQa = 2*VLᵀmat(qb)*Rmat(qb)*Rᵀmat(qa)*Rmat(UnitQuaternion(τ))*LVᵀmat(qa)

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end
@inline function ∂Fτ∂posb(joint::Rotational{T}, statea::State, stateb::State) where T
    xa, qa = posargsk(statea)
    xb, qb = posargsk(stateb)
    τ = joint.Fτ

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,3)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,3)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,3)
    τbXb = szeros(T,3,3)
    τbQb = 2*VRᵀmat(qb)*Rᵀmat(qa)*Rmat(UnitQuaternion(τ))*Rmat(qa)*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
@inline function ∂Fτ∂posb(joint::Rotational{T}, stateb::State) where T
    xb, qb = posargsk(stateb)
    τ = joint.Fτ

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,3)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,3)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,3)
    τbXb = szeros(T,3,3)
    τbQb = 2*VRᵀmat(qb)*Rmat(UnitQuaternion(τ))*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end


### Minimal coordinates
## Position and velocity offsets 
@inline function getPositionDelta(joint::Rotational, body1::AbstractBody, body2::Body, θ::SVector{N,T}) where {T,N}
    # axis angle representation
    θ = zerodimstaticadjoint(nullspacemat(joint)) * θ
    nθ = norm(θ)
    if nθ == 0
        q = one(UnitQuaternion{T})
    else
        q = UnitQuaternion(cos(nθ/2),(θ/nθ*sin(nθ/2))..., false)
    end
    
    Δq = q * joint.qoffset # in body1 frame
    return Δq
end
@inline function getVelocityDelta(joint::Rotational, body1::Body, body2::Body, ω::SVector)
    ω = zerodimstaticadjoint(nullspacemat(joint)) * ω
    Δω = vrotate(ω, inv(body2.state.qc)*body1.state.qc) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational, body1::Origin, body2::Body, ω::SVector)
    ω = zerodimstaticadjoint(nullspacemat(joint)) * ω
    Δω = vrotate(ω, inv(body2.state.qc)) # in body2 frame
    return Δω
end

## Minimal coordinate calculation
@inline function minimalCoordinates(joint::Rotational, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    # q = g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
    q = statea.qc \ stateb.qc / joint.qoffset
    return nullspacemat(joint) * rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational, body1::Origin, body2::Body)
    stateb = body2.state
    # q = g(joint, stateb.xc, stateb.qc)
    q = stateb.qc / joint.qoffset
    return nullspacemat(joint) * rotation_vector(q)
end
@inline function minimalVelocities(joint::Rotational, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return nullspacemat(joint) * (vrotate(stateb.ωc,statea.qc\stateb.qc) - statea.ωc) # in body1's frame
end
@inline function minimalVelocities(joint::Rotational, body1::Origin, body2::Body)
    stateb = body2.state
    return nullspacemat(joint) * vrotate(stateb.ωc,stateb.qc) # in body1's frame
end
