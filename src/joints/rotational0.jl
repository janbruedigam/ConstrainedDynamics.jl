mutable struct Rotational0{T,Nc} <: Joint{T,Nc}
    cid::Int64

    function Rotational0(body1::AbstractBody{T}, body2::Body{T}) where T
        Nc = 0
        cid = body2.id

        new{T,Nc}(cid), body1.id, body2.id
    end
end

function setForce!(joint::Rotational0, body1::AbstractBody, body2::Body, τ::SVector{0,T}, No) where T
    return
end

function setForce!(joint::Rotational0, body1::AbstractBody, body2::Body, τ::Nothing, No)
    return
end

function setForce!(joint::Rotational0, body1::Body, body2::Body{T}, τ::SVector{3,T}, No) where T
    τ1 = vrotate(-τ,body1.q[No])
    τ2 = -τ1

    body1.τ[No] = τ1
    body2.τ[No] = τ2
    return
end

function setForce!(joint::Rotational0, body1::Origin, body2::Body{T}, τ::SVector{3,T}, No) where T
    body2.τ[No] = τ
    return
end

@inline function minimalCoordinates(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, No)
    body2.q[No]
end
@inline g(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = g(joint)

@inline ∂g∂posa(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posa(joint)
@inline ∂g∂posb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posb(joint)
@inline ∂g∂vela(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂vela(joint)
@inline ∂g∂velb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂velb(joint)
