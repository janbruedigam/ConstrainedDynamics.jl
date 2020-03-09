mutable struct Translational0{T,Nc} <: Joint{T,Nc}
    cid::Int64

    function Translational0(body1::Origin{T}, body2::Body{T}) where T
        Nc = 0
        cid = body2.id

        new{T,Nc}(cid), body1.id, body2.id
    end
end

function setForce!(joint::Translational0, body1::AbstractBody, body2::Body, F::SVector{0,T}, No) where T
    return
end

function setForce!(joint::Translational0, body1::AbstractBody, body2::Body, F::Nothing, No)
    return
end

function setForce!(joint::Translational0, body1::Body, body2::Body{T}, F::SVector{3,T}, No) where T
    F1 = vrotate(-F,body1.q[No])
    F2 = -F1

    body2.F[No] = F1
    body2.F[No] = F2
    return
end

function setForce!(joint::Translational0, body1::Origin, body2::Body{T}, F::SVector{3,T}, No) where T
    body2.F[No] = F
    return
end

@inline function minimalCoordinates(joint::Translational0, body1::AbstractBody, body2::AbstractBody, No)
    body2.x[No]
end
@inline g(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = g(joint)

@inline ∂g∂posa(joint::Translational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posa(joint)
@inline ∂g∂posb(joint::Translational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posb(joint)
@inline ∂g∂vela(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂vela(joint)
@inline ∂g∂velb(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂velb(joint)
