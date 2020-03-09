mutable struct Translational0{T,Nc} <: Joint{T,Nc}
    cid::Int64

    function Translational0(body1::Origin{T}, body2::Body{T}) where T
        Nc = 0
        cid = body2.id

        new{T,Nc}(cid), body1.id, body2.id
    end
end

@inline function minimalCoordinates(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No)
    body2.x[No]
end
@inline g(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = g(joint)

@inline ∂g∂posa(joint::Translational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posa(joint)
@inline ∂g∂posb(joint::Translational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posb(joint)
@inline ∂g∂vela(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂vela(joint)
@inline ∂g∂velb(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂velb(joint)
