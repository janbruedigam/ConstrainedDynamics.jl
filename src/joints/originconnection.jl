mutable struct OriginConnection{T,Nc} <: Joint{T,Nc}
    cid::Int64

    function OriginConnection(body1::Origin{T},body2::Body{T}) where T
        Nc = 0
        cid = body2.id

        new{T,Nc}(cid), body1.id, body2.id
    end
end

@inline g(joint::OriginConnection,body1::AbstractBody,body2::AbstractBody,dt,No) = g(joint)

@inline ∂g∂posa(joint::OriginConnection,body1::AbstractBody,body2::AbstractBody,No) = ∂g∂posa(joint)
@inline ∂g∂posb(joint::OriginConnection,body1::AbstractBody,body2::AbstractBody,No) = ∂g∂posb(joint)
@inline ∂g∂vela(joint::OriginConnection,body1::AbstractBody,body2::AbstractBody,dt,No) = ∂g∂vela(joint)
@inline ∂g∂velb(joint::OriginConnection,body1::AbstractBody,body2::AbstractBody,dt,No) = ∂g∂velb(joint)
