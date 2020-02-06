mutable struct OriginConnection{T,Nc} <: Joint{T,Nc}
    cid::Int64

    function OriginConnection(link1::Origin{T},link2::Link{T}) where T
        Nc = 0
        cid = link2.id

        new{T,Nc}(cid), link1.id, link2.id
    end
end

@inline g(joint::OriginConnection,link1::AbstractLink,link2::AbstractLink,dt,No) = g(joint)

@inline ∂g∂posa(joint::OriginConnection,link1::AbstractLink,link2::AbstractLink,No) = ∂g∂posa(joint)
@inline ∂g∂posb(joint::OriginConnection,link1::AbstractLink,link2::AbstractLink,No) = ∂g∂posb(joint)
@inline ∂g∂vela(joint::OriginConnection,link1::AbstractLink,link2::AbstractLink,dt,No) = ∂g∂vela(joint)
@inline ∂g∂velb(joint::OriginConnection,link1::AbstractLink,link2::AbstractLink,dt,No) = ∂g∂velb(joint)
