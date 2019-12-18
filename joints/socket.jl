struct Socket{T,Nc} <: Joint{T,Nc}
    ps::Tuple{SVector{3,T},SVector{3,T}}
    link2id::Int64

    function Socket(link1::AbstractLink{T},link2::AbstractLink{T},p1::AbstractVector{T},p2::AbstractVector{T}) where T
        Nc = 3
        ps = (p1,p2)
        ids = SVector(link1.id,link2.id)

        new{T,Nc}(ps,link2.id), ids...
    end
end

function g(J::Socket,link1::Link,link2::Link)
    ps = J.ps
    getx3(link1) + rotate(ps[1],getq3(link1)) - (getx3(link2) + rotate(ps[2],getq3(link2)))
end

function ∂g∂posa(J::Socket{T},link1::Link,link2::Link) where T
    if link2.id == J.link2id
        X = SMatrix{3,3,T,9}(I)

        q = link1.q[link1.No]
        R = 2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(J.ps[1]))*Lmat(q)))

        return [X R]
    else
        return ∂g∂posa(J)
    end
end

function ∂g∂posb(J::Socket{T},link1::AbstractLink,link2::Link) where T
    if link2.id == J.link2id
        X = SMatrix{3,3,T,9}(-I)

        q = link2.q[link2.No]
        R = -2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(J.ps[2]))*Lmat(q)))

        return [X R]
    else
        return ∂g∂posb(J)
    end
end

function ∂g∂vela(J::Socket{T},link1::Link,link2::Link) where T
    if link2.id == J.link2id
        V = link1.dt*SMatrix{3,3,T,9}(I)

        q = link1.q[link1.No]
        Ω = 2*link1.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link1))*Rmat(Quaternion(J.ps[1])))*derivωbar(link1)

        return [V Ω]
    else
        return ∂g∂vela(J)
    end
end

function ∂g∂velb(J::Socket{T},link1::AbstractLink,link2::Link) where T
    if link2.id == J.link2id
        V = link2.dt*SMatrix{3,3,T,9}(-I)

        q = link2.q[link2.No]
        Ω = -2*link2.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link2))*Rmat(Quaternion(J.ps[2])))*derivωbar(link2)

        return [V Ω]
    else
        return ∂g∂velb(J)
    end
end


function g(J::Socket,link1::Origin,link2::Link)
    ps = J.ps
    ps[1] - (getx3(link2) + rotate(ps[2],getq3(link2)))
end

# function ∂g∂posa(J::Socket{T},link1::Origin,link2::Link) where T
#     if SVector{2,Int64}(link1.id,link2.id) == J.linkids
#         X = SMatrix{3,3,T,9}(I)
#
#         R = -2*skew(J.ps[1])
#
#         return [X R]
#     else
#         return ∂g∂posa(J)
#     end
# end
#
# function ∂g∂vela(J::Socket{T},link1::Origin,link2::Link) where T
#     if SVector{2,Int64}(link1.id,link2.id) == J.linkids
#         V = link2.dt*SMatrix{3,3,T,9}(I)
#
#         Ω = -2*link2.dt/2*skew(J.ps[1])
#
#         return [V Ω]
#     else
#         return ∂g∂vela(J)
#     end
# end
