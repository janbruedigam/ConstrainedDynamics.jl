mutable struct Socket{T,Nc} <: Joint{T,Nc}
    ps::Tuple{SVector{3,T},SVector{3,T}}
    link2id::Int64

    function Socket(link1::AbstractLink{T},link2::AbstractLink{T},p1::AbstractVector{T},p2::AbstractVector{T}) where T
        Nc = 3
        ps = (p1,p2)
        ids = SVector(link1.id,link2.id)

        new{T,Nc}(ps,link2.id), ids...
    end
end

@inline function g(J::Socket,link1::Link,link2::Link,dt,No)
    ps = J.ps
    getx3(link1,dt) + rotate(ps[1],getq3(link1,dt)) - (getx3(link2,dt) + rotate(ps[2],getq3(link2,dt)))
end

@inline function ∂g∂posa(J::Socket{T},link1::Link,link2::Link,No) where T
    if link2.id == J.link2id
        X = SMatrix{3,3,T,9}(I)

        q = link1.q[No]
        R = 2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(J.ps[1]))*Lmat(q)))

        return [X R]
    else
        return ∂g∂posa(J)
    end
end

@inline function ∂g∂posb(J::Socket{T},link1::AbstractLink,link2::Link,No) where T
    if link2.id == J.link2id
        X = SMatrix{3,3,T,9}(-I)

        q = link2.q[No]
        R = -2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(J.ps[2]))*Lmat(q)))

        return [X R]
    else
        return ∂g∂posb(J)
    end
end

@inline function ∂g∂vela(J::Socket{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == J.link2id
        V = SMatrix{3,3,T,9}(dt*I)

        q = link1.q[No]
        Ω = 2*dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link1,dt))*Rmat(Quaternion(J.ps[1])))*derivωbar(link1,dt)

        return [V Ω]
    else
        return ∂g∂vela(J)
    end
end

@inline function ∂g∂velb(J::Socket{T},link1::AbstractLink,link2::Link,dt,No) where T
    if link2.id == J.link2id
        V = SMatrix{3,3,T,9}(-dt*I)

        q = link2.q[No]
        Ω = -2*dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link2,dt))*Rmat(Quaternion(J.ps[2])))*derivωbar(link2,dt)

        return [V Ω]
    else
        return ∂g∂velb(J)
    end
end


@inline function g(J::Socket,link1::Origin,link2::Link,dt,No)
    ps = J.ps
    ps[1] - (getx3(link2,dt) + rotate(ps[2],getq3(link2,dt)))
end
