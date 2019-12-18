struct SocketYZ{T,Nc} <: Joint{T,Nc}
    ps::Tuple{SVector{3,T},SVector{3,T}}
    link2id::Int64

    function SocketYZ(link1::AbstractLink{T},link2::AbstractLink{T},p1::AbstractVector{T},p2::AbstractVector{T}) where T
        Nc = 2
        ps = (p1,p2)
        ids = SVector(link1.id,link2.id)

        new{T,Nc}(ps,link2.id), ids...
    end
end

function g(J::SocketYZ,link1::Link,link2::Link,dt,No)
    ps = J.ps
    (getx3(link1,dt) + rotate(ps[1],getq3(link1,dt)) - (getx3(link2,dt) + rotate(ps[2],getq3(link2,dt))))[SVector{2,Int64}(2,3)]
end

function ∂g∂posa(J::SocketYZ{T},link1::Link,link2::Link,No) where T
    if link2.id == J.link2id
        X = SMatrix{2,3,T,6}(0,0, 1,0, 0,1)

        q = link1.q[No]
        R = (2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(J.ps[1]))*Lmat(q))))[SVector{2,Int64}(2,3),:]

        return [X R]
    else
        return ∂g∂posa(J)
    end
end

function ∂g∂posb(J::SocketYZ{T},link1::AbstractLink,link2::Link,No) where T
    if link2.id == J.link2id
        X = SMatrix{2,3,T,6}(0,0, -1,0, 0,-1)

        q = link2.q[No]
        R = -(2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(J.ps[2]))*Lmat(q))))[SVector{2,Int64}(2,3),:]

        return [X R]
    else
        return ∂g∂posb(J)
    end
end

function ∂g∂vela(J::SocketYZ{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == J.link2id
        V = SMatrix{2,3,T,6}(0,0, dt,0, 0,dt)

        q = link1.q[No]
        Ω = (2*dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link1,dt))*Rmat(Quaternion(J.ps[1])))*derivωbar(link1,dt))[SVector{2,Int64}(2,3),:]

        return [V Ω]
    else
        return ∂g∂vel(J)
    end
end

function ∂g∂velb(J::SocketYZ{T},link1::AbstractLink,link2::Link,dt,No) where T
    if link2.id == J.link2id
        V = SMatrix{2,3,T,6}(0,0, -dt,0, 0,-dt)

        q = link2.q[No]
        Ω = -(2*dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link2,dt))*Rmat(Quaternion(J.ps[2])))*derivωbar(link2,dt))[SVector{2,Int64}(2,3),:]

        return [V Ω]
    else
        return ∂g∂velb(J)
    end
end


function g(J::SocketYZ,link1::Origin,link2::Link,dt,No)
    ps = J.ps
    (ps[1] - (getx3(link2,dt) + rotate(ps[2],getq3(link2,dt))))[SVector{2,Int64}(2,3)]
end
