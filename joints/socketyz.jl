struct SocketYZ{T,Nc} <: Joint{T,Nc}
    pids::SVector{2,Int64}

    function SocketYZ(link1::Link{T},link2::Link{T},pid1::Int64,pid2::Int64) where T
        Nc = 2
        pids = SVector(pid1,pid2)

        new{T,Nc}(pids), link1, link2
    end
end

function g(J::SocketYZ,link1::Link,link2::Link)
    pids = J.pids
    (getx3(link1) + rotate(link1.p[pids[1]],getq3(link1)) - (getx3(link2) + rotate(link2.p[pids[2]],getq3(link2))))[SVector{2,Int64}(2,3)]
end

function ∂g∂posa(C::SocketYZ{T},link1::Link,link2::Link) where T
    X = SMatrix{2,3,T,6}(0,0, 1,0, 0,1)

    q = link1.q[link1.No]
    R = (2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(link1.p[C.pids[1]]))*Lmat(q))))[SVector{2,Int64}(2,3),:]

    return [X R]
end

function ∂g∂posb(C::SocketYZ{T},link1::Link,link2::Link) where T
    X = SMatrix{2,3,T,6}(0,0, -1,0, 0,-1)

    q = link2.q[link2.No]
    R = -(2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(link2.p[C.pids[2]]))*Lmat(q))))[SVector{2,Int64}(2,3),:]

    return [X R]
end

function ∂g∂vela(C::SocketYZ{T},link1::Link,link2::Link) where T
    V = link1.dt*SMatrix{2,3,T,6}(0,0, 1,0, 0,1)

    q = link1.q[link1.No]
    Ω = (2*link1.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link1))*Rmat(Quaternion(link1.p[C.pids[1]])))*derivωbar(link1))[SVector{2,Int64}(2,3),:]

    return [V Ω]
end

function ∂g∂velb(C::SocketYZ{T},link1::Link,link2::Link) where T
    V = link2.dt*SMatrix{2,3,T,6}(0,0, -1,0, 0,-1)

    q = link2.q[link2.No]
    Ω = -(2*link2.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link2))*Rmat(Quaternion(link2.p[C.pids[2]])))*derivωbar(link2))[SVector{2,Int64}(2,3),:]

    return [V Ω]
end
