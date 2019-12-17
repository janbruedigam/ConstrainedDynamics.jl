struct Socket{T,Nc} <: Joint{T,Nc}
    pids::SVector{2,Int64}

    function Socket(link1::Link{T},link2::Link{T},pid1::Int64,pid2::Int64) where T
        Nc = 3
        pids = SVector(pid1,pid2)

        new{T,Nc}(pids), link1, link2
    end
end

function g(J::Socket,link1::Link,link2::Link)
    pids = J.pids
    getx3(link1) + rotate(link1.p[pids[1]],getq3(link1)) - (getx3(link2) + rotate(link2.p[pids[2]],getq3(link2)))
end

# function gŝ(J::Socket,link1::Link,link2::Link,data::NodeData)
#     pids = J.pids
#     data.ŝ = getx3(link1) + rotate(link1.p[pids[1]],getq3(link1)) - (getx3(link2) + rotate(link2.p[pids[2]],getq3(link2)))
#     return nothing
# end
#
# function gf(J::Socket,link1::Link,link2::Link,data::NodeData)
#     pids = J.pids
#     data.f = getx3(link1) + rotate(link1.p[pids[1]],getq3(link1)) - (getx3(link2) + rotate(link2.p[pids[2]],getq3(link2)))
#     return nothing
# end

function ∂g∂posa(J::Socket{T},link1::Link,link2::Link) where T
    X = SMatrix{3,3,T,9}(I)

    q = link1.q[link1.No]
    R = 2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(link1.p[J.pids[1]]))*Lmat(q)))

    return [X R]
end

function ∂g∂posb(J::Socket{T},link1::Link,link2::Link) where T
    X = SMatrix{3,3,T,9}(-I)

    q = link2.q[link2.No]
    R = -2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(link2.p[J.pids[2]]))*Lmat(q)))

    return [X R]
end

function ∂g∂vela(J::Socket{T},link1::Link,link2::Link) where T
    V = link1.dt*SMatrix{3,3,T,9}(I)

    q = link1.q[link1.No]
    Ω = 2*link1.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link1))*Rmat(Quaternion(link1.p[J.pids[1]])))*derivωbar(link1)

    return [V Ω]
end

function ∂g∂velb(J::Socket{T},link1::Link,link2::Link) where T
    V = link2.dt*SMatrix{3,3,T,9}(-I)

    q = link2.q[link2.No]
    Ω = -2*link2.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link2))*Rmat(Quaternion(link2.p[J.pids[2]])))*derivωbar(link2)

    return [V Ω]
end
