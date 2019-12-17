struct Axis{T,Nc} <: Joint{T,Nc}
    V12::SMatrix{2,3,T,6}
    linkids::SVector{2,Int64}

    function Axis(link1::Link{T},link2::Link{T},axis::AbstractVector{T}) where T
        Nc = 2
        V12 = (@SMatrix [1 0 0; 0 1 0])*svd(skew(axis)).Vt
        ids = SVector{2,Int64}(link1.id,link2.id)

        new{T,Nc}(V12,ids), link1, link2
    end
end

g(J::Axis,link1::Link,link2::Link) = J.V12*Vmat(LTmat(getq3(link1))*getq3(link2))
# gŝ(J::Axis,link1::Link,link2::Link,data::NodeData) = (data.ŝ=J.V12*Vmat(LTmat(getq3(link1))*getq3(link2));nothing)
# gf(J::Axis,link1::Link,link2::Link,data::NodeData) = (data.f=J.V12*Vmat(LTmat(getq3(link1))*getq3(link2));nothing)

function ∂g∂posa(J::Axis{T},link1::Link,link2::Link) where T
    if SVector{2,Int64}(link1.id,link2.id) == J.linkids
        X = @SMatrix zeros(T,2,3)

        No = link1.No
        R = -J.V12*Vmat(VTmat(Rmat(link2.q[No])*RTmat(link1.q[No])))

        return [X R]
    else
        return ∂g∂posa(J)
    end
end

function ∂g∂posb(J::Axis{T},link1::Link,link2::Link) where T
    if SVector{2,Int64}(link1.id,link2.id) == J.linkids
        X = @SMatrix zeros(T,2,3)

        No = link2.No
        R = J.V12*Vmat(VTmat(LTmat(link1.q[No])*Lmat(link2.q[No])))

        return [X R]
    else
        return ∂g∂posb(J)
    end
end

function ∂g∂vela(J::Axis{T},link1::Link,link2::Link) where T
    if SVector{2,Int64}(link1.id,link2.id) == J.linkids
        V = @SMatrix zeros(T,2,3)

        No = link1.No
        Ω = link1.dt^2/4*J.V12*Vmat(Rmat(ωbar(link2))*Rmat(link2.q[No])*RTmat(link1.q[No])*Tmat(T))*derivωbar(link1)

        return [V Ω]
    else
        return ∂g∂vela(J)
    end
end

function ∂g∂velb(J::Axis{T},link1::Link,link2::Link) where T
    if SVector{2,Int64}(link1.id,link2.id) == J.linkids
        V = @SMatrix zeros(T,2,3)

        No = link2.No
        Ω = link2.dt^2/4*J.V12*Vmat(LTmat(ωbar(link1))*LTmat(link1.q[No])*Lmat(link2.q[No]))*derivωbar(link2)

        return [V Ω]
    else
        return ∂g∂velb(J)
    end
end
