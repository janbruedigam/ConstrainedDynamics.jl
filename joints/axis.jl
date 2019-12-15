struct Axis{T,Nc} <: Joint{T,Nc}
    V12::SMatrix{2,3,T,6}

    function Axis(link1::AbstractLink{T},link2::AbstractLink{T},axis::AbstractVector{T}) where T
        Nc = 2
        V12 = (@SMatrix [1 0 0; 0 1 0])*svd(skew(axis)).Vt

        new{T,Nc}(V12), link1, link2
    end
end

@inline g(J::Axis,link1::Link,link2::Link) = J.V12*Vmat(LTmat(getq3(link1))*getq3(link2))
# @inline gŝ(J::Axis,link1::Link,link2::Link,data::NodeData) = (data.ŝ=J.V12*Vmat(LTmat(getq3(link1))*getq3(link2));nothing)
# @inline gf(J::Axis,link1::Link,link2::Link,data::NodeData) = (data.f=J.V12*Vmat(LTmat(getq3(link1))*getq3(link2));nothing)

@inline function ∂g∂posa(J::Axis{T},link1::Link,link2::Link) where T
    X = @SMatrix zeros(T,2,3)

    No = link1.No
    R = -J.V12*Vmat(VTmat(Rmat(link2.q[No])*RTmat(link1.q[No])))

    return [X R]
end

@inline function ∂g∂posb(J::Axis{T},link1::Link,link2::Link) where T
    X = @SMatrix zeros(T,2,3)

    No = link2.No
    R = J.V12*Vmat(VTmat(LTmat(link1.q[No])*Lmat(link2.q[No])))

    return [X R]
end

@inline function ∂g∂vela(J::Axis{T},link1::Link,link2::Link) where T
    V = @SMatrix zeros(T,2,3)

    No = link1.No
    Ω = link1.dt^2/4*J.V12*Vmat(Rmat(ωbar(link2))*Rmat(link2.q[No])*RTmat(link1.q[No])*Tmat(T))*derivωbar(link1)

    return [V Ω]
end

@inline function ∂g∂velb(J::Axis{T},link1::Link,link2::Link) where T
    V = @SMatrix zeros(T,2,3)

    No = link2.No
    Ω = link2.dt^2/4*J.V12*Vmat(LTmat(ωbar(link1))*LTmat(link1.q[No])*Lmat(link2.q[No]))*derivωbar(link2)

    return [V Ω]
end
