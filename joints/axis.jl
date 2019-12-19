mutable struct Axis{T,Nc} <: Joint{T,Nc}
    V12::SMatrix{2,3,T,6}
    link2id::Int64

    function Axis(link1::AbstractLink{T},link2::AbstractLink{T},axis::AbstractVector{T}) where T
        Nc = 2
        V12 = (@SMatrix [1 0 0; 0 1 0])*svd(skew(axis)).Vt
        ids = SVector(link1.id,link2.id)

        new{T,Nc}(V12,link2.id), ids...
    end
end

@inline g(J::Axis,link1::Link,link2::Link,dt,No) = J.V12*Vmat(LTmat(getq3(link1,dt))*getq3(link2,dt))

@inline function ∂g∂posa(J::Axis{T},link1::Link,link2::Link,No) where T
    if link2.id == J.link2id
        X = @SMatrix zeros(T,2,3)

        R = -J.V12*Vmat(VTmat(Rmat(link2.q[No])*RTmat(link1.q[No])))

        return [X R]
    else
        return ∂g∂posa(J)
    end
end

@inline function ∂g∂posb(J::Axis{T},link1::Link,link2::Link,No) where T
    if link2.id == J.link2id
        X = @SMatrix zeros(T,2,3)

        R = J.V12*Vmat(VTmat(LTmat(link1.q[No])*Lmat(link2.q[No])))

        return [X R]
    else
        return ∂g∂posb(J)
    end
end

@inline function ∂g∂vela(J::Axis{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == J.link2id
        V = @SMatrix zeros(T,2,3)

        Ω = dt^2/4*J.V12*Vmat(Rmat(ωbar(link2,dt))*Rmat(link2.q[No])*RTmat(link1.q[No])*Tmat(T))*derivωbar(link1,dt)

        return [V Ω]
    else
        return ∂g∂vela(J)
    end
end

@inline function ∂g∂velb(J::Axis{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == J.link2id
        V = @SMatrix zeros(T,2,3)

        Ω = dt^2/4*J.V12*Vmat(LTmat(ωbar(link1,dt))*LTmat(link1.q[No])*Lmat(link2.q[No]))*derivωbar(link2,dt)

        return [V Ω]
    else
        return ∂g∂velb(J)
    end
end


@inline g(J::Axis,link1::Origin,link2::Link,dt,No) = J.V12*Vmat(getq3(link2,dt))

@inline function ∂g∂posb(J::Axis{T},link1::Origin,link2::Link,No) where T
    if link2.id == J.link2id
        X = @SMatrix zeros(T,2,3)

        R = J.V12*Vmat(VTmat(Lmat(link2.q[No])))

        return [X R]
    else
        return ∂g∂posb(J)
    end
end

@inline function ∂g∂velb(J::Axis{T},link1::Origin,link2::Link,dt,No) where T
    if link2.id == J.link2id
        V = @SMatrix zeros(T,2,3)

        Ω = dt/2*J.V12*Vmat(Lmat(link2.q[No]))*derivωbar(link2,dt)

        return [V Ω]
    else
        return ∂g∂velb(J)
    end
end
