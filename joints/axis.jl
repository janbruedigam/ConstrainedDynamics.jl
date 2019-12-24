mutable struct Axis{T,Nc} <: Joint{T,Nc}
    V12::SMatrix{2,3,T,6}
    cid::Int64

    function Axis(link1::AbstractLink{T},link2::AbstractLink{T},axis::AbstractVector{T}) where T
        Nc = 2
        V12 = (@SMatrix [1 0 0; 0 1 0])*svd(skew(axis)).Vt
        cid = link2.id

        new{T,Nc}(V12,cid), link1.id, link2.id
    end
end

@inline g(joint::Axis,link1::Link,link2::Link,dt,No) = joint.V12*VLᵀmat(getq3(link1,dt))*getq3(link2,dt)

@inline function ∂g∂posa(joint::Axis{T},link1::Link,link2::Link,No) where T
    if link2.id == joint.cid
        X = @SMatrix zeros(T,2,3)

        R = -joint.V12*VRmat(link2.q[No])*RᵀVᵀmat(link1.q[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Axis{T},link1::Link,link2::Link,No) where T
    if link2.id == joint.cid
        X = @SMatrix zeros(T,2,3)

        R = joint.V12*VLᵀmat(link1.q[No])*LVᵀmat(link2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Axis{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = @SMatrix zeros(T,2,3)

        Ω = dt^2/4*joint.V12*VRmat(ωbar(link2,dt))*Rmat(link2.q[No])*Rᵀmat(link1.q[No])*Tmat(T)*derivωbar(link1,dt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Axis{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = @SMatrix zeros(T,2,3)

        Ω = dt^2/4*joint.V12*VLᵀmat(ωbar(link1,dt))*Lᵀmat(link1.q[No])*Lmat(link2.q[No])*derivωbar(link2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline g(joint::Axis,link1::Origin,link2::Link,dt,No) = joint.V12*Vmat(getq3(link2,dt))

@inline function ∂g∂posb(joint::Axis{T},link1::Origin,link2::Link,No) where T
    if link2.id == joint.cid
        X = @SMatrix zeros(T,2,3)

        R = joint.V12*VLmat(link2.q[No])*Vᵀmat(T)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Axis{T},link1::Origin,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = @SMatrix zeros(T,2,3)

        Ω = dt/2*joint.V12*VLmat(link2.q[No])*derivωbar(link2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
