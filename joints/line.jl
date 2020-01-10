mutable struct Line{T,Nc} <: Joint{T,Nc}
    vertices::NTuple{2,SVector{3,T}}
    V12::SMatrix{2,3,T,6}
    cid::Int64

    function Line(link1::AbstractLink{T},link2::AbstractLink{T},axis::AbstractVector{T},p1::AbstractVector{T},p2::AbstractVector{T}) where T
        Nc = 2
        vertices = (p1,p2)
        V12 = (@SMatrix [1 0 0; 0 1 0])*svd(skew(axis)).Vt
        cid = link2.id

        new{T,Nc}(vertices,V12,cid), link1.id, link2.id
    end
end

@inline function g(joint::Line,link1::Link,link2::Link,dt,No)
    q1 = getq3(link1,dt)
    joint.V12*VLᵀmat(q1)*RVᵀmat(q1)*(getx3(link2,dt)-getx3(link1,dt))
end

@inline function ∂g∂posa(joint::Line{T},link1::Link,link2::Link,No) where T
    if link2.id == joint.cid
        q1 = link1.q[No]

        X = -joint.V12*VLᵀmat(q1)*RVᵀmat(q1)

        R = 2*joint.V12*VLᵀmat(q1)*Lmat(Quaternion(link2.x[No]-link1.x[No]))*LVᵀmat(q1)

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Line{T},link1::Link,link2::Link,No) where T
    if link2.id == joint.cid
        q1 = link1.q[No]

        X = joint.V12*VLᵀmat(q1)*RVᵀmat(q1)

        R = @SMatrix zeros(T,2,3)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Line{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == joint.cid
        q1 = link1.q[No]
        ω1 = ωbar(link1,dt)

        V = -dt^3/4*joint.V12*VLᵀmat(ω1)*Rmat(ω1)*Lᵀmat(q1)*RVᵀmat(q1)

        Ω = dt^2/4*joint.V12*VLᵀmat(ω1)*Lᵀmat(q1)*Lmat(Quaternion(link2.x[No]-link1.x[No]))*Lmat(q1)*derivωbar(link1,dt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Line{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == joint.cid
        q1 = link1.q[No]
        ω1 = ωbar(link1,dt)

        V = dt^3/4*joint.V12*VLᵀmat(ω1)*Rmat(ω1)*Lᵀmat(q1)*RVᵀmat(q1)

        Ω = @SMatrix zeros(T,2,3)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline g(joint::Line,link1::Origin,link2::Link,dt,No) = joint.V12*getx3(link2,dt)

@inline function ∂g∂posb(joint::Line{T},link1::Origin,link2::Link,No) where T
    if link2.id == joint.cid
        X = joint.V12

        R = @SMatrix zeros(T,2,3)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Line{T},link1::Origin,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = dt*joint.V12

        Ω = @SMatrix zeros(T,2,3)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
