mutable struct MatchedOrientation{T,Nc} <: Joint{T,Nc}
    qoff::Quaternion{T} # link2's orientation in link1's frame
    cid::Int64

    function MatchedOrientation(link1::AbstractLink{T},link2::AbstractLink{T},qoff::Quaternion{T}) where T
        Nc = 3
        cid = link2.id

        new{T,Nc}(qoff,cid), link1.id, link2.id
    end
    function MatchedOrientation(link1::AbstractLink{T},link2::AbstractLink{T}) where T
        Nc = 3
        qoff = Quaternion{T}()
        cid = link2.id

        new{T,Nc}(qoff,cid), link1.id, link2.id
    end
end


@inline function g(joint::MatchedOrientation,link1::Link,link2::Link,dt,No)
    Vmat(getq3(link1,dt)*joint.qoff-getq3(link2,dt))
end

@inline function ∂g∂posa(joint::MatchedOrientation{T},link1::Link,link2::Link,No) where T
    if link2.id == joint.cid
        X = @SMatrix zeros(T,3,3)

        R = VRmat(joint.qoff)*LVᵀmat(link1.q[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::MatchedOrientation{T},link1::AbstractLink,link2::Link,No) where T
    if link2.id == joint.cid
        X = @SMatrix zeros(T,3,3)

        R = -Vmat(LVᵀmat(link2.q[No]))

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::MatchedOrientation{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = @SMatrix zeros(T,3,3)

        Ω = dt/2*VRmat(joint.qoff)*Lmat(link1.q[No])*derivωbar(link1,dt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::MatchedOrientation{T},link1::AbstractLink,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = @SMatrix zeros(T,3,3)

        Ω = -dt/2*VLmat(link2.q[No])*derivωbar(link2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end

@inline function g(joint::MatchedOrientation,link1::Origin,link2::Link,dt,No)
    Vmat(joint.qoff-getq3(link2,dt))
end
