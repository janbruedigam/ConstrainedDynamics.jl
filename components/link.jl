abstract type AbstractLink{T} <: Node{T} end

mutable struct Link{T} <: AbstractLink{T}
    id::Int64

    m::T
    J::SMatrix{3,3,T,9}

    x::Vector{SVector{3,T}}
    q::Vector{Quaternion{T}}

    F::Vector{SVector{3,T}}
    τ::Vector{SVector{3,T}}

    s0::SVector{6,T}
    s1::SVector{6,T}
    f::SVector{6,T}

    function Link(m::T,J::AbstractArray{T,2}) where T
        J = convert(SMatrix{3,3,T,9},J)

        x = [@SVector zeros(T,3)]
        q = [Quaternion{T}()]

        F = [@SVector zeros(T,3)]
        τ = [@SVector zeros(T,3)]

        s0 = @SVector zeros(T,6)
        s1 = @SVector zeros(T,6)
        f = @SVector zeros(T,6)

        new{T}(getGlobalID(),m,J,x,q,F,τ,s0,s1,f)
    end

    function Link(Box::Box)
        L = Link(Box.m,Box.J)
        push!(Box.linkids,L.id)
        return L
    end
end

Base.length(C::Link) = 6

struct Origin{T} <: AbstractLink{T}
    id::Int64

    Origin{T}() where T = new{T}(getGlobalID())
end

function setInit!(link::Link{T}; x::AbstractVector{T}=zeros(T,3), q::Quaternion{T}=Quaternion{T}(),
    F::AbstractVector{T}=zeros(T,3), τ::AbstractVector{T}=zeros(T,3)) where T

    link.x[1] = convert(SVector{3,T},x)
    link.q[1] = q
    link.F[1] = convert(SVector{3,T},F)
    link.τ[1] = convert(SVector{3,T},τ)

end

function setInit!(link1::Link{T}, link2::Link{T}, p1,p2,; q::Quaternion{T}=Quaternion{T}(),
    F::AbstractVector{T}=zeros(T,3), τ::AbstractVector{T}=zeros(T,3)) where T

    p1 = convert(SVector{3,T},p1)
    p2 = convert(SVector{3,T},p2)
    x2 = link1.x[1] + rotate(p1,link1.q[1]) - rotate(p2,q)

    setInit!(link2; x=x2, q=q, F=F, τ=τ)
end

function setInit!(link1::Origin{T}, link2::Link{T}, p1,p2,; q::Quaternion{T}=Quaternion{T}(),
    F::AbstractVector{T}=zeros(T,3), τ::AbstractVector{T}=zeros(T,3)) where T

    p1 = convert(SVector{3,T},p1)
    p2 = convert(SVector{3,T},p2)
    x2 = p1 - rotate(p2,q)

    setInit!(link2; x=x2, q=q, F=F, τ=τ)
end

getx3(link,dt) = getvnew(link)*dt + link.x[2]
getq3(link,dt) = Quaternion(dt/2*(Lmat(link.q[2])*ωbar(link,dt)))
getv1(link,dt) = (link.x[2]-link.x[1])/dt
function getω1(link,dt) # 2/link.dt*Vmat()*LTmat(link.q[1])*link.q[2]
    q1 = link.q[1]
    q2 = link.q[2]
    2/dt*(q1.w*q2.v-q2.w*q1.v-cross(q1.v,q2.v))
end

getvnew(link) = link.s1[SVector{3}(1:3)]
getωnew(link) = link.s1[SVector{3}(4:6)]

function derivωbar(link::Link{T},dt) where T
    ωnew = getωnew(link)
    msq = -sqrt(4/dt^2 - dot(ωnew,ωnew))
    [ωnew'/msq;SMatrix{3,3,T,9}(I)]
end
function ωbar(link,dt)
    ωnew = getωnew(link)
    Quaternion(sqrt(4/dt^2 - dot(ωnew,ωnew)),ωnew)
end


function dynamics(robot, link::Link{T}) where T
    No = robot.No
    dt = robot.dt
    ezg = SVector{3,T}(0,0,-robot.g)
    dynT = link.m*((getvnew(link) - getv1(link,dt))/dt + ezg) - link.F[No]

    J = link.J
    ω1 = getω1(link,dt)
    ωnew = getωnew(link)
    sq1 = sqrt(4/dt^2 - ω1'*ω1)
    sq2 = sqrt(4/dt^2 - ωnew'*ωnew)
    dynR = skewplusdiag(ωnew,sq2)*(J*ωnew) - skewplusdiag(ω1,sq1)*(J*ω1) - 2*link.τ[No]

    link.f = [dynT;dynR]

    # for cid in connections(robot.graph,link.id)
    #     # cid == -1 && break
    #     GtλTof!(robot,getconstraint(robot,cid),link)
    # end

    return link.f
end

function ∂dyn∂vel(robot, link::Link{T}) where T
    dt = robot.dt
    J = link.J
    ωnew = getωnew(link)
    sq = sqrt(4/dt^2 - ωnew'*ωnew)

    dynT = SMatrix{3,3,T,9}(link.m/dt*I)
    dynR = skewplusdiag(ωnew,sq)*J - J*ωnew*(ωnew'/sq) - skew(J*ωnew)

    Z = @SMatrix zeros(T,3,3)

    return [[dynT Z];[Z dynR]]
end
