# abstract type AbstractLink{T} end

mutable struct Link{T,N} <: Node{T,N}
    id::Int64

    #TODO remove
    g::T
    dt::T
    No::Int64

    m::T
    J::SMatrix{3,3,T,9}

    x::Vector{SVector{3,T}}
    q::Vector{Quaternion{T}}

    F::Vector{SVector{3,T}}
    τ::Vector{SVector{3,T}}

    p::Vector{SVector{3,T}}

    s0::SVector{N,T}
    s1::SVector{N,T}
    f::SVector{N,T}

    function Link{N}(m::T,J::Array{T,2},p::Vector{<:AbstractVector{T}}) where {T,N}
        J = convert(SMatrix{3,3,T,9},J)

        x = [@SVector zeros(T,3)]
        q = [Quaternion{T}()]

        F = [@SVector zeros(T,3)]
        τ = [@SVector zeros(T,3)]

        p = convert(Vector{SVector{3,T}},p)

        s0 = @SVector zeros(T,N)
        s1 = @SVector zeros(T,N)
        f = @SVector zeros(T,N)

        g = 0
        dt = 0
        No = 0

        new{T,N}(getGlobalID(),g,dt,No,m,J,x,q,F,τ,p,s0,s1,f)
    end

    Link(p::Vector{<:AbstractVector{T}}) where T = Link{0}(zero(T),zeros(T,3,3),p)
    Link{T}() where T = Link([zeros(T,3)])
end

# struct Origin{T} <: AbstractLink{T}
#     id::Int64
#     p::Vector{SVector{3,T}}
#
#     Origin(p::Vector{<:AbstractVector{T}}) where T = new{T}(getGlobalID(),p)
#     Origin{T}() where T = Origin([zeros(T,3)])
# end



function setInit!(link::Link{T}; x::AbstractVector{T}=zeros(T,3), q::Quaternion{T}=Quaternion{T}(),
    F::AbstractVector{T}=zeros(T,3), τ::AbstractVector{T}=zeros(T,3)) where T

    link.x[1] = convert(SVector{3,T},x)
    link.q[1] = q
    link.F[1] = convert(SVector{3,T},F)
    link.τ[1] = convert(SVector{3,T},τ)

end

function setInit!(link1::Link{T}, link2::Link{T}, pids::Vector{Int64}; q::Quaternion{T}=Quaternion{T}(),
    F::AbstractVector{T}=zeros(T,3), τ::AbstractVector{T}=zeros(T,3)) where T

    p1 = link1.p[pids[1]]
    p2 = link2.p[pids[2]]
    x2 = link1.x[1] + rotate(p1,link1.q[1]) - rotate(p2,q)

    setInit!(link2; x=x2, q=q, F=F, τ=τ)
end

getv1(link::Link) = (link.x[2]-link.x[1])/link.dt
getvnew(link) = link.s1[SVector{3}(1:3)]
function getω1(link) # 2/link.dt*Vmat()*LTmat(link.q[1])*link.q[2]
    q1 = link.q[1]
    q2 = link.q[2]
    2/link.dt*(q1.w*q2.v-q2.w*q1.v-cross(q1.v,q2.v))
end
getωnew(link) = link.s1[SVector{3}(4:6)]
getx3(link) = getvnew(link)*link.dt + link.x[2]
getq3(link) = Quaternion(link.dt/2*(Lmat(link.q[2])*ωbar(link)))
derivωbar(link::Link{T}) where T = [-(getωnew(link)/(ωbar(link)[1]))';SMatrix{3,3,T,9}(I)]
ωbar(link) = Quaternion(sqrt(4/link.dt^2 - getωnew(link)'*getωnew(link)),getωnew(link))


getv1(link::Link{T,0}) where T = @SVector zeros(T,3)
getvnew(link::Link{T,0}) where T = @SVector zeros(T,3)
getω1(link::Link{T,0}) where T = @SVector zeros(T,3)
getx3(link::Link{T,0}) where T = @SVector zeros(T,3)
getq3(link::Link{T,0}) where T = Quaternion{T}()
derivωbar(link::Link{T,0}) where T = @SMatrix zeros(T,4,3)
ωbar(link::Link{T,0}) where T = Quaternion{T}()

function dynamics(link::Link{T}) where T
    No = link.No
    dt = link.dt
    ezg = SVector{3,T}(0,0,-link.g)
    dynT = link.m*((getvnew(link) - getv1(link))/dt + ezg) - link.F[No]

    J = link.J
    ω1 = getω1(link)
    ωnew = getωnew(link)
    sq1 = sqrt(4/dt^2 - ω1'*ω1)
    sq2 = sqrt(4/dt^2 - ωnew'*ωnew)
    dynR = skewplusdiag(ωnew,sq2)*(J*ωnew) - skewplusdiag(ω1,sq1)*(J*ω1) - 2*link.τ[No]

    return [dynT;dynR]
end

function ∂dyn∂vel(link::Link{T}) where T
    dt = link.dt
    J = link.J
    ωnew = getωnew(link)
    sq = sqrt(4/dt^2 - ωnew'*ωnew)

    dynT = SMatrix{3,3,T,9}(link.m/dt*I)
    dynR = skewplusdiag(ωnew,sq)*J - J*ωnew*(ωnew'/sq) - skew(J*ωnew)

    Z = @SMatrix zeros(T,3,3)

    return [[dynT Z];[Z dynR]]
end
