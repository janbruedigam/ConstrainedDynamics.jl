mutable struct Link{T,N,Nc,N²,NNc} <: Node{T,N,Nc,N²,NNc}
    No::Int64

    dt::T
    g::T
    m::T
    J::SMatrix{3,3,T,9}

    x::Vector{SVector{3,T}}
    q::Vector{Quaternion{T}}

    F::Vector{SVector{3,T}}
    τ::Vector{SVector{3,T}}

    trajectoryX::Vector{SVector{3,T}}
    trajectoryQ::Vector{Quaternion{T}}
    trajectoryΦ::Vector{T}

    p::Vector{SVector{3,T}}

    data::NodeData{T,N,Nc,N²,NNc}
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, L::Link)
    summary(io, L); println(io)
    print(io, "\nx_",L.No,": ")
    show(io, mime, L.x[L.No])
    print(io, "\nq_",L.No,": ")
    show(io, mime, L.q[L.No])
end

Base.show(io::IO, L::Link) = summary(io, L)

# function Base.getproperty(L::Link{T,N,Nc,N²,NNc,No},x::Symbol) where {T,N,Nc,N²,NNc,No}
#     if x == :No
#         return No
#     else
#         return getfield(L,x)
#     end
# end


function Link(m::T,J::Array{T,2},p::Array{Array{T,1},1},dof::Int64;No=2) where T
    Nc = 6-dof
    N = 6
    N² = N^2
    NNc = N*Nc
    Np = length(p)

    dt = 0
    g = 0
    J = convert(SMatrix{3,3,T,9},J)

    x = repeat([@SVector zeros(T,3)],No)
    q = repeat([Quaternion{T}()],No)

    F = repeat([@SVector zeros(T,3)],No)
    τ = repeat([@SVector zeros(T,3)],No)

    trajX = [@SVector zeros(T,3)]
    trajQ = [Quaternion{T}()]
    trajΦ = [0]

    p = convert(Vector{SVector{3,T}},p)

    data = NodeData{T,N,Nc}()

    Link{T,N,Nc,N²,NNc}(No,dt,g,m,J,x,q,F,τ,trajX,trajQ,trajΦ,p,data)
end
Link(T::Type;No=2) = Link(zero(T),diagm(zeros(T,3)),[zeros(T,3)],0;No=No)


function setInit!(link::Link{T}; x::AbstractVector{T}=zeros(T,3), q::Quaternion{T}=Quaternion{T}(),
    F::AbstractVector{T}=zeros(T,3), τ::AbstractVector{T}=zeros(T,3)) where T
    No = link.No
    for i=1:No
        link.x[i] = convert(SVector{3,T},x)
        link.q[i] = q
        link.F[i] = convert(SVector{3,T},F)
        link.τ[i] = convert(SVector{3,T},τ)
    end
end

function setInit!(link2::Link{T}, link1::Link{T}, pids::Vector{Int64}; q::Quaternion{T}=Quaternion{T}(),
    F::AbstractVector{T}=zeros(T,3), τ::AbstractVector{T}=zeros(T,3)) where T
    No = link1.No
    p1 = link1.p[pids[1]]
    p2 = link2.p[pids[2]]
    x2 = link1.x[No] + rotate(p1,link1.q[No]) - rotate(p2,q)

    setInit!(link2; x=x2, q=q, F=F, τ=τ)
end

# function initialPosition(link1::Link{T},link2::Link{T},q2::Quaternion{T},pId::Vector{Int64}) where T
#     No = link1.No
#     p1 = link1.p[pId[1]]
#     p2 = link2.p[pId[2]]
#
#     x2 = link1.x[No] + rotate(p1,link1.q[No]) - rotate(p2,q2)
# end

@inline getv1(link) = (link.x[2]-link.x[1])/link.dt
@inline getvnew(link) = link.data.s1[SVector{3}(1:3)]
@inline function getω1(link) # 2/link.dt*Vmat()*LTmat(link.q[1])*link.q[2]
    q1 = link.q[1]
    q2 = link.q[2]
    2/link.dt*(q1.w*q2.v-q2.w*q1.v-cross(q1.v,q2.v))
end
@inline getωnew(link) = link.data.s1[SVector{3}(4:6)]
@inline getx3(link) = getvnew(link)*link.dt + link.x[2]
@inline getq3(link) = Quaternion(link.dt/2*(Lmat(link.q[2])*ωbar(link)))
@inline derivωbar(link::Link{T}) where T = [-(getωnew(link)/(ωbar(link)[1]))';SMatrix{3,3,T,9}(I)]
@inline ωbar(link) = Quaternion(sqrt(4/link.dt^2 - getωnew(link)'*getωnew(link)),getωnew(link))

@inline function dynamics(link::Link{T}) where T
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

@inline function ∂dyn∂vel(link::Link{T}) where T
    dt = link.dt
    J = link.J
    ωnew = getωnew(link)
    sq = sqrt(4/dt^2 - ωnew'*ωnew)

    dynT = SMatrix{3,3,T,9}(link.m/dt*I)
    dynR = skewplusdiag(ωnew,sq)*J - J*ωnew*(ωnew'/sq) - skew(J*ωnew)

    return dynT, dynR
end
