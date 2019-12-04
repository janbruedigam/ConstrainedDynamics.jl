abstract type Node{T,N1,N2,N1²,N1N2} end

mutable struct NodeData{T,N1,N2,N1²,N1N2}
    id::Int64
    s0::SVector{N1,T}
    s1::SVector{N1,T}
    ŝ::SVector{N1,T}
    f::SVector{N1,T}
    normf::T
    normDiff::T
    D::SMatrix{N1,N1,T,N1²}
    Dinv::SMatrix{N1,N1,T,N1²}
    JL::SMatrix{N2,N1,T,N1N2}
    JU::SMatrix{N1,N2,T,N1N2}
end

function NodeData{T,N1,N2}() where {T,N1,N2}
    N1² = N1^2
    N1N2 = N1*N2
    id = 0
    s0 = @SVector zeros(T,N1)
    s1 = @SVector zeros(T,N1)
    ŝ = @SVector zeros(T,N1)
    f = @SVector zeros(T,N1)
    normf = zero(T)
    normDiff = zero(T)
    D = @SMatrix zeros(T,N1,N1)
    Dinv = @SMatrix zeros(T,N1,N1)
    JL = @SMatrix zeros(T,N2,N1)
    JU = @SMatrix zeros(T,N1,N2)
    NodeData{T,N1,N2,N1²,N1N2}(id,s0,s1,ŝ,f,normf,normDiff,D,Dinv,JL,JU)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, N::NodeData)
    summary(io, N); println(io)
    print(io, "\nD: ")
    show(io, mime, N.D)
    print(io, "\nDinv: ")
    show(io, mime, N.Dinv)
    print(io, "\nJL: ")
    show(io, mime, N.JL)
    print(io, "\nJU: ")
    show(io, mime, N.JU)
end

Base.show(io::IO, N::Node) = summary(io, N)
@inline Base.foreach(f,itr::Vector{<:Node},arg) = (for x in itr; f(x,arg); end; nothing)

@inline update!(node) = (d = node.data; d.s1 = d.s0 - d.ŝ; nothing)

@inline s0tos1!(node) = (d = node.data; d.s1 = d.s0; nothing)
@inline s1tos0!(node) = (d = node.data; d.s0 = d.s1; nothing)

@inline function setNormDiff!(node)
    data = node.data
    diff = data.s1-data.s0
    data.normDiff = diff'*diff
    return nothing
end
