abstract type Node{T} end

#TODO do id differently?
CURRENTID = 1
getGlobalID() = (global CURRENTID+=1; return CURRENTID-1)
resetGlobalID() = (global CURRENTID=1; return)

Base.show(io::IO, N::Node) = summary(io, N)

# Base.length(::Node{T,N}) where {T,N} = N
@inline Base.foreach(f,itr::Vector{<:Node},arg...) = (for x in itr; f(x,arg...); end; return)

update!(node,diagonal) = (node.s1 = node.s0 - diagonal.ŝ; return)

s0tos1!(node) = (node.s1 = node.s0; return)
s1tos0!(node) = (node.s0 = node.s1; return)

@inline function normΔs(node::Node)
    diff = node.s1-node.s0
    return dot(diff,diff)
end
