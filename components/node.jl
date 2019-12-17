abstract type Node{T,N} end

#TODO do id differently?
CURRENTID = 1
getGlobalID() = (global CURRENTID+=1; return CURRENTID-1)
resetGlobalID() = (global CURRENTID=1; nothing)

Base.show(io::IO, N::Node) = summary(io, N)

Base.length(::Node{T,N}) where {T,N} = N
Base.foreach(f,itr::Vector{<:Node},arg...) = (for x in itr; f(x,arg...); end; nothing)

update!(node,diagonal) = (node.s1 = node.s0 - diagonal.ŝ; nothing)

s0tos1!(node) = (node.s1 = node.s0; nothing)
s1tos0!(node) = (node.s0 = node.s1; nothing)

function normΔs(node)
    diff = node.s1-node.s0
    return dot(diff,diff)
end
