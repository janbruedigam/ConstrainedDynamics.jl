abstract type Node{T,N} end

#TODO do id differently?
CURRENTID = 1
getGlobalID() = (global CURRENTID+=1; return CURRENTID-1)
resetGlobalID() = (global CURRENTID=1; nothing)

Base.show(io::IO, N::Node) = summary(io, N)

Base.length(::Node{T,N}) where {T,N} = N
Base.foreach(f,itr::Vector{<:Node},arg...) = (for x in itr; f(x,arg...); end; nothing)

function setentries!(robot,nodes::Vector{<:Node})
    graph = robot.graph
    ldu = robot.ldu

    for node in nodes
        id = node.id

        for cid in successors(graph,id)
            cid == -1 && break
            offdiagonal = getentry(ldu,(id,cid))
            cnode =  getnode(robot,cid)
            setJ!(robot,offdiagonal,node,cnode)
        end

        diagonal = getentry(ldu,id)
        setD!(diagonal,node)
        setSol!(diagonal,node,robot)
    end
end

update!(node,diagonal) = (node.s1 = node.s0 - diagonal.ŝ; return)

s0tos1!(node) = (node.s1 = node.s0; return)
s1tos0!(node) = (node.s0 = node.s1; return)

function normΔs(node::Node{T})::T where T
    diff = node.s1-node.s0
    return dot(diff,diff)
end
