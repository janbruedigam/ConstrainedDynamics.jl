struct Graph{N,F,NpF}
    root::Int64
    adjacency::Vector{SVector{N,Bool}}
    dfsgraph::Vector{SVector{N,Bool}}
    dfslist::SVector{N,Int64}
    parentlist::SVector{N,Int64}
    idlist::SVector{N,Int64} # matches ids to "nodes" index, 0 is for root

    pattern::Vector{SVector{N,Int64}} # includes fillins with ids of fillins at positions
    fillinid::SVector{NpF,Int64} # matches ids to "fillins" index in robot struct

    function Graph(origin::Link,links::Vector{<:Link},constraints::Vector{<:Constraint};offset::Int64=0,root::Int64=1)
        adjacency = adjacencyMatrix(constraints,offset=offset)
        dfsgraph, dfslist = dfs(adjacency,root)

        N = length(adjacency)
        parentlist = zeros(Int64,N)
        for i=1:N
            parentlist[i] = parent(dfsgraph,i)
        end

        idlist = zeros(Int64,N)
        idlist[origin.data.id] = 0
        Nl = length(links)
        for (i,link) in enumerate(links)
            idlist[link.data.id] = i
        end
        for (i,constraint) in enumerate(constraints)
            idlist[constraint.data.id] = i+Nl
        end

        # TODO do fillinid propertly without offset
        pat,F = pattern(dfsgraph,root,offset=N)
        NpF = N+F
        fillinid = zeros(Int64,NpF)
        for i=N+1:NpF
            fillinid[i]=i-N
        end

        new{N,F,NpF}(root,adjacency,dfsgraph,dfslist,parentlist,idlist,pat,fillinid)
    end
end

@inline Base.size(g::Graph{N,F}) where {N,F} = N,F
@inline Base.length(g::Graph{N}) where N = N


function adjacencyMatrix(constraints::Vector{<:Constraint};offset::Int64=0)
    A = zero(Bool)
    n = 1

    for (i,constraint) in enumerate(constraints)
        for linkid in constraint.linkids#linkids(constraint)
            cid = constraint.data.id
            m = maximum([linkid;cid])
            if m>n
                A = [A zeros(Bool,n,m-n); zeros(Bool,m-n,n) zeros(Bool,m-n,m-n)]
                n = m
            end
            A[cid,linkid] = true
        end
    end

    A=A.|A'
    Avec = repeat([@SVector zeros(Bool,n)],n)
    for i=1:n
        Avec[i] = convert(SVector{n,Bool},A[i,:])
    end

    return Avec
end

function dfs(adjacency::Vector{SVector{N,T}},rootid) where {N,T}
    dfsgraph = zeros(Bool,N,N)
    dfslist = zeros(Int64,N)
    visited = zeros(Bool,N)
    index = N

    dfslist[index] = rootid
    visited[rootid] = true
    dfs!(adjacency,dfsgraph,rootid,dfslist,visited,N)

    dfsgraphvec = repeat([@SVector zeros(Bool,N)],N)
    for i=1:N
        dfsgraphvec[i] = convert(SVector{N,Bool},dfsgraph[i,:])
    end
    return dfsgraphvec, convert(SVector{N},dfslist)
end

function dfs!(A::Vector{SVector{N,T}},Adfs::Matrix,nodeid::Int64,list::Vector,visited::Vector,index::Int64) where {N,T}
    for j=1:N
        if A[nodeid][j] && !visited[j]
            index-=1
            visited[j] = true
            list[index] = j
            Adfs[nodeid,j] = true
            index = dfs!(A,Adfs,j,list,visited,index)
        end
    end
    return index
end

@inline children(graph::Graph,n) = graph.dfsgraph[n]
@inline enumchildren(graph::Graph,n) = enumerate(graph.dfsgraph[n])
@inline parent(graph::Graph,n) = graph.parentlist[n]
@inline connected(graph::Graph,n) = graph.adjacency[n]
@inline enumconnected(graph::Graph,n) = enumerate(graph.adjacency[n])
@inline isroot(graph::Graph,n) = n==graph.root ? true : false

function parent(dfsgraph::Vector{SVector{N,T}},n) where {N,T}
    for i=1:N
        dfsgraph[i][n] && (return i)
    end
    return 0
end

#TODO include loops
function pattern(dfsgraph::Vector{SVector{N,T}},rootid;offset::Int64=0) where {N,T}
    pat = zeros(Int64,N,N)
    id = 0
    for i=1:N
        if i!=rootid
            for j=1:N
                j!=rootid && dfsgraph[i][j]==true && (id+=1;pat[i,j]=id+offset)
            end
        end
    end

    patternvec = repeat([@SVector zeros(Int64,N)],N)
    for i=1:N
        patternvec[i] = convert(SVector{N,Int64},pat[i,:])
    end
    return patternvec, id
end

# function createfillins(graph::Graph,origin::Link{T},nodes::Vector) where T
#     idlist = graph.idlist
#     fillins = Vector{FillIn}(undef,0)
#     for n in graph.dfslist
#         if !isroot(graph,n)
#             p = parent(graph,n)
#
#             node = nodes[idlist[n]]
#             isroot(graph,p) ? parentnode=origin : parentnode=nodes[idlist[p]]
#             push!(fillins,FillIn{T,length(node),length(parentnode)}(parentnode.data.id,node.data.id))
#         end
#     end
#     return fillins
# end

function createfillins(graph::Graph,origin::Link{T},nodes::Vector) where T
    idlist = graph.idlist
    fillins = Vector{FillIn}(undef,0)
    for (i,row) in enumerate(graph.pattern)
        for (j,id) in enumerate(row)
            if id!=0
                childnode = nodes[idlist[i]]
                parentnode = nodes[idlist[j]]

                push!(fillins,FillIn{T,length(childnode),length(parentnode)}(id,childnode.data.id,parentnode.data.id))
            end
        end
    end
    return fillins
end
