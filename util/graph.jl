struct Graph{N,F,NpF}
    root::Int64
    adjacency::Vector{SVector{N,Bool}}
    dfsgraph::Vector{SVector{N,Bool}}
    dfslist::SVector{N,Int64}
    parentlist::SVector{N,Int64}
    idlist::SVector{N,Int64} # matches ids to "nodes" index, 0 is for root

    pattern::Vector{SVector{N,Int64}} # includes fillins with ids of fillins at positions
    fillinid::SVector{NpF,Int64} # matches ids to "fillins" index in robot struct

    loops::Vector{Vector{Int64}}

    function Graph(origin::Link,links::Vector{<:Link},constraints::Vector{<:Constraint};offset::Int64=0,root::Int64=1)
        adjacency = adjacencyMatrix(constraints,offset=offset)
        dfsgraph, dfslist, loops = dfs(adjacency,root)

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
        pat,F = pattern(dfsgraph,root,loops,offset=N)
        NpF = N+F
        fillinid = zeros(Int64,NpF)
        for i=N+1:NpF
            fillinid[i]=i-N
        end

        new{N,F,NpF}(root,adjacency,dfsgraph,dfslist,parentlist,idlist,pat,fillinid,loops)
    end
end

@inline Base.size(g::Graph{N,F}) where {N,F} = N,F
@inline Base.length(g::Graph{N}) where N = N


function adjacencyMatrix(constraints::Vector{<:Constraint};offset::Int64=0)
    A = zero(Bool)
    n = 1

    for (i,constraint) in enumerate(constraints)
        for linkid in linkids(constraint)
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
    loops = Vector{Vector{Int64}}(undef,0)
    index = N

    dfslist[index] = rootid
    visited[rootid] = true
    dfs!(adjacency,dfsgraph,rootid,dfslist,visited,N,0,loops)
    loops = loops[sortperm(sort.(loops))[1:2:length(loops)]] # removes double entries of loop connections and keeps the first found pair

    dfsgraphvec = repeat([@SVector zeros(Bool,N)],N)
    for i=1:N
        dfsgraphvec[i] = convert(SVector{N,Bool},dfsgraph[i,:])
    end
    return dfsgraphvec, convert(SVector{N},dfslist), loops
end

function dfs!(A::Vector{SVector{N,T}},Adfs::Matrix,nodeid::Int64,list::Vector,visited::Vector,index::Int64,parentid::Int64,loops::Vector{Vector{Int64}}) where {N,T}
    for j=1:N
        if A[nodeid][j] && parentid!=j
            if visited[j]
                push!(loops,[j,nodeid])
            else
                index-=1
                visited[j] = true
                list[index] = j
                Adfs[nodeid,j] = true
                index = dfs!(A,Adfs,j,list,visited,index,nodeid,loops)
            end
        end
    end
    return index
end

@inline children(graph::Graph,n) = graph.dfsgraph[n]
@inline enumchildren(graph::Graph,n) = enumerate(graph.dfsgraph[n])
@inline enumpattern(graph::Graph,n) = enumerate(graph.pattern[n])
@inline parent(graph::Graph,n) = graph.parentlist[n]
@inline connected(graph::Graph,n) = graph.adjacency[n]
@inline enumconnected(graph::Graph,n) = enumerate(graph.adjacency[n])
@inline isroot(graph::Graph,n) = n==graph.root ? true : false
function ispredecessor(graph::Graph{N},p,c) where N # true if p is predecessor of c
    while c!=0
        c = parent(graph,c)
        p == c && (return true)
    end
    return false
end
issuccessor(graph::Graph,c,p) = ispredecessor(graph,p,c) # true if c is successor of p

function parent(dfsgraph::Vector{SVector{N,T}},n) where {N,T}
    for i=1:N
        dfsgraph[i][n] && (return i)
    end
    return 0
end

#TODO include loops
function pattern(dfsgraph::Vector{SVector{N,T}},rootid,loops;offset::Int64=0) where {N,T}
    pat = zeros(Int64,N,N)
    id = 0
    for i=1:N
        if i!=rootid
            for j=1:N
                j!=rootid && dfsgraph[i][j]==true && (id+=1;pat[i,j]=id+offset)
            end
        end
    end

    for loop in loops # loop = [index of starting constraint, index of last link in loop]
        pid = loop[1]
        cid = loop[2]
        while parent(dfsgraph,cid)!=loop[1]
            id+=1
            pat[pid,cid]=id+offset
            cid = parent(dfsgraph,cid)
        end
    end


    patternvec = repeat([@SVector zeros(Int64,N)],N)
    for i=1:N
        patternvec[i] = convert(SVector{N,Int64},pat[i,:])
    end
    return patternvec, id
end

function createfillins(graph::Graph,origin::Link{T},nodes::Vector) where T
    idlist = graph.idlist
    fillins = Vector{FillIn}(undef,0)
    for (i,row) in enumerate(graph.pattern)
        for (j,id) in enumerate(row)
            if id!=0
                parentnode = nodes[idlist[i]]
                childnode = nodes[idlist[j]]


                push!(fillins,FillIn{T,length(childnode),length(parentnode)}(id,childnode.data.id,parentnode.data.id))
            end
        end
    end
    return fillins
end
