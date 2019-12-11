struct Graph{N}
    adjacency::Vector{SVector{N,Bool}}
    dfsgraph::Vector{SVector{N,Bool}}
    dfslist::SVector{N,Int64}
    parentlist::SVector{N,Int64}
    root::Int64

    sizes::SVector{N,Int64}

    function Graph(constraints::Vector{<:Constraint};offset::Int64=0,root::Int64=1)
        adjacency = adjacencyMatrix(constraints,offset=offset)
        dfsgraph, dfslist = dfs(adjacency,root)

        N = length(adjacency)
        parentlist = zeros(Int64,N)
        for i=1:N
            parentlist[i] = parent(dfsgraph,i)
        end

        # sizes = zeros(Int64,currentGlobalID())
        sizes = zeros(Int64,N)
        for constraint in constraints
            sizes[constraint.data.id] = length(constraint)
            links = getlinks(constraint)
            for link in links
                sizes[link.data.id] = length(link)
            end
        end

        new{N}(adjacency,dfsgraph,dfslist,parentlist,root,sizes)
    end
end


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
