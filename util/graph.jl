struct Graph{N}
    adjacency::Vector{SVector{N,Bool}}
    dfsgraph::Vector{SVector{N,Bool}}
    pattern::Vector{SVector{N,Bool}} # includes fillins

    dfslist::SVector{N,Int64}

    dict::Dict{Int64,Int64}

    function Graph(origin::Link,links::Vector{<:Link},constraints::Vector{<:Constraint})
        adjacency, dict = adjacencyMatrix(constraints,offset=offset)
        dfsgraph, dfslist, loops = dfs(adjacency,dict,origin.id)
        pat = pattern(dfsgraph,dict,loops)
        deleteat!(adjacency,dict[origin.id])
        deleteat!(dfsgraph,dict[origin.id])
        deleteat!(pat,dict[origin.id])
        dfslist = deleteat!(dfslist,dict[origin.id])
        pop!(dict,origin.id)

        N = length(dict)
        adjacency = convert(SVector{N,SVector{N,Bool}},adjacency)
        dfsgraph = convert(SVector{N,SVector{N,Bool}},dfsgraph)
        pat = convert(SVector{N,SVector{N,Bool}},pat)

        new{N}(adjacency,dfsgraph,pat,dfslist,dict)
    end
end

function adjacencyMatrix(constraints::Vector{<:Constraint},links::Vector{<:Link})
    A = zeros(Bool,0,0)
    dict = Dict{Int64,Int64}()
    n = 0

    for constraint in constraints
        cid = constraint.id
        dict[cid] = n+=1
        A = [A zeros(Bool,n,1); zeros(Bool,1,n) zero(Bool)]
        for linkid in linkids(constraint)
            if !haskey(dict,linkid)
                dict[linkid] = n+=1
                A = [A zeros(Bool,n,1); zeros(Bool,1,n) zero(Bool)]
            end
            A[dict[cid],dict[linkid]] = true
        end
    end
    for link in links # add unconnected links
        if !haskey(dict,link.id)
            dict[link.id] = n+=1
            A = [A zeros(Bool,n,1); zeros(Bool,1,n) zero(Bool)]
        end
    end

    A = A.|A'
    A, dict
end

function dfs(adjacency::Matrix,dict::Dict,rootid::Int64)
    N = size(adjacency)[1]
    dfsgraph = zeros(Bool,N,N)
    dfslist = zeros(Int64,N)
    visited = zeros(Bool,N)
    loops = Vector{Vector{Int64}}(undef,0)
    index = N

    dfslist[index] = rootid
    visited[dict[rootid]] = true
    dfs!(adjacency,dfsgraph,dict,rootid,dfslist,visited,N,0,loops)
    loops = loops[sortperm(sort.(loops))[1:2:length(loops)]] # removes double entries of loop connections and keeps the first found pair

    return dfsgraph, convert(SVector{N},dfslist), loops
end

function dfs!(A::Matrix,Adfs::Matrix,dict::Dict,list::Vector,visited::Vector,loops::Vector{Vector{Int64}},listindex::Int64,currentid::Int64,parentid::Int64) where {N,T}
    i = dict[currentid]
    for (j,childid) in pairs(dict)
        if A[i,j] && parentid != childid # connection from i to j in adjacency && not a direct connection back to the parent
            if visited[j]
                push!(loops,[childid,currentid]) # childid is actually a predecessor of currentid since it's a loop
            else
                index-=1
                list[index] = childid
                visited[j] = true
                Adfs[i,j] = true
                index = dfs!(A,Adfs,dict,list,visited,loops,index,childid,currentid)
            end
        end
    end
    return index
end

function pattern(dfsgraph::Matrix,dict::Dict,loops::Vector{Vector{Int64}})
    pat = dfsgraph

    for loop in loops # loop = [index of starting constraint, index of last link in loop]
        startid = loop[1]
        endid = loop[2]
        currentid = endid
        nextid = parent(dfsgraph,dict,currentid)
        while nextid != startid
            pat[parentid,current]=true
            currentid = nextid
            nextid = parent(dfsgraph,dict,next)
        end
    end

    return pat
end

function parent(dfsgraph::Vector{SVector{N,T}},dict::Dict,childid::Int64) where {N,T}
    j = dict[childid]
    for (i,parentid) in pairs(dict)
        dfsgraph[i][j] && (return parentid)
    end
    return -1
end
