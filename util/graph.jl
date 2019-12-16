struct Graph{N}
    adjacency::Vector{SVector{N,Bool}}
    dfsgraph::Vector{SVector{N,Bool}}
    pattern::Vector{SVector{N,Bool}} # includes fillins
    fillins::Vector{SVector{N,Bool}}

    dfslist::SVector{N,Int64}

    dict::Dict{Int64,Int64}
    rdict::Dict{Int64,Int64}

    function Graph(origin::Link,links::Vector{<:Link},constraints::Vector{<:Constraint})
        adjacency, dict = adjacencyMatrix(constraints,links)
        dfsgraph, dfslist, loops = dfs(adjacency,dict,origin.id)
        pat = pattern(dfsgraph,dict,loops)
        fillins = convert(Matrix{Bool},dfsgraph .⊻ pat) # xor so only fillins remain

        adjacency = deleteat(adjacency,dict[origin.id])
        dfsgraph = deleteat(dfsgraph,dict[origin.id])
        pat = deleteat(pat,dict[origin.id])
        fillins = deleteat(fillins,dict[origin.id])
        dfslist = StaticArrays.deleteat(dfslist,length(dfslist))

        for (id,ind) in dict
            ind>dict[origin.id] && (dict[id] = ind-1)
        end
        pop!(dict,origin.id)
        rdict = Dict(ind => id for (id, ind) in dict)

        N = length(dict)
        #TODO make properly to convert
        adjacency = convert(SVector{N,SVector{N,Bool}},adjacency)
        dfsgraph = convert(SVector{N,SVector{N,Bool}},dfsgraph)
        pat = convert(SVector{N,SVector{N,Bool}},pat)
        fillins = convert(SVector{N,SVector{N,Bool}},fillins)

        new{N}(adjacency,dfsgraph,pat,fillins,dfslist,dict,rdict)
    end
end

function adjacencyMatrix(constraints::Vector{<:Constraint},links::Vector{<:Link})
    A = zeros(Bool,0,0)
    dict = Dict{Int64,Int64}()
    n = 0

    for constraint in constraints
        cid = constraint.id
        A = [A zeros(Bool,n,1); zeros(Bool,1,n) zero(Bool)]
        dict[cid] = n+=1
        for linkid in constraint.linkids
            if !haskey(dict,linkid)
                A = [A zeros(Bool,n,1); zeros(Bool,1,n) zero(Bool)]
                dict[linkid] = n+=1
            end
            A[dict[cid],dict[linkid]] = true
        end
    end
    for link in links # add unconnected links
        if !haskey(dict,link.id)
            A = [A zeros(Bool,n,1); zeros(Bool,1,n) zero(Bool)]
            dict[link.id] = n+=1
        end
    end

    A = A.|A'
    convert(Matrix{Bool},A), dict
end

function dfs(adjacency::Matrix,dict::Dict,originid::Int64)
    N = size(adjacency)[1]
    dfsgraph = zeros(Bool,N,N)
    dfslist = zeros(Int64,N)
    visited = zeros(Bool,N)
    loops = Vector{Vector{Int64}}(undef,0)
    index = N

    dfslist[index] = originid
    visited[dict[originid]] = true
    dfs!(adjacency,dfsgraph,dict,dfslist,visited,loops,N,originid,-1)
    loops = loops[sortperm(sort.(loops))[1:2:length(loops)]] # removes double entries of loop connections and keeps the first found pair

    return dfsgraph, convert(SVector{N},dfslist), loops
end

function dfs!(A::Matrix,Adfs::Matrix,dict::Dict,list::Vector,visited::Vector,loops::Vector{Vector{Int64}},index::Int64,currentid::Int64,parentid::Int64) where {N,T}
    i = dict[currentid]
    for (childid,j) in dict
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
    pat = deepcopy(dfsgraph)

    for loop in loops # loop = [index of starting constraint, index of last link in loop]
        startid = loop[1]
        endid = loop[2]
        currentid = endid
        nextid = parent(dfsgraph,dict,currentid)
        while nextid != startid
            pat[dict[startid],dict[currentid]]=true
            currentid = nextid
            nextid = parent(dfsgraph,dict,nextid)
        end
    end

    return pat
end

function parent(dfsgraph::Matrix,dict::Dict,childid::Int64) where {N,T}
    j = dict[childid]
    for (parentid,i) in dict
        dfsgraph[i,j] && (return parentid)
    end
    return -1
end
