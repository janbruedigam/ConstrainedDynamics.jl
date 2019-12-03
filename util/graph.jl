function adjacencyMat(constraints::Vector{<:Constraint},N::Int64)
    A = zeros(Bool,N,N)

    for constraint in constraints
        for linkid in linkids(constraint)
            A[constraint.data.id,linkid] = true
        end
    end

    A=A.|A'
    Avec = repeat([@SVector zeros(Bool,N)],N)
    for i=1:N
        Avec[i] = convert(SVector{N,Bool},A[i,:])
    end

    return Avec
end

function dfs(A::Vector{SVector{N,T}};root=1) where {N,T}
    Adfs = zeros(Bool,N,N)
    list = zeros(Int64,N)
    visited = zeros(Bool,N)
    index = N

    list[index] = root
    visited[root] = true
    dfs!(A,Adfs,root,list,visited,N)

    Adfsvec = repeat([@SVector zeros(Bool,N)],N)
    for i=1:N
        Adfsvec[i] = convert(SVector{N,Bool},Adfs[i,:])
    end
    return Adfsvec, convert(SVector{N},list)
end

function dfs!(A::Vector{SVector{N,T}},Adfs::Matrix,i::Int64,list::Vector,visited::Vector,index::Int64) where {N,T}
    for j=1:N
        if A[i][j] && !visited[j]
            index-=1
            visited[j] = true
            list[index] = j
            Adfs[i,j] = true
            index = dfs!(A,Adfs,j,list,visited,index)
        end
    end
    return index
end

function parentNode(dfsgraph::Vector{SVector{N,T}},n) where {N,T}
    for i=1:N
        dfsgraph[i][n] ? (return i) : nothing
    end
    return 0
end
