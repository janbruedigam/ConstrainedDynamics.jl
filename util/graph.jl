function adjacencyMat(constraints::Vector{<:Constraint},N::Int64)
    A = zeros(Bool,N,N)

    for constraint in constraints
        for linkid in linkids(constraint)
            A[constraint.data.id+1,linkid+1] = true
        end
    end

    A=A.|A'
    Avec = repeat([@SVector zeros(Bool,N)],N)
    for i=1:N
        Avec[i] = convert(SVector{N,Bool},A[i,:])
    end

    return Avec
end

function dfs(A::Vector{SVector{N,T}},rootid) where {N,T}
    rootid+=1
    Adfs = zeros(Bool,N,N)
    list = zeros(Int64,N)
    visited = zeros(Bool,N)
    index = N

    list[index] = rootid
    visited[rootid] = true
    dfs!(A,Adfs,rootid,list,visited,N)

    Adfsvec = repeat([@SVector zeros(Bool,N)],N)
    for i=1:N
        Adfsvec[i] = convert(SVector{N,Bool},Adfs[i,:])
    end
    list.-=1
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

SVdeleteat(a::SVector{N,T},i) where {T,N} = convert(SVector{N-1,T},deleteat!(convert(Vector,a),i))
