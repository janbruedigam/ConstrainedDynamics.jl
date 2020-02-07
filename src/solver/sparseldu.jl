abstract type Entry{T} end

mutable struct DiagonalEntry{T,N,N²} <: Entry{T}
    D::SMatrix{N,N,T,N²}
    Dinv::SMatrix{N,N,T,N²}
    ŝ::SVector{N,T}

    function DiagonalEntry{T,N}() where {T,N}
        N² = N^2
        D = @SMatrix zeros(T,N,N)
        Dinv = @SMatrix zeros(T,N,N)
        ŝ = @SVector zeros(T,N)

        new{T,N,N²}(D,Dinv,ŝ)
    end
end

mutable struct OffDiagonalEntry{T,N1,N2,N1N2} <: Entry{T}
    L::SMatrix{N2,N1,T,N1N2}
    U::SMatrix{N1,N2,T,N1N2}

    function OffDiagonalEntry{T,N1,N2}() where {T,N1,N2}
        L = @SMatrix zeros(T,N2,N1)
        U = @SMatrix zeros(T,N1,N2)
        new{T,N1,N2,N1*N2}(L,U)
    end
end

struct SparseLDU{T}
    diagonals::UnitDict{Base.OneTo{Int64},DiagonalEntry{T}}
    offdiagonals::Dict{Tuple{Int64,Int64},OffDiagonalEntry{T}}

    function SparseLDU(graph::Graph{N},bodies::Vector{Body{T}},constraints::Vector{<:Constraint{T}},ldict::Dict,cdict::Dict) where {T,N}
        diagonals = Vector{DiagonalEntry{T}}(undef,0)
        for body in bodies
            push!(diagonals,DiagonalEntry{T,length(body)}())
        end
        for constraint in constraints
            push!(diagonals,DiagonalEntry{T,length(constraint)}())
        end
        diagonals = UnitDict(diagonals)

        offdiagonals = Dict{Tuple{Int64,Int64},OffDiagonalEntry{T}}()
        for id in graph.dfslist
            haskey(ldict,id) ? node=bodies[ldict[id]] : node=constraints[cdict[id]]
            N1 = length(node)

            for cid in successors(graph,id)
                haskey(ldict,cid) ? cnode=bodies[ldict[cid]] : cnode=constraints[cdict[cid]]
                N2 = length(cnode)

                offdiagonals[(id,cid)] = OffDiagonalEntry{T,N2,N1}()
            end
        end

        new{T}(diagonals,offdiagonals)
    end
end

@inline getentry(ldu::SparseLDU,id::Int64) = ldu.diagonals[id]
@inline getentry(ldu::SparseLDU,ids::Tuple{Int64,Int64}) = ldu.offdiagonals[ids]
