abstract type Entry{T} end

mutable struct DiagonalEntry{T,N,N²} <: Entry{T}
    D::SMatrix{N,N,T,N²}
    Dinv::SMatrix{N,N,T,N²}
    Δs::SVector{N,T}

    function DiagonalEntry{T,N}() where {T,N}
        N² = N^2
        D = @SMatrix zeros(T, N, N)
        Dinv = @SMatrix zeros(T, N, N)
        Δs = @SVector zeros(T, N)

        new{T,N,N²}(D, Dinv, Δs)
    end
end

mutable struct OffDiagonalEntry{T,N1,N2,N1N2} <: Entry{T}
    L::SMatrix{N2,N1,T,N1N2}
    U::SMatrix{N1,N2,T,N1N2}

    function OffDiagonalEntry{T,N1,N2}() where {T,N1,N2}
        L = @SMatrix zeros(T, N2, N1)
        U = @SMatrix zeros(T, N1, N2)
        
        new{T,N1,N2,N1 * N2}(L, U)
    end
end

mutable struct InequalityEntry{T,N} <: Entry{T}
    Δs::SVector{N,T}
    Δγ::SVector{N,T}

    function InequalityEntry{T,N}() where {T,N}
        Δs = @SVector zeros(T, N)
        Δγ = @SVector zeros(T, N)

        new{T,N}(Δs, Δγ)
    end
end

struct SparseLDU{T}
    diagonals::UnitDict{Base.OneTo{Int64},DiagonalEntry{T}}
    offdiagonals::Dict{Tuple{Int64,Int64},OffDiagonalEntry{T}}
    inequalities::UnitDict{UnitRange{Int64},InequalityEntry{T}}

    function SparseLDU(graph::Graph{N},bodies::Vector{Body{T}},eqcs::Vector{<:EqualityConstraint{T}},
        ineqcs::Vector{<:InequalityConstraint{T}},bdict::Dict,eqdict::Dict,ineqdict::Dict) where {T,N}

        diagonals = DiagonalEntry{T}[]
        for body in bodies
            push!(diagonals, DiagonalEntry{T,length(body)}())
        end
        for constraint in eqcs
            push!(diagonals, DiagonalEntry{T,length(constraint)}())
        end
        diagonals = UnitDict(diagonals)

        offdiagonals = Dict{Tuple{Int64,Int64},OffDiagonalEntry{T}}()
        for id in graph.dfslist
            haskey(bdict, id) ? node = bodies[bdict[id]] : node = eqcs[eqdict[id]]
            N1 = length(node)

            for cid in successors(graph, id)
                haskey(bdict, cid) ? cnode = bodies[bdict[cid]] : cnode = eqcs[eqdict[cid]]
                N2 = length(cnode)

                offdiagonals[(id, cid)] = OffDiagonalEntry{T,N2,N1}()
            end
        end

        inequalities = InequalityEntry{T}[]
        for ineq in ineqcs
            push!(inequalities, InequalityEntry{T,length(ineq)}())
        end
        if !isempty(inequalities)
            inequalities = UnitDict((ineqcs[1].id):(ineqcs[length(inequalities)].id), inequalities)
        else
            inequalities = UnitDict(0:0, inequalities)
        end

        new{T}(diagonals, offdiagonals, inequalities)
    end
end

@inline getentry(ldu::SparseLDU, id::Integer) = ldu.diagonals[id]
@inline getentry(ldu::SparseLDU, ids::Tuple{Integer,Integer}) = ldu.offdiagonals[ids]
@inline getineqentry(ldu::SparseLDU, id::Integer) = ldu.inequalities[id]
