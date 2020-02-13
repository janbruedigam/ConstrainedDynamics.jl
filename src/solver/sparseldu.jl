abstract type Entry{T} end

mutable struct DiagonalEntry{T,N,N²} <: Entry{T}
    D::SMatrix{N,N,T,N²}
    Dinv::SMatrix{N,N,T,N²}
    ŝ::SVector{N,T}
    # sl::Float64
    # ga::Float64

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

mutable struct InequalityEntry{T,N,N²} <: Entry{T}
    sl::Float64
    ga::Float64
    slf::Float64
    psi::Float64
    b::Vector{Float64}

    function InequalityEntry{T,N}() where {T,N}
        N² = N^2
        sl = 0
        ga = 0

        new{T,N,N²}(sl,ga,0,0,zeros(2))
    end
end

struct SparseLDU{T}
    diagonals::UnitDict{Base.OneTo{Int64},DiagonalEntry{T}}
    offdiagonals::Dict{Tuple{Int64,Int64},OffDiagonalEntry{T}}
    inequalities::UnitDict{UnitRange{Int64},InequalityEntry{T}}

    function SparseLDU(graph::Graph{N},bodies::Vector{Body{T}},eqconstraints::Vector{<:EqualityConstraint{T}},
        ineqconstraints::Vector{<:InequalityConstraint{T}},bdict::Dict,eqdict::Dict,ineqdict::Dict) where {T,N}

        diagonals = Vector{DiagonalEntry{T}}(undef,0)
        for body in bodies
            push!(diagonals,DiagonalEntry{T,length(body)}())
        end
        for constraint in eqconstraints
            push!(diagonals,DiagonalEntry{T,length(constraint)}())
        end
        diagonals = UnitDict(diagonals)

        offdiagonals = Dict{Tuple{Int64,Int64},OffDiagonalEntry{T}}()
        for id in graph.dfslist
            haskey(bdict,id) ? node=bodies[bdict[id]] : node=eqconstraints[eqdict[id]]
            N1 = length(node)

            for cid in successors(graph,id)
                haskey(bdict,cid) ? cnode=bodies[bdict[cid]] : cnode=eqconstraints[eqdict[cid]]
                N2 = length(cnode)

                offdiagonals[(id,cid)] = OffDiagonalEntry{T,N2,N1}()
            end
        end

        inequalities = Vector{InequalityEntry{T}}(undef,0)
        for ineq in ineqconstraints
            push!(inequalities,InequalityEntry{T,length(ineq)}())
        end
        if !isempty(inequalities)
            inequalities = UnitDict((ineqconstraints[1].id):(ineqconstraints[length(inequalities)].id),inequalities)
        else
            inequalities = UnitDict(0:0,inequalities)
        end

        new{T}(diagonals,offdiagonals,inequalities)
    end
end

@inline getentry(ldu::SparseLDU,id::Int64) = ldu.diagonals[id]
@inline getentry(ldu::SparseLDU,ids::Tuple{Int64,Int64}) = ldu.offdiagonals[ids]
@inline getineq(ldu::SparseLDU,id::Int64) = ldu.inequalities[id]
