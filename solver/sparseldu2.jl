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
    JL::SMatrix{N2,N1,T,N1N2}
    JU::SMatrix{N1,N2,T,N1N2}

    function OffDiagonalEntry{T,N1,N2}() where {T,N1,N2}
        JL = @SMatrix zeros(T,N2,N1)
        JU = @SMatrix zeros(T,N1,N2)
        new{T,N1,N2,N1*N2}(JL,JU)
    end
end

struct SparseLDU{T}
    diagonals::Vector{DiagonalEntry{T}}
    offdiagonals::Vector{OffDiagonalEntry{T}}

    ddict::Dict{Int64,Int64}
    odict::Dict{Tuple{Int64,Int64},Int64}

    function SparseLDU(graph::Graph{N},links::Vector{<:Link{T}},constraints::Vector{<:Constraint{T}},ldict::Dict,cdict::Dict) where {T,N}
        diagonals = Vector{DiagonalEntry{T}}(undef,0)
        ddict = Dict{Int64,Int64}()
        for link in links
            push!(diagonals,DiagonalEntry{T,length(link)}())
            ddict[link.id] = length(diagonals)
        end
        Nl = length(links)
        for constraint in constraints
            push!(diagonals,DiagonalEntry{T,length(constraint)}())
            ddict[constraint.id] = length(diagonals)
        end

        offdiagonals = Vector{OffDiagonalEntry{T}}(undef,0)
        odict = Dict{Tuple{Int64,Int64},Int64}()
        for id in graph.dfslist
            haskey(ldict,id) ? node=links[ldict[id]] : node=constraints[cdict[id]]
            N1 = length(node)

            for cid in successors(graph,id)
                haskey(ldict,cid) ? cnode=links[ldict[cid]] : cnode=constraints[cdict[cid]]
                N2 = length(cnode)

                push!(offdiagonals,OffDiagonalEntry{T,N2,N1}())
                odict[(id,cid)] = length(offdiagonals)
            end
        end

        new{T}(diagonals,offdiagonals,ddict,odict)
    end
end

@inline getentry(ldu::SparseLDU,id::Int64) = ldu.diagonals[ldu.ddict[id]]
@inline getentry(ldu::SparseLDU,ids::Tuple{Int64,Int64}) = ldu.offdiagonals[ldu.odict[ids]]


@inline function updateJ1!(offdiagonal::OffDiagonalEntry,d::DiagonalEntry,gc::OffDiagonalEntry,cgc::OffDiagonalEntry)
    offdiagonal.JL -= gc.JL*d.D*cgc.JU
    offdiagonal.JU -= cgc.JL*d.D*gc.JU
    return
end

@inline function updateJ2!(offdiagonal::OffDiagonalEntry,d::DiagonalEntry)
    offdiagonal.JL = offdiagonal.JL*d.Dinv
    offdiagonal.JU = d.Dinv*offdiagonal.JU
    return
end

@inline function updateD!(diagonal::DiagonalEntry,c::DiagonalEntry,f::OffDiagonalEntry)
    diagonal.D -= f.JL*c.D*f.JU
    return
end

invertD!(diagonal::DiagonalEntry) = (diagonal.Dinv = inv(diagonal.D); return)

@inline function LSol!(diagonal::DiagonalEntry,child::DiagonalEntry,fillin::OffDiagonalEntry)
    diagonal.ŝ -= fillin.JL*child.ŝ
    return
end

DSol!(diagonal::DiagonalEntry) = (diagonal.ŝ = diagonal.Dinv*diagonal.ŝ; return)

@inline function USol!(diagonal::DiagonalEntry,parent::DiagonalEntry,fillin::OffDiagonalEntry)
    diagonal.ŝ -= fillin.JU*parent.ŝ
    return
end


function factor!(graph::Graph,ldu::SparseLDU)
    for id in graph.dfslist
        # Keep this for now
        # for cid in successors(graph,id)
        #     offdiagonal = getentry(ldu,(id,cid))
        #     for gcid in successors(graph,cid)
        #         if hassuccessor(graph,id,gcid) # is actually a loop child
        #             updateJ1!(offdiagonal,getentry(ldu,gcid),getentry(ldu,(id,gcid)),getentry(ldu,(cid,gcid)))
        #         end
        #     end
        #     updateJ2!(offdiagonal,getentry(ldu,cid))
        # end
        sucs = successors(graph,id)
        for cid in sucs
            offdiagonal = getentry(ldu,(id,cid))
            for gcid in sucs
                gcid == cid && break
                if hasdirectchild(graph,cid,gcid)
                    updateJ1!(offdiagonal,getentry(ldu,gcid),getentry(ldu,(id,gcid)),getentry(ldu,(cid,gcid)))
                end
            end
            updateJ2!(offdiagonal,getentry(ldu,cid))
        end

        diagonal = getentry(ldu,id)

        for cid in successors(graph,id)
            updateD!(diagonal,getentry(ldu,cid),getentry(ldu,(id,cid)))
        end
        invertD!(diagonal)
    end
end

function solve!(graph::Graph,ldu::SparseLDU)
    dfslist = graph.dfslist

    for id in dfslist
        diagonal = getentry(ldu,id)

        for cid in successors(graph,id)
            LSol!(diagonal,getentry(ldu,cid),getentry(ldu,(id,cid)))
        end
    end

    for id in reverse(dfslist)
        diagonal = getentry(ldu,id)

        DSol!(diagonal)

        for pid in predecessors(graph,id)
            USol!(diagonal,getentry(ldu,pid),getentry(ldu,(pid,id)))
        end
    end
end

@inline update!(node::Node,ldu::SparseLDU) = update!(node,getentry(ldu,node.id))
