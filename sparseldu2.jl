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

# struct SparseLDU{T}
#     diagonals::Vector{DiagonalEntry{T}}
#     offdiagonals::Vector{OffDiagonalEntry{T}}
#
#     ddict::Dict{Int64,Int64}
#     odict::Dict{Tuple{Int64,Int64},Int64}
#
#     function SparseLDU(graph::Graph{N},links::Vector{<:Link{T}},constraints::Vector{<:Constraint{T}},ldict::Dict,cdict::Dict) where {T,N}
#         diagonals = Vector{DiagonalEntry{T}}(undef,0)
#         ddict = Dict{Int64,Int64}()
#         for (ind,link) in enumerate(links)
#             push!(diagonals,DiagonalEntry{T,6}())
#             ddict[link.id] = length(diagonals)
#         end
#         Nl = length(links)
#         for (ind,constraint) in enumerate(constraints)
#             push!(diagonals,DiagonalEntry{T,constraint.Nc}())
#             ddict[constraint.id] = length(diagonals)
#         end
#
#         offdiagonals = Vector{OffDiagonalEntry{T}}(undef,0)
#         odict = Dict{Tuple{Int64,Int64},Int64}()
#         for (parentid,i) in pairs(graph.dict)
#             for (childid,j) in pairs(graph.dict)
#                 if graph.pattern[i][j]
#                     haskey(ldict,parentid) ? n1=links[ldict[parentid]] : n1=constraints[cdict[parentid]]
#                     haskey(ldict,childid) ? n2=links[ldict[childid]] : n2=constraints[cdict[childid]]
#
#                     typeof(n1)<:Link ? N1=6 : N1=n1.Nc
#                     typeof(n2)<:Link ? N2=6 : N2=n2.Nc
#                     push!(offdiagonals,OffDiagonalEntry{T,N1,N2}())
#                     odict[(parentid,childid)] = length(offdiagonals)
#                 end
#             end
#         end
#         new{T}(diagonals,offdiagonals,ddict,odict)
#     end
# end

# @inline getentry(ldu::SparseLDU,id::Int64) = ldu.diagonals[ldu.ddict[id]]
# @inline getentry(ldu::SparseLDU,ids::Tuple{Int64,Int64}) = ldu.offdiagonals[ldu.odict[ids]]
#
# @inline function setD!(entry::DiagonalEntry,link::Link,dt)
#     entry.data.D = ∂dyn∂vel(link,dt)
#     return nothing
# end
#
# @inline function setD!(entry::DiagonalEntry{T,Nc},::Constraint) where {T,Nc}
#     entry.D = @SMatrix zeros(T,Nc,Nc)
#     return nothing
# end
#
# @inline function setSol!(entry::DiagonalEntry,link::Link,dt,g,No)
#     entry.ŝ = dynamics(link,dt,g,No)
#     return nothing
# end
#
# @inline function setSol!(entry::DiagonalEntry,constraint::Constraint,dt)
#     entry.ŝ = g(constraint,dt)
#     return nothing
# end
#
# @inline function setJ!(entry::OffDiagonalEntry,link::Link,constraint::Constraint,dt,No)
#     entry.JL = ∂g∂vel(constraint,link,dt,No)
#     entry.JU = -∂g∂pos(constraint,link,No)'
#     return nothing
# end
#
# @inline function setJ!(entry::OffDiagonalEntry,constraint::Constraint,link::Link,dt,No)
#     entry.JL = -∂g∂pos(constraint,link,No)'
#     entry.JU = ∂g∂vel(constraint,link,dt,No)
#     return nothing
# end
#
# @inline function setJ!(entry::OffDiagonalEntry{T,N1,N2}) where {T,N1,N2}
#     entry.JL = @SMatrix zeros(T,N2,N1)
#     entry.JU = @SMatrix zeros(T,N1,N2)
#     return nothing
# end
#
# # function factor(ldu::SparseLDU, graph::Graph)
# #
# # end
# #
# # function solve(ldu::SparseLDU, graph::Graph)
# #
# # end
#
# function test(robot)
#     d = robot.ldu.diagonals[1]
#     setD!(d,robot.links[1],robot.dt)
#     nothing
# end
