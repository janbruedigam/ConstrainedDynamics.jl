# #TODO vectorize constraints and links
# struct Combined{T,N,Nc,Nl} <: Node{T,N}
#     id::Int64
#     linkids::SVector{Nl,Int64}
#
#     constr::Vector{Joint{T}}
#     links::Vector{Link{T}}
#     data::NodeData{T,N}
#
#     function Combined{T}(jointdata...) where T #c[i]=(constri,li1,li2)
#         constr = Vector{Joint{T}}(undef,0)
#         links = Vector{Link{T}}(undef,0)
#         N = 0
#         for set in jointdata
#             push!(constr,set[1])
#             push!(links,set[2])
#             push!(links,set[3])
#             N+=set[1].Nc
#         end
#         links = unique(links)
#
#         Nc = length(constr)
#         Nl = length(links)
#
#         id = getGlobalID()
#         linkids = [links[i].id for i=1:Nl]
#
#         data = NodeData{T,N}()
#
#         new{T,N,Nc,Nl}(id,linkids,constr,links,data)
#     end
# end
#
# # @inline g(C::Combined{T,N,Nc}) where {T,N,Nc} = [g(C.constr[i]) for i=1:Nc]
# @generated g(C::Combined{T,N,Nc}) where {T,N,Nc} = :(vcat($([:(g(C.constr[$i])) for i=1:Nc]...)))
#
#
# # @inline function ∂g∂pos(C::Combined2,L::Link)
# #     id = L.id
# #     ids = C.linkids
# #     if id == ids[1]
# #         return [∂g∂posa(C.constr1,L,C.link2);∂g∂posa(C.constr2,L,C.link2)]
# #     elseif id == ids[2]
# #         return [∂g∂posb(C.constr1,C.link1,L);∂g∂posb(C.constr2,C.link1,L)]
# #     else
# #         return [∂g∂posb(C.constr1);∂g∂posb(C.constr2)] #[0,0]
# #     end
# # end
# #
# # @inline function ∂g∂vel(C::Combined2,L::Link)
# #     id = L.id
# #     ids = C.linkids
# #     if id == ids[1]
# #         return [∂g∂vela(C.constr1,L,C.link2);∂g∂vela(C.constr2,L,C.link2)]
# #     elseif id == ids[2]
# #         return [∂g∂velb(C.constr1,C.link1,L);∂g∂velb(C.constr2,C.link1,L)]
# #     else
# #         return [∂g∂velb(C.constr1);∂g∂velb(C.constr2)] #[0,0]
# #     end
# # end
# #
# # getNc(C::Combined2) = C.constr1.Nc+C.constr2.Nc
# # # @generated function linkids(C::Combined2{T,Nc,Nc²,Nl}) where {T,Nc,Nc²,Nl}
# # #     ids = [:(C.links[$i].data.id) for i=1:Nl]
# # #     :(SVector{Nl,Int64}($(ids...)))
# # # end
