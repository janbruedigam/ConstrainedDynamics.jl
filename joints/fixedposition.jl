# #TODO update to make a 2-link joint
#
# struct FixedPosition{T,Nc,Nl,L} <: Joint{T,Nc,Nl}
#     link::L
#     pid::Int64
#
#     function FixedPosition(link::Link{T},pid) where T
#         Nc = 3
#         Nl = 1
#         FixedPosition{T,Nc,Nl,typeof(link)}(link,pid)
#     end
# end
#
#
# function g(C::FixedPosition)
#     link = C.link
#     getx3(link) + rotate(link.p[C.pid],getq3(link))
# end
#
# function ∂g∂posa(C::FixedPosition{T}) where T
#     link = C.link
#     X = SMatrix{3,3,T,9}(I)
#
#     q = link.q[link.No]
#     R = 2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(link.p[C.pid]))*Lmat(q)))
#
#     return [X R]
# end
#
# function ∂g∂vela(C::FixedPosition{T}) where T
#     link = C.link
#     V = SMatrix{3,3,T,9}(link.dt*I)
#
#     q = link.q[link.No]
#     Ω = 2*link.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link))*Rmat(Quaternion(link.p[C.pid])))*derivωbar(link)
#
#     return [V Ω]
# end
