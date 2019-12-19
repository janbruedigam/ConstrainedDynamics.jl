# #TODO update to make a 2-link joint
#
# struct FixedOrientation{T,Nc,Nl,L} <: Joint{T,Nc,Nl}
#     link::L
#
#     function FixedOrientation(link::Link{T}) where T
#         Nc = 3
#         Nl = 1
#         FixedOrientation{T,Nc,Nl,typeof(link)}(link)
#     end
# end
#
#
# g(C::FixedOrientation{T}) where T = Vmat(getq3(C.link))
#
# function ∂g∂posa(C::FixedOrientation{T}) where T
#     link = C.link
#     X = @SMatrix zeros(T,3,3)
#
#     R = Vmat(VTmat(Lmat(link.q[link.No])))
#
#     return [X R]
# end
#
# function ∂g∂vela(C::FixedOrientation{T}) where T
#     link = C.link
#     V = @SMatrix zeros(T,3,3)
#
#     Ω = link.dt/2*Vmat(Lmat(link.q[link.No])*derivωbar(link))
#
#     return [V Ω]
# end
