# function linearize_velkp1(mechanism::Mechanism{T}) where T
#     bodies = mechanism.bodies
    
#     velT = [SMatrix{3,3,T,9}(I) SMatrix{3,3,T,9}(Δt * I) SMatrix{3,3,T,9}(Δt^2/body.m * I)]
#     n = 0
#     for body in bodies
#         n += length(body)
#     end

#     A = zeros(T,n,n)

#     for body in bodies
#         id = body.id
#         A[(id-1)*6+1:(id-1)*6+6,(id-1)*6+1:(id-1)*6+6] = ∂dyn∂velm1(body, mechanism.Δt)
#     end

#     return A
# end