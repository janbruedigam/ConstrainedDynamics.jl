# function newton_ip!(mechanism::Mechanism,body::Body)
#     nv = 6
#     ns = 1
#     nλ = 0
#     nγ = 1
#
#     v = body.s0
#     # λ = zeros(nλ)
#     γ = 0.#zeros(nγ)
#     s = rand()#rand(ns)
#     # e = ones(ns)
#     μ = 1.
#     τ = 0.995
#     σ = 0.2
#
#     dt = mechanism.dt
#
#     for iter=1:100
#         φ = body.x[2][3]+dt*body.s1[3]
#
#         dddv = ∂dyn∂vel(body, dt)
#         Nx = [0;0;1.;0;0;0]'
#         Nv = dt*Nx
#         Γ = γ#diagm(γ)
#         S = s#diagm(s)
#         In = Matrix{Float64}(I,ns,ns)
#
#         Mat = [
#             dddv         zeros(nv,ns)  -Nx'
#             zeros(ns,nv) Γ              S
#             Nv          -In             zeros(nγ,nγ)
#         ]
#
#         vec = [
#             dynamics(body, mechanism) - Nx'γ
#             S*γ .- μ
#             φ
#         ]
#
#         sol = Mat\vec
#
#         Δv = sol[1:nv]
#         Δs = sol[nv+1]
#         Δγ = sol[nv+2]
#
#         αs = ones(ns)
#         αsmax = 1.
#         αγ = ones(nγ)
#         αγmax = 1.
#
#         for i=1:ns
#             if Δs[i] > 0
#                 αs[i] = minimum([1.;τ*s[i]/Δs[i]])
#             end
#             αsmax = minimum(αs)
#
#             if Δγ[i] > 0
#                 αγ[i] = minimum([1.;τ*γ[i]/Δγ[i]])
#             end
#             αγmax = minimum(αγ)
#         end
#
#         v -= αsmax.*Δv
#         s -= αsmax.*Δs
#         γ -= αγmax.*Δγ
#
#         body.s1 = v
#         s1tos0!(body)
#
#         norm1 = norm(vec[1:nv])
#         norm2 = norm(vec[nv+1])
#         norm4 = norm(vec[nv+2])
#
#         E = maximum([norm1;norm2;norm4])
#         μ = σ*μ
#
#         if E<μ
#             display(iter)
#             return
#         end
#     end
#     display("failed")
#     return
# end

function newton_ip!(mechanism::Mechanism,body::Body)
    nv = 6
    ns = 1
    nλ = 0
    nγ = 1

    v = body.s0
    γ = rand()
    s = rand()
    μ = 1.
    τ = 0.995
    σ = 0.2

    dt = mechanism.dt

    for iter=1:100
        φ = body.x[2][3]+dt*body.s1[3]

        Nx = [0;0;1.;0;0;0]'
        Nv = dt*Nx

        Σ = γ/s
        Σm = s/γ

        Mat = dddv = ∂dyn∂vel(body, dt) + Nx'*Σ*Nv

        Γ = γ#diagm(γ)
        S = s#diagm(s)


        vec = dynamics(body, mechanism) + Nx'*(Σ*φ - γ - μ/s)

        sol = Mat\vec

        Δv = sol
        Δγ = Σ*(φ - Nv*Δv) - μ/s
        Δs = Σm*(γ - Δγ) - μ/γ


        αs = ones(ns)
        αsmax = 1.
        αγ = ones(nγ)
        αγmax = 1.

        for i=1:ns
            if Δs[i] > 0
                αs[i] = minimum([1.;τ*s[i]/Δs[i]])
            end
            αsmax = minimum(αs)

            if Δγ[i] > 0
                αγ[i] = minimum([1.;τ*γ[i]/Δγ[i]])
            end
            αγmax = minimum(αγ)
        end

        v -= αsmax.*Δv
        s -= αsmax.*Δs
        γ -= αγmax.*Δγ

        body.s1 = v
        s1tos0!(body)

        norm1 = norm(vec)

        E = norm1
        μ = σ*μ

        if E<1e-5
            display(iter)
            return
        end
    end
    display("failed")
    return
end
