# function newton_ip!(mechanism::Mechanism,body::Body)
#     nv = 6
#     ns = 1
#     nλ = 0
#     nγ = 1
#
#     v = body.s0
#     body.ga0 = rand()
#     body.ga1 = body.ga0
#     body.sl0 = rand()
#     body.sl1 = body.sl0
#
#     mechanism.μ = 1.
#     τ = 0.995
#     σ = 0.2
#
#     dt = mechanism.dt
#
#     for iter=1:100
#         γ = body.ga1
#         s = body.sl1
#         μ = mechanism.μ
#         φ = body.x[2][3]+dt*body.s1[3]
#
#         Nx = [0;0;1.;0;0;0]'
#         Nv = dt*Nx
#
#         Σ = γ/s
#         Σm = s/γ
#
#         Mat = dddv = ∂dyn∂vel(body, dt)
#
#         vec = dynamics(body, mechanism)
#
#
#         sol = Mat\vec
#
#         setentries!(mechanism)
#         factor!(mechanism.graph,mechanism.ldu)
#         solve!(mechanism.graph,mechanism.ldu) # x̂1 for each body and constraint
#
#         computeΔ!(body,getentry(mechanism.ldu,body.id),mechanism)
#
#         # Δv = getentry(mechanism.ldu,body.id).ŝ
#         # Δγ = Σ*(φ - Nv*Δv) - μ/s
#         # Δs = Σm*(γ - Δγ) - μ/γ
#
#         Δv = getentry(mechanism.ldu,body.id).ŝ
#         Δγ = getentry(mechanism.ldu,body.id).ga
#         Δs = getentry(mechanism.ldu,body.id).sl
#
#         # getentry(mechanism.ldu,body.id).ŝ = Δv
#         # getentry(mechanism.ldu,body.id).sl = Δs
#         # getentry(mechanism.ldu,body.id).ga = Δγ
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
#         mechanism.αsmax = αsmax
#         mechanism.αγmax = αγmax
#
#         update!(body,mechanism.ldu,mechanism.αsmax,mechanism.αγmax)
#
#         s1tos0!(body)
#
#         norm1 = norm(vec)
#
#         E = norm1
#         mechanism.μ = σ*mechanism.μ
#
#         if E<1e-5
#             display(iter)
#             return
#         end
#     end
#     display("failed")
#     return
# end


function newton_ip!(mechanism::Mechanism{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, lineIter=10, warning::Bool=false) where {T,Nl}
    # n = 1
    bodies = mechanism.bodies
    constraints = mechanism.equalityconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    dt = mechanism.dt

    for body in mechanism.bodies
        body.sl1 = rand()
        body.sl0 = rand()
        body.ga1 = rand()
        body.ga0 = rand()
    end
    mechanism.μ = 1.
    σ = 0.2

    normf0 = normf(mechanism)
    for n=Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph,ldu)
        solve!(graph,ldu) # x̂1 for each body and constraint

        if mechanism.inequalityconstraints != nothing
            for ineqs in mechanism.inequalityconstraints
                body = getbody(mechanism,ineqs.pid)
                computeΔ!(body,getentry(ldu,body.id),mechanism)
            end
        end
        computeα!(mechanism)

        foreach(update!,bodies,ldu,mechanism.αsmax,mechanism.αγmax)
        foreach(update!,constraints,ldu,mechanism.αsmax,mechanism.αγmax)

        normf1 = normf(mechanism)
        # normf1>normf0 && lineSearch!(mechanism,normf0;iter=lineIter, warning=warning)

        # normΔs not changed yet !!!
        if normΔs(mechanism) < ε && normf1 < ε
            foreach(s1tos0!,bodies)
            foreach(s1tos0!,constraints)
            # display(n)
            return
        else
            foreach(s1tos0!,bodies)
            foreach(s1tos0!,constraints)
            mechanism.μ = σ*mechanism.μ
            normf0=normf1
        end
    end

    if warning
        display(string("WARNING:  newton! did not converge. n = ",newtonIter,", tol = ",normf0,"."))
    end

    return
end
