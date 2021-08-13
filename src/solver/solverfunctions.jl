@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, body::Body)
    matrix_entry.value = ∂dyn∂vel(mechanism, body)
    vector_entry.value = -dynamics(mechanism, body)
    return
end

# @inline function extendDandΔs!(mechanism::Mechanism, diagonal::DiagonalEntry, body::Body, ineqc::InequalityConstraint)
#     diagonal.D += schurD(ineqc, body, mechanism.Δt) # + SMatrix{6,6,Float64,36}(1e-5*I)
#     diagonal.Δs -= schurf(mechanism, ineqc, body)
#     return
# end

@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, eqc::EqualityConstraint)
    μ = 1e-10 # 0 for no regularization
    matrix_entry.value = I*μ
    vector_entry.value = -g(mechanism, eqc)
    return
end

@inline function setLU!(mechanism::Mechanism, matrix_entry_L::Entry, matrix_entry_U::Entry, bodyid::Integer, eqc::EqualityConstraint)
    matrix_entry_L.value = -∂g∂ʳpos(mechanism, eqc, bodyid)'
    matrix_entry_U.value = ∂g∂ʳvel(mechanism, eqc, bodyid)
    return
end

@inline function setLU!(mechanism::Mechanism, matrix_entry_L::Entry, matrix_entry_U::Entry, eqc::EqualityConstraint, bodyid::Integer)
    matrix_entry_L.value = ∂g∂ʳvel(mechanism, eqc, bodyid)
    matrix_entry_U.value = -∂g∂ʳpos(mechanism, eqc, bodyid)'
    return
end

@inline function setLU!(mechanism::Mechanism, matrix_entry_L::Entry, matrix_entry_U::Entry, eqc::EqualityConstraint, body1id::Integer, body2id::Integer)
    D = -offdiagonal∂damper∂ʳvel(mechanism, eqc, body1id, body2id)
    matrix_entry_L.value = D
    matrix_entry_U.value = D'
    return
end

@inline function setLU!(matrix_entry_L::Entry, matrix_entry_U::Entry)
    matrix_entry_L.value *= 0
    matrix_entry_U.value *= 0
    return
end


function feasibilityStepLength!(mechanism::Mechanism)
    ldu = mechanism.ldu

    τ = 0.995
    mechanism.α = 1.

    for ineqc in mechanism.ineqconstraints
        feasibilityStepLength!(mechanism, ineqc, getineqentry(ldu, ineqc.id), τ)
    end

    return
end

# function feasibilityStepLength!(mechanism, ineqc::InequalityConstraint{T,N}, ineqentry::InequalityEntry, τ) where {T,N}
#     s1 = ineqc.ssol[2]
#     γ1 = ineqc.γsol[2]
#     Δs = ineqentry.Δs
#     Δγ = ineqentry.Δγ

#     for i = 1:N
#         αmax = τ * s1[i] / Δs[i]
#         (αmax > 0) && (αmax < mechanism.α) && (mechanism.α = αmax)
#         αmax = τ * γ1[i] / Δγ[i]
#         (αmax > 0) && (αmax < mechanism.α) && (mechanism.α = αmax)
#     end
#     return
# end


function setentries!(mechanism::Mechanism)
    system = mechanism.system

    for (id, body) in pairs(mechanism.bodies)        
        for childid in children(system, id)
            setLU!(mechanism, getentry(system, id, childid), getentry(system, childid, id), id, geteqconstraint(mechanism, childid))
        end

        # for grandchildid in dampergrandchildren(graph, id)
        #     for parentid in predecessors(graph, grandchildid) # Maybe predecessors works out for loop closure?
        #         setLU!(mechanism, getentry(system, id, grandchildid), getentry(system, grandchildid, id), geteqconstraint(mechanism, parentid), id, grandchildid)
        #     end
        # end 

        setDandΔs!(mechanism, getentry(system, id, id), getentry(system, id), body)
        # for childid in ineqchildren(graph, id)
        #     ineqc = getineqconstraint(mechanism, childid)
        #     calcFrictionForce!(mechanism, ineqc)
        #     extendDandΔs!(mechanism, diagonal, body, ineqc)
        # end
    end

    for eqc in mechanism.eqconstraints
        id = eqc.id
        
        for cyclic_children in system.cycles[id]
            for childid in cyclic_children
                setLU!(getentry(system, id, childid), getentry(system, childid, id))
            end
        end

        for childid in children(system, id)
            setLU!(mechanism, getentry(system, id, childid), getentry(system, childid, id), eqc, childid)
        end

        setDandΔs!(mechanism, getentry(system, id, id), getentry(system, id), eqc)
    end

    return 
end

# function eliminatedsolve!(mechanism::Mechanism, ineqentry::InequalityEntry, diagonal::DiagonalEntry, bodyid::Integer, ineqc::InequalityConstraint)
#     μ = mechanism.μ

#     φ = g(mechanism, ineqc)

#     Nv = ∂g∂ʳvel(mechanism, ineqc, bodyid)

#     γ1 = ineqc.γsol[2]
#     s1 = ineqc.ssol[2]

#     Δv = diagonal.Δs
#     ineqentry.Δγ = -γ1 ./ s1 .* φ + μ ./ s1 - γ1 ./ s1 .* (Nv * Δv)
#     ineqentry.Δs = -s1 + μ ./ γ1 - s1 ./ γ1 .* ineqentry.Δγ

#     return
# end


@inline function updatesolution!(body::Body)
    body.state.vsol[1] = body.state.vsol[2]
    body.state.ωsol[1] = body.state.ωsol[2]
    return
end

@inline function updatesolution!(eqc::EqualityConstraint)
    eqc.λsol[1] = eqc.λsol[2]
    return
end

@inline function updatesolution!(ineqc::InequalityConstraint)
    ineqc.ssol[1] = ineqc.ssol[2]
    ineqc.γsol[1] = ineqc.γsol[2]
    return
end

@inline function normΔs(body::Body)
    d1 = body.state.vsol[2] - body.state.vsol[1]
    d2 = body.state.ωsol[2] - body.state.ωsol[1]
    return dot(d1, d1) + dot(d2, d2)
end

@inline function normΔs(eqc::EqualityConstraint)
    d = eqc.λsol[2] - eqc.λsol[1]
    return dot(d, d)
end

@inline function normΔs(ineqc::InequalityConstraint)
    d1 = ineqc.ssol[2] - ineqc.ssol[1]
    d2 = ineqc.γsol[2] - ineqc.γsol[1]
    return dot(d1, d1) + dot(d2, d2)
end


