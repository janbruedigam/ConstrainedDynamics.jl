@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, component::Component)
    matrix_entry.value = ∂g∂ʳself(mechanism, component)
    vector_entry.value = -g(mechanism, component)
    return
end
@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, ineqc::InequalityConstraint)
    matrix_entry.value = [[diagm(ineqc.γsol[2]);-I] [diagm(ineqc.ssol[2]);I*0]]
    vector_entry.value = [-complementarityμ(mechanism, ineqc);-gs(mechanism, ineqc)]
    return
end

@inline function setLU!(mechanism::Mechanism, matrix_entry_L::Entry, matrix_entry_U::Entry, componenta::Component, componentb::Component)
    L, U = ∂gab∂ʳba(mechanism, componenta, componentb)
    matrix_entry_L.value = L
    matrix_entry_U.value = U
    return
end


# @inline function setLU!(mechanism::Mechanism, matrix_entry_L::Entry, matrix_entry_U::Entry, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}, body::Body) where {T,N,Nc,Cs,N½}
#     # A = ∂g∂ʳvel(mechanism, ineqc, id)
#     # Z = A*0
#     matrix_entry_L.value = SA[0.0 0.0;0.0 ineqc.constraints[1].cf]#[Z;A]
#     matrix_entry_U.value = SA[0.0 0.0;0.0 0.0]#[Z;-∂g∂ʳpos(mechanism, ineqc, id)]'
#     return
# end




@inline function zeroLU!(matrix_entry_L::Entry{ET1}, matrix_entry_U::Entry{ET2}) where {ET1,ET2}
    matrix_entry_L.value = zero(ET1)
    matrix_entry_U.value = zero(ET2)
    return
end


function feasibilityStepLength!(mechanism::Mechanism)
    system = mechanism.system
    mechanism.α = 1.0

    for ineqc in mechanism.ineqconstraints
        feasibilityStepLength!(mechanism, ineqc, getentry(system, ineqc.id))
    end

    return
end

function feasibilityStepLength!(mechanism, ineqc::InequalityConstraint{T,N}, vector_entry::Entry) where {T,N}
    τ = mechanism.τ
    sγ = [ineqc.ssol[2];ineqc.γsol[2]]
    Δsγ = vector_entry.value

    for i = 1:N
        αmax = τ * sγ[i] / Δsγ[i]
        (αmax > 0) && (αmax < mechanism.α) && (mechanism.α = αmax)
    end
    return
end


function setentries!(mechanism::Mechanism)
    system = mechanism.system

    for id in reverse(system.dfs_list)
        for childid in system.cyclic_children[id]
            zeroLU!(getentry(system, id, childid), getentry(system, childid, id))
        end

        component = getcomponent(mechanism, id)
        setDandΔs!(mechanism, getentry(system, id, id), getentry(system, id), component)

        for childid in children(system,id)
            setLU!(mechanism, getentry(system, id, childid), getentry(system, childid, id), component, getcomponent(mechanism, childid))
        end
    end

    return 
end


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

@inline function updatesolution!(fric::Friction)
    fric.βsol[1] = fric.βsol[2]
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


