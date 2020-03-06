function setentries!(mechanism::Mechanism)
    graph = mechanism.graph
    ldu = mechanism.ldu

    for (id, body) in pairs(mechanism.bodies)
        for cid in directchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)), id, geteqconstraint(mechanism, cid), mechanism)
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(diagonal, body, mechanism)
        for cid in ineqchildren(graph, id)
            extendDandΔs!(diagonal, body, getineqconstraint(mechanism, cid), mechanism)
        end
    end

    for node in mechanism.eqconstraints
        id = node.id

        for cid in directchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)), node, cid, mechanism)
        end

        for cid in loopchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)))
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(diagonal, node, mechanism)
    end
end

@inline function normf(body::Body{T}, mechanism::Mechanism) where T
    f = dynamics(body, mechanism)
    return dot(f, f)
end

@inline function normf(c::EqualityConstraint, mechanism::Mechanism)
    f = g(c, mechanism)
    return dot(f, f)
end

@inline function normf(ineqc::InequalityConstraint, mechanism::Mechanism)
    f = gs(ineqc, mechanism)
    d = h(ineqc)
    return dot(f, f) + dot(d, d)
end

@inline function normfμ(ineqc::InequalityConstraint, mechanism::Mechanism)
    f = gs(ineqc, mechanism)
    d = hμ(ineqc, mechanism.μ)
    return dot(f, f) + dot(d, d)
end

@inline function GtλTof!(body::Body, eqc::EqualityConstraint, mechanism)
    body.f -= ∂g∂pos(eqc, body.id, mechanism)' * eqc.s1
    return
end

@inline function NtγTof!(body::Body, ineqc::InequalityConstraint, mechanism)
    body.f -= ∂g∂pos(ineqc, body, mechanism)' * ineqc.γ1
    return
end

@inline function normf(mechanism::Mechanism)
    mechanism.normf = 0

    # Allocates otherwise
    for body in mechanism.bodies
        mechanism.normf += normf(body, mechanism)
    end
    foreach(addNormf!, mechanism.eqconstraints, mechanism)
    foreach(addNormf!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function meritf(mechanism::Mechanism)
    mechanism.normf = 0

    # Allocates otherwise
    for body in mechanism.bodies
        mechanism.normf += normf(body, mechanism)
    end
    foreach(addNormf!, mechanism.eqconstraints, mechanism)
    foreach(addNormfμ!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function normΔs(mechanism::Mechanism)
    mechanism.normΔs = 0

    # Allocates otherwise
    mechanism.normΔs += mapreduce(normΔs, +, mechanism.bodies)
    foreach(addNormΔs!, mechanism.eqconstraints, mechanism)
    foreach(addNormΔs!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normΔs)
end

@inline function addNormf!(ineqc::InequalityConstraint, mechanism::Mechanism)
    mechanism.normf += normf(ineqc, mechanism)
    return
end

@inline function addNormfμ!(ineqc::InequalityConstraint, mechanism::Mechanism)
    mechanism.normf += normfμ(ineqc, mechanism)
    return
end

@inline function addNormf!(eqc::EqualityConstraint, mechanism::Mechanism)
    mechanism.normf += normf(eqc, mechanism)
    return
end

@inline function addNormΔs!(component::Component, mechanism::Mechanism)
    mechanism.normΔs += normΔs(component)
    return
end

function feasibilityStepLength!(mechanism::Mechanism)
    ldu = mechanism.ldu

    τ = 0.995
    mechanism.α = 1.

    for ineqc in mechanism.ineqconstraints
        feasibilityStepLength!(ineqc, getineqentry(ldu, ineqc.id), τ, mechanism)
    end

    return
end

function feasibilityStepLength!(ineqc::InequalityConstraint{T,N}, ineqentry::InequalityEntry, τ, mechanism) where {T,N}
    s1 = ineqc.s1
    γ1 = ineqc.γ1
    Δs = ineqentry.Δs
    Δγ = ineqentry.Δγ

    for i = 1:N
        αmax = τ * s1[i] / Δs[i]
        (αmax > 0) && (αmax < mechanism.α) && (mechanism.α = αmax)
        αmax = τ * γ1[i] / Δγ[i]
        (αmax > 0) && (αmax < mechanism.α) && (mechanism.α = αmax)
    end
    return
end


function saveToStorage!(mechanism::Mechanism, t)
    No = mechanism.No
    for (ind, body) in enumerate(mechanism.bodies)
        mechanism.storage.x[ind][t] = body.x[No]
        mechanism.storage.q[ind][t] = body.q[No]
    end
    for (ind, eqc) in enumerate(mechanism.eqconstraints)
        mechanism.storage.eqmultipliers[ind][t] = eqc.s1
    end
end

@inline function updatePos!(body::Body, Δt)
    body.x[1] = body.x[2]
    body.x[2] = getx3(body, Δt)
    body.q[1] = body.q[2]
    body.q[2] = getq3(body, Δt)
    return
end

function verifyConstraints!(mechanism::Mechanism)
    for eqc in mechanism.eqconstraints
        if norm(g(eqc,mechanism)) > 1e-3
            @info string("Probably disconnected bodies at constraint: ", eqc.id)
        end
    end
end


function simulate!(mechanism::Mechanism;save::Bool = false,debug::Bool = false)
    debug && verifyConstraints!(mechanism)
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    Δt = mechanism.Δt
    foreach(s0tos1!, bodies)
    foreach(s0tos1!, eqcs)
    foreach(s0tos1!, ineqcs)

    for i = mechanism.steps
        newton!(mechanism, warning = debug)
        save && saveToStorage!(mechanism, i)
        foreach(updatePos!, bodies, Δt)

        # debug && (i*Δt)%1<Δt*(1.0-.1) && display(i*Δt)
    end
    return
end


# function inputcontrol!(mechanism::Mechanism,controller::Controller)
#     for body in mechanism.bodies
#         inputcontrol!(body,controller)
#     end
# end

function simulate!(mechanism::Mechanism,control!::Function;save::Bool = false,debug::Bool = false)
    debug && verifyConstraints!(mechanism)
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    Δt = mechanism.Δt
    foreach(s0tos1!, bodies)
    foreach(s0tos1!, eqcs)
    foreach(s0tos1!, ineqcs)

    for i = mechanism.steps
        control!(mechanism,(i-1)*Δt)
        newton!(mechanism, warning = debug)
        save && saveToStorage!(mechanism, i)
        foreach(updatePos!, bodies, Δt)
    end
    return
end

function simulate!(mechanism::Mechanism,controller::Controller;save::Bool = false,debug::Bool = false)
    debug && verifyConstraints!(mechanism)
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    Δt = mechanism.Δt
    foreach(s0tos1!, bodies)
    foreach(s0tos1!, eqcs)
    foreach(s0tos1!, ineqcs)

    for i = mechanism.steps
        control!(mechanism,controller)
        newton!(mechanism, warning = debug)
        save && saveToStorage!(mechanism, i)
        foreach(updatePos!, bodies, Δt)
    end
    return
end