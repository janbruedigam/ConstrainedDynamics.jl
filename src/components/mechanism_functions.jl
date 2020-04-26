function setentries!(mechanism::Mechanism)
    graph = mechanism.graph
    ldu = mechanism.ldu

    for (id, body) in pairs(mechanism.bodies)
        for cid in directchildren(graph, id)
            setLU!(mechanism, getentry(ldu, (id, cid)), id, geteqconstraint(mechanism, cid))
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(mechanism, diagonal, body)
        for cid in ineqchildren(graph, id)
            extendDandΔs!(mechanism, diagonal, body, getineqconstraint(mechanism, cid))
        end
    end

    for node in mechanism.eqconstraints
        id = node.id

        for cid in directchildren(graph, id)
            setLU!(mechanism, getentry(ldu, (id, cid)), node, cid)
        end

        for cid in loopchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)))
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(mechanism, diagonal, node)
    end
end

@inline function normf(mechanism::Mechanism, body::Body{T}) where T
    f = dynamics(mechanism, body)
    return dot(f, f)
end

@inline function normf(mechanism::Mechanism, c::EqualityConstraint)
    f = g(mechanism, c)
    return dot(f, f)
end

@inline function normf(mechanism::Mechanism, ineqc::InequalityConstraint)
    f = gs(mechanism, ineqc)
    d = h(ineqc)
    return dot(f, f) + dot(d, d)
end

@inline function normfμ(mechanism::Mechanism, ineqc::InequalityConstraint)
    f = gs(mechanism, ineqc)
    d = hμ(ineqc, mechanism.μ)
    return dot(f, f) + dot(d, d)
end

@inline function GtλTof!(mechanism::Mechanism, body::Body, eqc::EqualityConstraint)
    body.f -= ∂g∂pos(mechanism, eqc, body.id)' * eqc.s1
    return
end

@inline function NtγTof!(mechanism::Mechanism, body::Body, ineqc::InequalityConstraint)
    body.f -= ∂g∂pos(mechanism, ineqc, body)' * ineqc.γ1
    return
end

@inline function normf(mechanism::Mechanism)
    mechanism.normf = 0

    # Allocates otherwise
    for body in mechanism.bodies
        mechanism.normf += normf(mechanism, body)
    end
    foreach(addNormf!, mechanism.eqconstraints, mechanism)
    foreach(addNormf!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function meritf(mechanism::Mechanism)
    mechanism.normf = 0

    # Allocates otherwise
    for body in mechanism.bodies
        mechanism.normf += normf(mechanism, body)
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
    mechanism.normf += normf(mechanism, ineqc)
    return
end

@inline function addNormfμ!(ineqc::InequalityConstraint, mechanism::Mechanism)
    mechanism.normf += normfμ(mechanism, ineqc)
    return
end

@inline function addNormf!(eqc::EqualityConstraint, mechanism::Mechanism)
    mechanism.normf += normf(mechanism, eqc)
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
        feasibilityStepLength!(mechanism, ineqc, getineqentry(ldu, ineqc.id), τ)
    end

    return
end

function feasibilityStepLength!(mechanism, ineqc::InequalityConstraint{T,N}, ineqentry::InequalityEntry, τ) where {T,N}
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

function setFrictionForce!(mechanism)
    foreach(setFrictionForce!, mechanism.ineqconstraints, mechanism)
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
        if norm(g(eqc, mechanism)) > 1e-3
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
    end
    return
end


function simulate!(mechanism::Mechanism, control!::Function;save::Bool = false,debug::Bool = false)
    debug && verifyConstraints!(mechanism)
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    Δt = mechanism.Δt
    foreach(s0tos1!, bodies)
    foreach(s0tos1!, eqcs)
    foreach(s0tos1!, ineqcs)

    for i = mechanism.steps
        control!(mechanism, (i - 1) * Δt)
        newton!(mechanism, warning = debug)
        save && saveToStorage!(mechanism, i)
        foreach(updatePos!, bodies, Δt)
    end
    return
end

function simulate!(mechanism::Mechanism, controller::Controller;save::Bool = false,debug::Bool = false)
    debug && verifyConstraints!(mechanism)
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    Δt = mechanism.Δt
    foreach(s0tos1!, bodies)
    foreach(s0tos1!, eqcs)
    foreach(s0tos1!, ineqcs)

    for i = mechanism.steps
        control!(mechanism, controller)
        newton!(mechanism, warning = debug)
        save && saveToStorage!(mechanism, i)
        foreach(updatePos!, bodies, Δt)
    end
    return
end