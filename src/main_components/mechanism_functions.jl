"""
    getbody!(mechanism, id)

Gets the body with ID `id` from `mechanism` if it exists. If `id = nothing`, the origin will be returned.
"""
@inline getbody(mechanism::Mechanism, id::Integer) = mechanism.bodies[id]
@inline getbody(mechanism::Mechanism, id::Nothing) = mechanism.origin

"""
    getbody!(mechanism, name)

Gets the body with name `name` from `mechanism` if it exists.
"""
function getbody(mechanism::Mechanism, name::String)
    if mechanism.origin.name == name
        return mechanism.origin
    else
        for body in mechanism.bodies
            if body.name == name
                return body
            end
        end
    end
    return
end

"""
    geteqconstraint!(mechanism, id)

Gets the equality constraint with ID `id` from `mechanism` if it exists.
"""
@inline geteqconstraint(mechanism::Mechanism, id::Integer) = mechanism.eqconstraints[id]

"""
    geteqconstraint!(mechanism, name)

Gets the equality constraint with name `name` from `mechanism` if it exists.
"""
function geteqconstraint(mechanism::Mechanism, name::String)
    for eqc in mechanism.eqconstraints
        if eqc.name == name
            return eqc
        end
    end
    return
end

@inline getineqconstraint(mechanism::Mechanism, id::Integer) = mechanism.ineqconstraints[id]
function getineqconstraint(mechanism::Mechanism, name::String)
    for ineqc in mechanism.eqconstraints
        if ineqc.name == name
            return ineqc
        end
    end
    return
end

"""
    getcomponent!(mechanism, id)

Gets the component (body or equality constraint) with ID `id` from `mechanism` if it exists.
"""
function getcomponent(mechanism::Mechanism, id::Integer)
    if haskey(mechanism.bodies, id)
        return getbody(mechanism, id)
    elseif haskey(mechanism.eqconstraints, id)
        return geteqconstraint(mechanism, id)
    elseif haskey(mechanism.ineqconstraints, id)
        return getineqconstraint(mechanism, id)
    else
        return
    end
end
getcomponent(mechanism::Mechanism, id::Nothing) = mechanism.origin

"""
    getcomponent!(mechanism, name)

Gets the component (body or equality constraint) with name `name` from `mechanism` if it exists.
"""
function getcomponent(mechanism::Mechanism, name::String)
    component = getbody(mechanism,name)
    if component === nothing
        component = geteqconstraint(mechanism,name)
    end
    if component === nothing
        component = getineqconstraint(mechanism,name)
    end
    return component
end

@inline function normf(mechanism::Mechanism, body::Body)
    f = dynamics(mechanism, body)
    return dot(f, f)
end

@inline function normf(mechanism::Mechanism, eqc::EqualityConstraint)
    f = g(mechanism, eqc)
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

@inline function normf(mechanism::Mechanism)
    mechanism.normf = 0
    
    foreach(addNormf!, mechanism.bodies, mechanism)
    foreach(addNormf!, mechanism.eqconstraints, mechanism)
    foreach(addNormf!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function normf(mechanism::LinearMechanism)
    mechanism.normf = norm([fdynamics(mechanism); fconstraints(mechanism)])

    return sqrt(mechanism.normf)
end

@inline function meritf(mechanism::Mechanism)
    mechanism.normf = 0

    foreach(addNormf!, mechanism.bodies, mechanism)
    foreach(addNormf!, mechanism.eqconstraints, mechanism)
    foreach(addNormfμ!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function normΔs(mechanism::Mechanism)
    mechanism.normΔs = 0

    foreach(addNormΔs!, mechanism.bodies, mechanism)
    foreach(addNormΔs!, mechanism.eqconstraints, mechanism)
    foreach(addNormΔs!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normΔs)
end

@inline function addNormf!(component::Component, mechanism::Mechanism)
    mechanism.normf += normf(mechanism, component)
    return
end

@inline function addNormfμ!(ineqc::InequalityConstraint, mechanism::Mechanism)
    mechanism.normf += normfμ(mechanism, ineqc)
    return
end

@inline function addNormΔs!(component::Component, mechanism::Mechanism)
    mechanism.normΔs += normΔs(component)
    return
end

@inline function discretizestate!(mechanism::Mechanism)
    foreach(discretizestate!, mechanism.bodies, mechanism.Δt)
    return
end

@inline function fdynamics(mechanism::LinearMechanism)
    A = mechanism.A
    Bu = mechanism.Bu
    Bλ = mechanism.Bλ
    z0 = mechanism.z
    z1 = mechanism.zsol[2]
    λ1 = mechanism.λsol[2]
    u = mechanism.u

    return A*z0 + Bu*u + Bλ*λ1 - z1
end

@inline function fconstraints(mechanism::LinearMechanism)
    G = mechanism.G
    z1 = mechanism.zsol[2]

    return G*z1
end

function disassemble(mechanism::Mechanism{T}) where T
    origin = mechanism.origin
    bodies = mechanism.bodies.values
    eqconstraints = mechanism.eqconstraints.values
    ineqconstraints = mechanism.ineqconstraints.values

    # Flip component ids
    for body in bodies
        body.id *= -1
    end
    for eqc in eqconstraints
        eqc.id *= -1
        if eqc.parentid !== nothing
            eqc.parentid *= -1
        end
        eqc.childids *= -1
    end
    for ineqc in ineqconstraints
        ineqc.id *= -1
        if ineqc.parentid !== nothing
            ineqc.parentid *= -1
        end
        ineqc.childids *= -1
    end

    # Set CURRENTID
    global CURRENTID = -1
    for body in bodies
        if body.id <= CURRENTID
            CURRENTID = body.id-1
        end
    end
    for eqc in eqconstraints
        if eqc.id <= CURRENTID
            CURRENTID = eqc.id-1
        end
    end
    for ineqc in ineqconstraints
        if ineqc.id <= CURRENTID
            CURRENTID = ineqc.id-1
        end
    end

    # Set origin to next id
    # oldoid = origin.id
    origin.id = getGlobalID()
    for eqc in eqconstraints
        if eqc.parentid === nothing
            eqc.parentid = origin.id
        end
    end
    for ineqc in ineqconstraints
        if ineqc.parentid === nothing
            ineqc.parentid = origin.id
        end
    end

    return origin, bodies, eqconstraints, ineqconstraints
end
