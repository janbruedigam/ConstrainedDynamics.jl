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
    for ineqc in mechanism.ineqconstraints
        if ineqc.name == name
            return ineqc
        end
    end
    return
end

@inline getfriction(mechanism::Mechanism, id::Integer) = mechanism.frictions[id]
function getfriction(mechanism::Mechanism, name::String)
    for friction in mechanism.frictions
        if friction.name == name
            return friction
        end
    end
    return
end

"""
    getcomponent(mechanism, id)

Gets the component (body or equality constraint) with ID `id` from `mechanism` if it exists.
"""
function getcomponent(mechanism::Mechanism{T,Nn,Nb,Ne,Nf}, id::Integer) where {T,Nn,Nb,Ne,Nf}
    if id <= Ne
        return geteqconstraint(mechanism, id)
    elseif id <= Ne+Nb
        return getbody(mechanism, id)
    elseif id <= Ne+Nb+Nf
        return getineqconstraint(mechanism, id)
    else
        return getfriction(mechanism, id)
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
    if component === nothing
        component = getfriction(mechanism,name)
    end
    return component
end

@inline function normf(mechanism::Mechanism, body::Component)
    f = g(mechanism, body)
    return dot(f, f)
end

@inline function normf(mechanism::Mechanism, ineqc::InequalityConstraint)
    f1 = complementarity(mechanism, ineqc)
    f2 = gs(mechanism, ineqc)
    return dot(f1, f1) + dot(f2, f2)
end

@inline function normfμ(mechanism::Mechanism, ineqc::InequalityConstraint)
    f1 = complementarityμ(mechanism, ineqc)
    f2 = gs(mechanism, ineqc)
    return dot(f1, f1) + dot(f2, f2)
end

@inline function normf(mechanism::Mechanism)
    mechanism.normf = 0
    
    foreach(addNormf!, mechanism.bodies, mechanism)
    foreach(addNormf!, mechanism.eqconstraints, mechanism)
    foreach(addNormf!, mechanism.ineqconstraints, mechanism)

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
