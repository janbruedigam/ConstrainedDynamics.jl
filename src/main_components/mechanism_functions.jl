@inline getbody(mechanism::Mechanism, id::Integer) = mechanism.bodies[id]
@inline getbody(mechanism::Mechanism, id::Nothing) = mechanism.origin
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
@inline geteqconstraint(mechanism::Mechanism, id::Integer) = mechanism.eqconstraints[id]
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

getcomponent(mechanism::Mechanism, id::Nothing) = mechanism.origin
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

function getshape(shapes::Vector{<:Shape}, id)
    for shape in shapes
        for bodyid in shape.bodyids
            if bodyid == id
                return shape
            end
        end
    end

    return
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
    
    foreachactive(addNormf!, mechanism.bodies, mechanism)
    foreachactive(addNormf!, mechanism.eqconstraints, mechanism)
    foreachactive(addNormf!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function normf(mechanism::LinearMechanism)
    mechanism.normf = norm([fdynamics(mechanism); fconstraints(mechanism)])

    return sqrt(mechanism.normf)
end

@inline function meritf(mechanism::Mechanism)
    mechanism.normf = 0

    foreachactive(addNormf!, mechanism.bodies, mechanism)
    foreachactive(addNormf!, mechanism.eqconstraints, mechanism)
    foreachactive(addNormfμ!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function normΔs(mechanism::Mechanism)
    mechanism.normΔs = 0

    foreachactive(addNormΔs!, mechanism.bodies, mechanism)
    foreachactive(addNormΔs!, mechanism.eqconstraints, mechanism)
    foreachactive(addNormΔs!, mechanism.ineqconstraints, mechanism)

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

@inline function currentasknot!(mechanism::Mechanism)
    foreach(currentasknot!, mechanism.bodies)
    return
end

@inline function activate!(mechanism::Mechanism, id::Integer)
    component = getcomponent(mechanism, id)
    activate!(component)
    activate!(mechanism.graph,id)
    return
end

@inline function deactivate!(mechanism::Mechanism, id::Integer)
    component = getcomponent(mechanism, id)
    deactivate!(component)
    deactivate!(mechanism.graph,id)
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

function disassemble(mechanism::Mechanism{T}; shapes::Vector{<:Shape} = Shape{T}[]) where T
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
    for shape in shapes
        for (i,bodyid) in enumerate(shape.bodyids)
            if bodyid != origin.id 
                shape.bodyids[i] *= -1
            end
        end
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
    oldoid = origin.id
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
    for shape in shapes
        for (i,bodyid) in enumerate(shape.bodyids)
            if bodyid == oldoid 
                shape.bodyids[i] = origin.id
            end
        end
    end

    return origin, bodies, eqconstraints, ineqconstraints, shapes
end

@inline function applyFτ!(mechanism::Mechanism{T}; clear::Bool = true) where T
    for eqc in mechanism.eqconstraints
        parentid = eqc.parentid
        for (i,cid) in enumerate(eqc.childids)
            applyFτ!(eqc.constraints[i],getbody(mechanism,parentid),getbody(mechanism,cid))
            clear && (eqc.constraints[i].Fτ = szeros(T,3))
        end
    end

    return
end
