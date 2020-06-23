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

@inline function normf(mechanism::Mechanism)
    graph = mechanism.graph
    mechanism.normf = 0
    
    foreachactive(addNormf!, graph, mechanism.bodies, mechanism)
    foreachactive(addNormf!, graph, mechanism.eqconstraints, mechanism)
    foreachactive(addNormf!, graph, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function meritf(mechanism::Mechanism)
    graph = mechanism.graph
    mechanism.normf = 0

    foreachactive(addNormf!, graph, mechanism.bodies, mechanism)
    foreachactive(addNormf!, graph, mechanism.eqconstraints, mechanism)
    foreachactive(addNormfμ!, graph, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function normΔs(mechanism::Mechanism)
    graph = mechanism.graph
    mechanism.normΔs = 0

    foreachactive(addNormΔs!, graph, mechanism.bodies, mechanism)
    foreachactive(addNormΔs!, graph, mechanism.eqconstraints, mechanism)
    foreachactive(addNormΔs!, graph, mechanism.ineqconstraints, mechanism)

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
    activate!(mechanism.graph,id)
    return
end

@inline function deactivate!(mechanism::Mechanism, id::Integer)
    deactivate!(mechanism.graph,id)
    return
end