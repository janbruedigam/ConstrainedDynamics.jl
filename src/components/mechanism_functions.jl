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
    return nothing
end
@inline geteqconstraint(mechanism::Mechanism, id::Integer) = mechanism.eqconstraints[id]
function geteqconstraint(mechanism::Mechanism, name::String)
    for eqc in mechanism.eqconstraints
        if eqc.name == name
            return eqc
        end
    end
    return nothing
end
@inline getineqconstraint(mechanism::Mechanism, id::Integer) = mechanism.ineqconstraints[id]
function getineqconstraint(mechanism::Mechanism, name::String)
    for ineqc in mechanism.eqconstraints
        if ineqc.name == name
            return ineqc
        end
    end
    return nothing
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
        return nothing
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

function getshape(shapes::Vector{Shape}, id)
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


# function plotθ(mechanism::Mechanism{T}, id) where T
#     n = length(mechanism.bodies)
#     θ = zeros(T, n, length(mechanism.steps))
#     for i = 1:n
#         qs = mechanism.storage.q[i]
#         for (t, q) in enumerate(qs)
#             θ[i,t] = angleaxis(q)[1] * sign(angleaxis(q)[2][1])
#         end
#     end

#     p = plot(collect(0:mechanism.Δt:mechanism.tend - mechanism.Δt), θ[id[1],:])
#     for ind in Iterators.rest(id, 2)
#         plot!(collect(0:mechanism.Δt:mechanism.tend - mechanism.Δt), θ[ind,:])
#     end
#     return p
# end

# function plotλ(mechanism::Mechanism{T}, id) where T
#     n = sum(length.(mechanism.eqconstraints))
#     λ = zeros(T, n, length(mechanism.steps))
#     startpos = 1
#     endpos = 0
#     for i = 1:length(mechanism.eqconstraints)
#         endpos = startpos + length(mechanism.eqconstraints[i]) - 1

#         λs = mechanism.storage.λ[i]
#         for (t, val) in enumerate(λs)
#             λ[startpos:endpos,t] = val
#         end

#         startpos = endpos + 1
#     end

#     p = plot(collect(0:mechanism.Δt:mechanism.tend - mechanism.Δt), λ[id[1],:])
#     for ind in Iterators.rest(id, 2)
#         plot!(collect(0:mechanism.Δt:mechanism.tend - mechanism.Δt), λ[ind,:])
#     end
#     return p
# end