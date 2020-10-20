function saveToStorage!(mechanism::Mechanism, storage::Storage, i)
    for (id, body) in pairs(mechanism.bodies)
        state = body.state
        storage.x[id][i] = state.xc
        storage.q[id][i] = state.qc
        storage.v[id][i] = state.vc
        storage.ω[id][i] = state.ωc
    end
    return
end

function saveToStorage!(mechanism::LinearMechanism, storage::Storage, i)
    qvm = QuatVecMap()
    for (id, body) in pairs(mechanism.bodies)
        storage.x[id][i] = mechanism.xd[id] + mechanism.z[offsetrange(id,3,12,1)]
        storage.q[id][i] = Rotations.add_error(mechanism.qd[id],RotationError(SA[mechanism.z[offsetrange(id,3,12,3)]...],qvm))
        storage.v[id][i] = mechanism.vd[id] + mechanism.z[offsetrange(id,3,12,2)]
        storage.ω[id][i] = mechanism.ωd[id] + mechanism.z[offsetrange(id,3,12,4)]
    end
    return
end

function initializeSimulation!(mechanism::Mechanism, debug::Bool)
    discretizestate!(mechanism)
    debug && verifyConstraints!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    return
end
function initializeSimulation!(mechanism::LinearMechanism, debug::Bool)
    # debug && verifyConstraints!(mechanism)

    qvm = QuatVecMap()
    for (id, body) in pairs(mechanism.bodies)
        state = body.state
        mechanism.z[offsetrange(id,3,12,1)] = state.xc - mechanism.xd[id]
        mechanism.z[offsetrange(id,3,12,3)] = Rotations.rotation_error(state.qc,mechanism.qd[id],qvm)
        mechanism.z[offsetrange(id,3,12,2)] = state.vc - mechanism.vd[id]
        mechanism.z[offsetrange(id,3,12,4)] = state.ωc - mechanism.ωd[id]
    end
    mechanism.zsol[2] = mechanism.z
    return
end


# with control function 
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, control!::Function;
        ε = 1e-10, newtonIter = 100, lineIter = 10,
        record::Bool = false,debug::Bool = false
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, k)
        applyFτ!(mechanism)
        newton!(mechanism, ε = ε, newtonIter = newtonIter, lineIter = lineIter, warning = debug)
        foreachactive(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

function simulate!(mechanism::LinearMechanism, steps::AbstractUnitRange, storage::Storage, control!::Function;
        ε = 1e-10, newtonIter = 100, lineIter = 10,
        record::Bool = false,debug::Bool = false
    )

    initializeSimulation!(mechanism, debug)

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, k)
        applyFτ!(mechanism)
        newton!(mechanism, ε = ε, newtonIter = newtonIter, lineIter = lineIter, warning = debug)
        mechanism.z = mechanism.zsol[2]
        mechanism.λ = mechanism.λsol[2]
    end
    record ? (return storage) : (return) 
end

# with controller
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, controller::Controller;
        ε = 1e-10, newtonIter = 100, lineIter = 10,
        record::Bool = false,debug::Bool = false
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    control! = controller.control!

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, controller, k)
        applyFτ!(mechanism)
        newton!(mechanism, ε = ε, newtonIter = newtonIter, lineIter = lineIter, warning = debug)
        foreachactive(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

# without control
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage;
        ε = 1e-10, newtonIter = 100, lineIter = 10,
        record::Bool = false,debug::Bool = false
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
   
    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        newton!(mechanism, ε = ε, newtonIter = newtonIter, lineIter = lineIter, warning = debug)
        foreachactive(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

function simulate!(mechanism::LinearMechanism, steps::AbstractUnitRange, storage::Storage;
        ε = 1e-10, newtonIter = 100, lineIter = 10,
        record::Bool = false,debug::Bool = false
    )

    initializeSimulation!(mechanism, debug)

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        newton!(mechanism, ε = ε, newtonIter = newtonIter, lineIter = lineIter, warning = debug)
        mechanism.z = mechanism.zsol[2]
        mechanism.λ = mechanism.λsol[2]
    end
    record ? (return storage) : (return) 
end

function simulate!(mechanism::AbstractMechanism{T}, tend::Real, args...;
        ε = 1e-10, newtonIter = 100, lineIter = 10,
        record::Bool = false,debug::Bool = false
    ) where T

    steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
    record ? storage = Storage{T}(steps,length(mechanism.bodies)) : storage = Storage{T}()        
    storage = simulate!(mechanism, steps, storage, args...; ε = ε, newtonIter = newtonIter, lineIter = lineIter, record = record, debug = debug)
    return storage # can be "nothing"
end

function simulate!(mechanism::AbstractMechanism, storage::Storage{T,N}, args...;
        ε = 1e-10, newtonIter = 100, lineIter = 10,
        record::Bool = true, debug::Bool = false
    ) where {T,N}
    
    steps = Base.OneTo(N)
    storage = simulate!(mechanism, steps, storage, args...; ε = ε, newtonIter = newtonIter, lineIter = lineIter, record = record, debug = debug)
    return storage
end