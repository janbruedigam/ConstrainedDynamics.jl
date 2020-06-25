function saveToStorage!(mechanism::Mechanism, storage::Storage, i)
    Δt = mechanism.Δt
    for (id, body) in pairs(mechanism.bodies)
        state = body.state
        storage.x[id][i] = state.xc
        storage.q[id][i] = state.qc
        storage.v[id][i] = state.vc
        storage.ω[id][i] = state.ωc
    end
    return
end

function initializeSimulation!(mechanism::Mechanism, debug::Bool)
    discretizestate!(mechanism)
    debug && verifyConstraints!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    return
end

# with control function 
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, control!::Function;
        record::Bool = false,debug::Bool = false
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, k)
        newton!(mechanism, warning = debug)
        foreachactive(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

# with controller
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, controller::Controller;
        record::Bool = false,debug::Bool = false
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    control! = controller.control!

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, controller, k)
        newton!(mechanism, warning = debug)
        foreachactive(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

# without control
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage;
        record::Bool = false,debug::Bool = false
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
   
    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        newton!(mechanism, warning = debug)
        foreachactive(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

function simulate!(mechanism::Mechanism{T}, tend::Real, args...;
        record::Bool = false,debug::Bool = false
    ) where T

    steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
    record ? storage = Storage{T}(steps,length(mechanism.bodies)) : storage = Storage{T}()        
    storage = simulate!(mechanism, steps, storage, args...;record=record,debug=debug)
    return storage # can be "nothing"
end

function simulate!(mechanism::Mechanism, storage::Storage{T,N}, args...;
        record::Bool = true,debug::Bool = false
    ) where {T,N}
    
    steps = Base.OneTo(N)
    storage = simulate!(mechanism, steps, storage, args...;record=record,debug=debug)
    return storage
end