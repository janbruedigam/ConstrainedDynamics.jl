function saveToStorage!(mechanism::Mechanism, storage::Storage, i)
    Δt = mechanism.Δt
    for (ind, body) in enumerate(mechanism.bodies)
        state = body.state
        storage.x[ind][i] = state.xc
        storage.q[ind][i] = state.qc
        storage.v[ind][i] = state.vc
        storage.ω[ind][i] = state.ωc
    end
    return
end

function verifyConstraints!(mechanism::Mechanism)
    for eqc in mechanism.eqconstraints
        if norm(g(mechanism, eqc)) > 1e-3
            @info string("Bad constraint satisfaction at constraint: ", eqc.id, ", |g| = ", norm(g(mechanism, eqc)))
        end
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
function simulate!(mechanism::Mechanism{T}, steps::AbstractUnitRange, storage::Storage{T}, control!::Function;record::Bool = false,debug::Bool = false) where T
    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    graph = mechanism.graph
    bodies = mechanism.bodies

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, k)
        newton!(mechanism, warning = debug)
        foreachactive(updatestate!, graph, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

# with controller
function simulate!(mechanism::Mechanism{T}, steps::AbstractUnitRange, storage::Storage{T}, controller::Controller;record::Bool = false,debug::Bool = false) where T
    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    graph = mechanism.graph
    bodies = mechanism.bodies

    control! = controller.control!

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, controller, k)
        newton!(mechanism, warning = debug)
        foreachactive(updatestate!, graph, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

# without control
function simulate!(mechanism::Mechanism{T}, steps::AbstractUnitRange, storage::Storage{T}; record::Bool = false,debug::Bool = false) where T
    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    graph = mechanism.graph
    bodies = mechanism.bodies
   
    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        newton!(mechanism, warning = debug)
        foreachactive(updatestate!, graph, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

function simulate!(mechanism::Mechanism{T}, tend::T, args...; record::Bool = false,debug::Bool = false) where T
    steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
    record ? storage = Storage{T}(steps,length(mechanism.bodies)) : storage = Storage{T}()        
    storage = simulate!(mechanism, steps, storage, args...;record=record,debug=debug)
    return storage # can be "nothing"
end

function simulate!(mechanism::Mechanism{T}, storage::Storage{T,N}, args...; record::Bool = true,debug::Bool = false) where {T,N}
    steps = Base.OneTo(N)
    storage = simulate!(mechanism, steps, storage, args...;record=record,debug=debug)
    return storage
end