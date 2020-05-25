function saveToStorage!(mechanism::Mechanism, storage::Storage, i)
    Δt = mechanism.Δt
    for (ind, body) in enumerate(mechanism.bodies)
        storage.x[ind][i] = getx1(body)
        storage.q[ind][i] = getq1(body)
        storage.v[ind][i] = getv1(body)
        storage.ω[ind][i] = getω1(body)
    end
end

@inline function updatePos!(body::Body, Δt)
    body.state.xd[1] = body.state.xd[2]
    body.state.xd[2] = getx2(body, Δt)
    body.state.qd[1] = body.state.qd[2]
    body.state.qd[2] = getq2(body, Δt)
    body.state.vc[1] = body.solv
    body.state.ωc[1] = body.solω
    return
end

function verifyConstraints!(mechanism::Mechanism)
    for eqc in mechanism.eqconstraints
        if norm(g(mechanism, eqc)) > 1e-3
            @info string("Probably disconnected bodies at constraint: ", eqc.id)
        end
    end
end

function initializeSimulation!(mechanism::Mechanism, debug::Bool)
    discretizestate!(mechanism)
    debug && verifyConstraints!(mechanism)
    foreach(setsol!, mechanism.bodies)
    # foreach(s0tos1!, mechanism.eqconstraints)
    # foreach(s0tos1!, mechanism.ineqconstraints)
    return
end

# with control function 
function simulate!(mechanism::Mechanism{T}, steps::AbstractUnitRange, storage::Storage{T}, control!::Function;record::Bool = false,debug::Bool = false) where T
    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, k)
        newton!(mechanism, warning = debug)
        foreach(updatePos!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

# with controller
function simulate!(mechanism::Mechanism{T}, steps::AbstractUnitRange, storage::Storage{T}, controller::Controller;record::Bool = false,debug::Bool = false) where T
    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    control! = controller.control!

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, controller, k)
        newton!(mechanism, warning = debug)
        foreach(updatePos!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

# without control
function simulate!(mechanism::Mechanism{T}, steps::AbstractUnitRange, storage::Storage{T}; record::Bool = false,debug::Bool = false) where T
    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
   
    for k = steps
        newton!(mechanism, warning = debug)
        record && saveToStorage!(mechanism, storage, k)
        foreach(updatePos!, bodies, Δt)
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