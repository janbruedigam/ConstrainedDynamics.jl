function saveToStorage!(mechanism::Mechanism, storage::Storage, i)
    No = mechanism.No
    for (ind, body) in enumerate(mechanism.bodies)
        storage.x[ind][i] = body.x[No]
        storage.q[ind][i] = body.q[No]
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
        if norm(g(mechanism, eqc)) > 1e-3
            @info string("Probably disconnected bodies at constraint: ", eqc.id)
        end
    end
end

function initializeSimulation!(mechanism::Mechanism, debug::Bool)
    debug && verifyConstraints!(mechanism)
    foreach(s0tos1!, mechanism.bodies)
    foreach(s0tos1!, mechanism.eqconstraints)
    foreach(s0tos1!, mechanism.ineqconstraints)
    return
end

# with control function 
function simulate!(mechanism::Mechanism{T}, steps::AbstractUnitRange, storage::Storage{T}, control!::Function;record::Bool = false,debug::Bool = false) where T
    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    for i = steps
        control!(mechanism, i)
        newton!(mechanism, warning = debug)
        record && saveToStorage!(mechanism, storage, i)
        foreach(updatePos!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

# function simulate!(mechanism::Mechanism{T}, tend::T, control!::Function;record::Bool = false,debug::Bool = false) where T
#     steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
#     record ? storage = Storage{T}(steps,length(mechanism.bodies)) : storage = Storage{T}()        
#     simulate!(mechanism, steps, storage, control!;record=record,debug=debug)
#     return
# end

# function simulate!(mechanism::Mechanism{T}, storage::Storage{T,N}, control!::Function;record::Bool = true,debug::Bool = false) where {T,N}
#     steps = Base.OneTo(N)
#     simulate!(mechanism, steps, storage, control!;record=record,debug=debug)
#     return
# end

# with controller
function simulate!(mechanism::Mechanism{T}, steps::AbstractUnitRange, storage::Storage{T}, controller::Controller;record::Bool = false,debug::Bool = false) where T
    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    control! = controller.control!

    for i = steps
        control!(mechanism, controller, i)
        newton!(mechanism, warning = debug)
        record && saveToStorage!(mechanism, storage, i)
        foreach(updatePos!, bodies, Δt)
    end
    record ? (return storage) : (return) 
end

# function simulate!(mechanism::Mechanism{T}, tend::T, controller::Controller;record::Bool = false,debug::Bool = false) where T
#     steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
#     record ? storage = Storage{T}(steps,length(mechanism.bodies)) : storage = Storage{T}()        
#     simulate!(mechanism, steps, storage, controller;record=record,debug=debug)
#     return
# end

# function simulate!(mechanism::Mechanism{T}, storage::Storage{T,N}, controller::Controller;record::Bool = true,debug::Bool = false) where {T,N}
#     steps = Base.OneTo(N)
#     simulate!(mechanism, steps, storage, controller;record=record,debug=debug)
#     return
# end

# without control
function simulate!(mechanism::Mechanism{T}, steps::AbstractUnitRange, storage::Storage{T}; record::Bool = false,debug::Bool = false) where T
    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
   
    for i = steps
        newton!(mechanism, warning = debug)
        record && saveToStorage!(mechanism, storage, i)
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