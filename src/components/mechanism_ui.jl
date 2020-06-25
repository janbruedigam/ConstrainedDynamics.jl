function verifyConstraints!(mechanism::Mechanism)
    for eqc in mechanism.eqconstraints
        if norm(g(mechanism, eqc)) > 1e-3
            @info string("Bad constraint satisfaction at constraint: ", eqc.id, ", |g| = ", norm(g(mechanism, eqc)))
        end
    end
    return
end

function activate!(mechanism::Mechanism)
    graph = mechanism.graph

    for (id,body) in pairs(mechanism.bodies)
        activate!(body)
        activate!(mechanism.graph,id)
    end
    for (id,eqc) in pairs(mechanism.eqconstraints)
        activate!(eqc)
        activate!(mechanism.graph,id)
    end
    for (id,ineqc) in pairs(mechanism.ineqconstraints)
        activate!(ineqc)
        activate!(mechanism.graph,id)
    end

    return
end

function deactivate!(mechanism::Mechanism)
    graph = mechanism.graph

    for (id,body) in pairs(mechanism.bodies)
        deactivate!(body)
        deactivate!(mechanism.graph,id)
    end
    for (id,eqc) in pairs(mechanism.eqconstraints)
        deactivate!(eqc)
        deactivate!(mechanism.graph,id)
    end
    for (id,ineqc) in pairs(mechanism.ineqconstraints)
        deactivate!(ineqc)
        deactivate!(mechanism.graph,id)
    end
    
    return
end

function nameiddict(mechanism::Mechanism)
    dict = Dict{String,Int64}()
    for (id,body) in pairs(mechanism.bodies)
        if body.name != ""
            dict[body.name] = id
        end
    end
    for (id,eqc) in pairs(mechanism.eqconstraints)
        if eqc.name != ""
            dict[eqc.name] = id
        end
    end
    for (id,ineqc) in pairs(mechanism.ineqconstraints)
        if ineqc.name != ""
            dict[ineqc.name] = id
        end
    end

    return dict
end

function minimalCoordinates(mechanism::Mechanism)
    keys = mechanism.eqconstraints.keys
    values = Vector{SVector}()

    for eqc in mechanism.eqconstraints
        push!(values, minimalCoordinates(mechanism, eqc))
    end

    return UnitDict(keys, values)

end

function setPosition!(mechanism::Mechanism, dict::UnitDict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        setPosition!(mechanism, eqc, dict[id])
    end

    return
end

function setVelocity!(mechanism::Mechanism, dict::UnitDict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        setVelocity!(mechanism, eqc, dict[id])
    end

    return
end

function setForce!(mechanism::Mechanism, dict::UnitDict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        setVelocity!(mechanism, eqc, dict[id])
    end

    return
end
