function formAMatrix(mechanism::Mechanism{T}) where T
    bodies = mechanism.bodies
    eqconstraints = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu

    n = 0
    for body in bodies
        n += length(body)
    end
    for eqc in eqconstraints
        n += length(eqc)
    end

    A = zeros(T,n,n)

    rangeDict = Dict{Int64,UnitRange}()
    n1 = 1
    n2 = 0

    for id in graph.dfslist
        component = getcomponent(mechanism, id)
        n2 += length(component)
        rangeDict[id] = n1:n2


        diagonal = getentry(ldu,id)
        A[n1:n2,n1:n2] = diagonal.D

        for cid in successors(graph, id)
            offdiagonal = getentry(ldu, (id, cid))
            nc1 = first(rangeDict[cid])
            nc2 = last(rangeDict[cid])

            A[n1:n2,nc1:nc2] = offdiagonal.L
            A[nc1:nc2,n1:n2] = offdiagonal.U
        end

        n1 = n2+1
    end    

    return A
end

function formbVector(mechanism::Mechanism{T}) where T
    bodies = mechanism.bodies
    eqconstraints = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu

    n = 0
    for body in bodies
        n += length(body)
    end
    for eqc in eqconstraints
        n += length(eqc)
    end

    b = zeros(T,n)

    n1 = 1
    n2 = 0

    for id in graph.dfslist
        component = getcomponent(mechanism, id)
        n2 += length(component)

        if typeof(component) <:Body
            b[n1:n2] = dynamics(mechanism, component)
        else
            b[n1:n2] = g(mechanism, component)
        end
        

        n1 = n2+1
    end    

    return b
end

function formΔsVector(mechanism::Mechanism{T}) where T
    bodies = mechanism.bodies
    eqconstraints = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu

    n = 0
    for body in bodies
        n += length(body)
    end
    for eqc in eqconstraints
        n += length(eqc)
    end

    Δs = zeros(T,n)

    n1 = 1
    n2 = 0

    for id in graph.dfslist
        component = getcomponent(mechanism, id)
        n2 += length(component)

        diagonal = getentry(ldu,id)
        Δs[n1:n2] = diagonal.Δs

        n1 = n2+1
    end    

    return Δs
end
