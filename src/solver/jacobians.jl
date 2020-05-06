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

function ∂dyn∂pos(mechanism::Mechanism{T}) where T
    bodies = mechanism.bodies
    
    n = 0
    for body in bodies
        n += length(body)
    end

    A = zeros(T,n,n)

    for body in bodies
        id = body.id
        A[(id-1)*6+1:(id-1)*6+6,(id-1)*6+1:(id-1)*6+6] = ∂dyn∂pos(body, mechanism.Δt)
    end

    return A
end

function ∂dyn∂vel(mechanism::Mechanism{T}) where T
    bodies = mechanism.bodies
    
    n = 0
    for body in bodies
        n += length(body)
    end

    A = zeros(T,n,n)

    for body in bodies
        id = body.id
        A[(id-1)*6+1:(id-1)*6+6,(id-1)*6+1:(id-1)*6+6] = ∂dyn∂vel(body, mechanism.Δt)
    end

    return A
end

function ∂dyn∂con(mechanism::Mechanism{T}) where T
    bodies = mechanism.bodies
    
    n = 0
    for body in bodies
        n += length(body)
    end

    A = zeros(T,n,n)

    for body in bodies
        id = body.id
        A[(id-1)*6+1:(id-1)*6+6,(id-1)*6+1:(id-1)*6+6] = ∂dyn∂con(body, mechanism.Δt)
    end

    return A
end

function ∂g∂pos(mechanism::Mechanism{T}) where T
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    
    nb = 0
    nc = 0
    for body in bodies
        nb += length(body)
    end
    for eqc in eqcs
        nc += length(eqc)
    end

    A = zeros(T,nc,nb)

    n1 = 1
    n2 = 0
    for eqc in eqcs
        n2 += length(eqc)

        pid = eqc.pid
        bodyids = unique(eqc.bodyids)

        pid!=nothing && (A[n1:n2,(pid-1)*6+1:(pid-1)*6+6] = ∂g∂pos(mechanism,eqc,pid))
        for id in bodyids
            A[n1:n2,(id-1)*6+1:(id-1)*6+6] = ∂g∂pos(mechanism,eqc,id)
        end

        n1 = n2+1
    end

    return A
end

function ∂g∂vel(mechanism::Mechanism{T}) where T
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    
    nb = 0
    nc = 0
    for body in bodies
        nb += length(body)
    end
    for eqc in eqcs
        nc += length(eqc)
    end

    A = zeros(T,nc,nb)

    n1 = 1
    n2 = 0
    for eqc in eqcs
        n2 += length(eqc)

        pid = eqc.pid
        bodyids = unique(eqc.bodyids)

        pid!=nothing && (A[n1:n2,(pid-1)*6+1:(pid-1)*6+6] = ∂g∂vel(mechanism,eqc,pid))
        for id in bodyids
            A[n1:n2,(id-1)*6+1:(id-1)*6+6] = ∂g∂vel(mechanism,eqc,id)
        end

        n1 = n2+1
    end

    return A
end

function ∂g∂con(mechanism::Mechanism{T}) where T
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    
    nb = 0
    nc = 0
    for body in bodies
        nb += length(body)
    end
    for eqc in eqcs
        nc += length(eqc)
    end

    A = zeros(T,nc,nb)

    return A
end