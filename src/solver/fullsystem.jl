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

function linearizeMechanism(mechanism::Mechanism{T,N,Nb}, xc, vc, Fk, qc, ωc, τk) where {T,N,Nb}
    foreach(deactivate!,mechanism.eqconstraints)
    A, B = linearizeSystem(mechanism, xc, vc, Fk, qc, ωc, τk)
    foreach(activate!,mechanism.eqconstraints)
    G = linearizeConstraints(mechanism, xc, vc, Fk, qc, ωc, τk)

    return A, B, G
end

function linearizeSystem(mechanism::Mechanism{T,N,Nb}, xc, vc, Fk, qc, ωc, τk) where {T,N,Nb}
    Δt = mechanism.Δt

    statesold = [State{T}() for i=1:Nb]
    fsold = [(@SVector zeros(T,6)) for i=1:Nb]

    # store old state and set new initial state
    for (i,body) in enumerate(mechanism.bodies)
        col3 = (i-1)*3+1:i*3
        col4 = (i-1)*4+1:i*4
        stateold, fold = settempvars!(body, xc[col3], vc[col3], Fk[col3], qc[col4], ωc[col3], τk[col3], zeros(T,6))
        statesold[i] = stateold
        fsold[i] = fold
    end

    # calculate next state
    discretizestate!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    newton!(mechanism) 

    # get state linearization 
    A = zeros(T,12*Nb,12*Nb)
    B = zeros(T,12*Nb,6*Nb)

    for (i,body) in enumerate(mechanism.bodies)
        col6 = (i-1)*6+1:i*6
        col12 = (i-1)*12+1:i*12
        Ai, Bi = ∂zp1∂z(mechanism, body, Δt)
        A[col12,col12] = Ai
        B[col12,col6] = Bi
    end 

    # restore old state
    for (i,body) in enumerate(mechanism.bodies)
        resettempvars!(body, statesold[i], fsold[i])
    end

    return A, B
end

function linearizeConstraints(mechanism::Mechanism{T,N,Nb}, xc, vc, Fk, qc, ωc, τk) where {T,N,Nb}
    Δt = mechanism.Δt

    statesold = [State{T}() for i=1:Nb]
    fsold = [(@SVector zeros(T,6)) for i=1:Nb]

    # store old state and set new initial state
    for (i,body) in enumerate(mechanism.bodies)
        col3 = (i-1)*3+1:i*3
        col4 = (i-1)*4+1:i*4
        stateold, fold = settempvars!(body, xc[col3], vc[col3], Fk[col3], qc[col4], ωc[col3], τk[col3], zeros(T,6))
        statesold[i] = stateold
        fsold[i] = fold
    end

    # set current state to knot points
    currentasknot!(mechanism)

    # get constraint linearization
    # is ∂g∂vela = 0 correct?
    ne = 0
    for eqc in mechanism.eqconstraints
        ne += length(eqc)
    end
    G = zeros(T,ne,12*Nb)

    ind1 = 1
    ind2 = 0
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        ind2 += length(eqc)

        pid = eqc.parentid
        if pid !== nothing
            colpa = (pid-1)*6+1:pid*6
            colpb = pid*6+1:(pid+1)*6
            G[ind1:ind2,colpa] = ∂g∂posa(mechanism, eqc, pid)
            # G[ind1:ind2,colpb] = ∂g∂vela(mechanism, eqc, pid)
        end
        for cid in eqc.childids
            colca = (cid-1)*12+1:(cid-1)*12+6
            colcb = (cid-1)*12+7:(cid-1)*12+12

            G[ind1:ind2,colca] = ∂g∂posa(mechanism, eqc, cid)
            # G[ind1:ind2,colcb] = ∂g∂vela(mechanism, eqc, cid)
        end

        ind1 = ind2 + 1
    end

    # restore old state
    for (i,body) in enumerate(mechanism.bodies)
        resettempvars!(body, statesold[i], fsold[i])
    end

    return G
end