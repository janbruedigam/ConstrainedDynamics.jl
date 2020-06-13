function densesystem(mechanism::Mechanism{T,N,Nb}) where {T,N,Nb}
    bodies = mechanism.bodies
    eqconstraints = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu

    n = 6 * Nb
    for eqc in eqconstraints
        n += length(eqc)
    end

    A = zeros(T,n,n)
    x = zeros(T,n)
    b = zeros(T,n)
    
    rangeDict = Dict{Int64,UnitRange}()
    ind1 = 1
    ind2 = 0

    for id in graph.dfslist
        component = getcomponent(mechanism, id)
        ind2 += length(component)
        range = ind1:ind2
        rangeDict[id] = range


        # A
        diagonal = getentry(ldu,id)
        A[range,range] = diagonal.D

        for childid in successors(graph, id)
            offdiagonal = getentry(ldu, (id, childid))
            nc1 = first(rangeDict[childid])
            nc2 = last(rangeDict[childid])

            A[range,nc1:nc2] = offdiagonal.L
            A[nc1:nc2,range] = offdiagonal.U
        end

        # x
        x[range] = diagonal.Δs

        # b
        if typeof(component) <:Body
            b[range] = dynamics(mechanism, component)
        else
            b[range] = g(mechanism, component)
        end


        ind1 = ind2+1
    end    
    
    return A, x, b
end

# TODO only works for 1DOF active constaints (eqcids)
function linearsystem(mechanism::Mechanism{T,N,Nb}, xd, vd, qd, ωd, Fτ, bodyids, eqcids) where {T,N,Nb}
    statesold = [State{T}() for i=1:Nb]

    # store old state and set new initial state
    for (i,id) in enumerate(bodyids)
        stateold = settempvars!(getbody(mechanism, id), xd[i], vd[i], zeros(T,3), qd[i], ωd[i], zeros(T,3), zeros(T,6))
        statesold[i] = stateold
    end
    for (i,id) in enumerate(eqcids)
        setForce!(mechanism, geteqconstraint(mechanism, id), [Fτ[i]])
    end

    A, B = lineardynamics(mechanism, eqcids) # TODO check again for Fk , τk

    # reset new initial state
    for (i,id) in enumerate(bodyids)
        settempvars!(getbody(mechanism, id), xd[i], vd[i], zeros(T,3), qd[i], ωd[i], zeros(T,3), zeros(T,6))
    end
    for (i,id) in enumerate(eqcids)
        setForce!(mechanism, geteqconstraint(mechanism, id), [Fτ[i]])
    end

    G = linearconstraints(mechanism) # TODO check again for Fk , τk

    # restore old state
    for (i,id) in enumerate(bodyids)
        body = getbody(mechanism, id)
        body.state = statesold[i]
    end

    return A, B, G
end

function lineardynamics(mechanism::Mechanism{T,N,Nb}, eqcids) where {T,N,Nb}
    Δt = mechanism.Δt

    # calculate next state
    discretizestate!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    newton!(mechanism) 

    # get state linearization 
    A = zeros(T,12*Nb,12*Nb)
    Bbody = zeros(T,12*Nb,6*Nb)
    Bcontrol = zeros(T,6*Nb,length(eqcids))

    for (id,body) in enumerate(mechanism.bodies)
        col6 = (id-1)*6+1:id*6
        col12 = (id-1)*12+1:id*12
        Ai, Bi = ∂zp1∂z(mechanism, body, Δt)
        A[col12,col12] = Ai
        Bbody[col12,col6] = Bi
    end 
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        parentid = eqc.parentid
        if parentid !== nothing
            col6 = (parentid-1)*6+1:parentid*6
            Bcontrol[col6,i] = ∂Fτ∂ua(mechanism, eqc, parentid)
        end
        for childid in eqc.childids
            col6 = (childid-1)*6+1:childid*6
            Bcontrol[col6,i] = ∂Fτ∂ub(mechanism, eqc, childid)
        end
    end

    B = Bbody * Bcontrol

    return A, B
end

function linearconstraints(mechanism::Mechanism{T,N,Nb}) where {T,N,Nb}
    Δt = mechanism.Δt

    # set current state to knot points
    currentasknot!(mechanism)

    # TODO
    # get constraint linearization
    # is ∂g∂vel = 0 correct?
    neqcs = 0
    for eqc in mechanism.eqconstraints
        neqcs += length(eqc)
    end
    G = zeros(T,neqcs,12*Nb)

    ind1 = 1
    ind2 = 0
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        ind2 += length(eqc)

        parentid = eqc.parentid
        if parentid !== nothing
            colpa = (parentid-1)*6+1:parentid*6
            colpb = parentid*6+1:(parentid+1)*6
            
            G[ind1:ind2,colpa] = ∂g∂posa(mechanism, eqc, parentid)
            # G[ind1:ind2,colpb] = ∂g∂vela(mechanism, eqc, parentid)
        end
        for childid in eqc.childids
            colca = (childid-1)*12+1:(childid-1)*12+6
            colcb = (childid-1)*12+7:(childid-1)*12+12

            G[ind1:ind2,colca] = ∂g∂posb(mechanism, eqc, childid)
            # G[ind1:ind2,colcb] = ∂g∂velb(mechanism, eqc, childid)
        end

        ind1 = ind2 + 1
    end

    return G
end