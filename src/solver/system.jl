function densesystem(mechanism::Mechanism{T,N,Nb}) where {T,N,Nb}
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu

    n = 6 * Nb
    for eqc in eqcs
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


function ∂g∂ʳextension(mechanism::Mechanism{T,N,Nb}) where {T,N,Nb}
    Δt = mechanism.Δt
    eqcs = mechanism.eqconstraints

    nc = 0
    for eqc in eqcs
        nc += length(eqc)
    end

    Gl = zeros(T,nc,Nb*13)
    Gr = zeros(T,nc,Nb*13)

    oneindc = 0
    for eqc in eqcs
        ind1 = 1
        ind2 = 0

        parentid = eqc.parentid
        if parentid !== nothing
            pbody = getbody(mechanism,parentid)
            pstate = pbody.state
            
            for (i,childid) in enumerate(eqc.childids)
                cbody = getbody(mechanism,childid)
                cstate = cbody.state

                ind2 += getN(eqc.constraints[i])
                range = oneindc+ind1:oneindc+ind2
                pcol3a = offsetrange(parentid,3,13,1)
                pcol3b = offsetrange(parentid,3,13,2)
                pcol3c = offsetrange(parentid,3,13,3)
                pcol3c = first(pcol3c):last(pcol3c)+1
                pcol3d = offsetrange(parentid,3,13,4)
                pcol3d = first(pcol3d)+1:last(pcol3d)+1

                ccol3a = offsetrange(childid,3,13,1)
                ccol3b = offsetrange(childid,3,13,2)
                ccol3c = offsetrange(childid,3,13,3)
                ccol3c = first(ccol3c):last(ccol3c)+1
                ccol3d = offsetrange(childid,3,13,4)
                ccol3d = first(ccol3d)+1:last(ccol3d)+1

                pXl, pQl = ∂g∂posa(eqc.constraints[i], posargsnext(pstate, Δt)..., posargsnext(cstate, Δt)...) # x3
                pGr = ∂g∂ʳposa(eqc.constraints[i], posargssol(pstate)..., posargssol(cstate)...) # x2
                pXr, pQr = pGr[:,1:3], pGr[:,4:6]
                cXl, cQl =  ∂g∂posb(eqc.constraints[i], posargsnext(pstate, Δt)..., posargsnext(cstate, Δt)...) # x3
                cGr =  ∂g∂ʳposb(eqc.constraints[i], posargssol(pstate)..., posargssol(cstate)...) # x2
                cXr, cQr = cGr[:,1:3], cGr[:,4:6]

                if typeof(eqc.constraints[i])<:Joint{T,2}
                    pXl = eqc.constraints[i].V12 * pXl
                    pQl = eqc.constraints[i].V12 * pQl
                    pXr = eqc.constraints[i].V12 * pXr
                    pQr = eqc.constraints[i].V12 * pQr
                    cXl = eqc.constraints[i].V12 * cXl
                    cQl = eqc.constraints[i].V12 * cQl
                    cXr = eqc.constraints[i].V12 * cXr
                    cQr = eqc.constraints[i].V12 * cQr
                elseif typeof(eqc.constraints[i])<:Joint{T,3}
                    pXl = convert(SMatrix{3,3,T,9}, pXl)
                    pQl = convert(SMatrix{3,4,T,12}, pQl)
                    cXl = convert(SMatrix{3,3,T,9}, cXl)
                    cQl = convert(SMatrix{3,4,T,12}, cQl)
                else
                    @error "not supported"
                end

                pGlx = pXl
                pGlq = pQl
                pGrx = pXr
                pGrq = pQr

                Gl[range,pcol3a] = pGlx
                Gl[range,pcol3b] = pGlx*Δt
                Gl[range,pcol3c] = pGlq*Rmat(ωbar(pstate.ωsol[2],Δt))
                Gl[range,pcol3d] = pGlq*Lmat(pstate.qsol[2])*derivωbar(pstate.ωsol[2],Δt)
                Gr[range,pcol3b] = -pGrx
                Gr[range,pcol3d] = -pGrq

                cGlx = cXl
                cGlq = cQl
                cGrx = cXr
                cGrq = cQr

                Gl[range,ccol3a] = cGlx
                Gl[range,ccol3b] = cGlx*Δt
                Gl[range,ccol3c] = cGlq*Rmat(ωbar(cstate.ωsol[2],Δt))
                Gl[range,ccol3d] = cGlq*Lmat(cstate.qsol[2])*derivωbar(cstate.ωsol[2],Δt)
                Gr[range,ccol3b] = cGrx
                Gr[range,ccol3d] = cGrq

                ind1 = ind2+1
            end
        else
            for (i,childid) in enumerate(eqc.childids)
                cbody = getbody(mechanism,childid)
                cstate = cbody.state

                ind2 += getN(eqc.constraints[i])
                range = oneindc+ind1:oneindc+ind2

                ccol3a = offsetrange(childid,3,13,1)
                ccol3b = offsetrange(childid,3,13,2)
                ccol3c = offsetrange(childid,3,13,3)
                ccol3c = first(ccol3c):last(ccol3c)+1
                ccol3d = offsetrange(childid,3,13,4)
                ccol3d = first(ccol3d)+1:last(ccol3d)+1


                cXl, cQl =  ∂g∂posb(eqc.constraints[i], posargsnext(cstate, Δt)...) # x3
                cGr =  ∂g∂ʳposb(eqc.constraints[i], posargssol(cstate)...) # x2
                cXr, cQr = cGr[:,1:3], cGr[:,4:6]

                if typeof(eqc.constraints[i])<:Joint{T,2}
                    cXl = eqc.constraints[i].V12 * cXl
                    cQl = eqc.constraints[i].V12 * cQl
                    cXr = eqc.constraints[i].V12 * cXr
                    cQr = eqc.constraints[i].V12 * cQr
                elseif typeof(eqc.constraints[i])<:Joint{T,3}
                    cXl = convert(SMatrix{3,3,T,9}, cXl)
                    cQl = convert(SMatrix{3,4,T,12}, cQl)
                else
                    @error "not supported"
                end

                cGlx = cXl
                cGlq = cQl
                cGrx = cXr
                cGrq = cQr

                Gl[range,ccol3a] = cGlx
                Gl[range,ccol3b] = cGlx*Δt
                Gl[range,ccol3c] = cGlq*Rmat(ωbar(cstate.ωsol[2],Δt))
                Gl[range,ccol3d] = cGlq*Lmat(cstate.qsol[2])*derivωbar(cstate.ωsol[2],Δt)
                Gr[range,ccol3b] = cGrx
                Gr[range,ccol3d] = cGrq

                ind1 = ind2+1
            end
        end

        oneindc += length(eqc)
    end

    return Gl, -Gr'
end


# TODO only works for 1DOF active constaints (eqcids)
function linearsystem(mechanism::Mechanism{T,N,Nb}, xd, vd, qd, ωd, Fτd, bodyids, eqcids) where {T,N,Nb}
    statesold = [State{T}() for i=1:Nb]

    # store old state and set new initial state
    for (i,id) in enumerate(bodyids)
        stateold = settempvars!(getbody(mechanism, id), xd[i], vd[i], zeros(T,3), qd[i], ωd[i], zeros(T,3), zeros(T,6))
        statesold[i] = stateold
    end
    for (i,id) in enumerate(eqcids)
        setForce!(mechanism, geteqconstraint(mechanism, id), Fτd[i])
    end

    A, B = lineardynamics(mechanism, eqcids) # TODO check again for Fk , τk
    λsol2 = [eqc.λsol[2] for eqc in mechanism.eqconstraints]

    # restore old state
    for (i,id) in enumerate(bodyids)
        body = getbody(mechanism, id)
        body.state = statesold[i]
    end

    return A, B, λsol2
end


function lineardynamics(mechanism::Mechanism{T,N,Nb}, eqcids) where {T,N,Nb}
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    nc = 0
    for eqc in eqcs
        nc += length(eqc)
    end

    # calculate next state
    discretizestate!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    newton!(mechanism) 

    # get state linearization 
    Fz = zeros(T,Nb*13+nc,Nb*12+nc)
    Ffz = zeros(T,Nb*13+nc,Nb*13+nc)
    Dinv = zeros(T,Nb*12+nc,Nb*13+nc)
    Bbody = zeros(T,Nb*13+nc,Nb*6)
    Bcontrol = zeros(T,Nb*6,length(eqcids))

    for (id,body) in pairs(bodies)
        col6 = offsetrange(id,6)
        col12 = offsetrange(id,12)
        col13 = offsetrange(id,13)
        col18 = offsetrange(id,18)

        Ai, Bi = ∂F∂z(mechanism, body, Δt)
        Fz[col13,col12] = Ai
        Bbody[col13,col6] = Bi

        Ai, Dinvi = ∂F∂fz(mechanism, body, Δt)
        Ffz[col13,col13] = Ai
        Dinv[col12,col13] = Dinvi
    end
    for (i,id) in enumerate(eqcids)
        eqc = geteqconstraint(mechanism, id)
        parentid = eqc.parentid
        if parentid !== nothing
            col6 = offsetrange(parentid,6)
            Bcontrol[col6,i] = ∂Fτ∂ua(mechanism, eqc, parentid)
        end
        for childid in eqc.childids
            col6 = offsetrange(childid,6)
            Bcontrol[col6,i] = ∂Fτ∂ub(mechanism, eqc, childid)
        end
    end

    # Fz = Fz # no addition necessary

    Gl, mGrt = ∂g∂ʳextension(mechanism)
    Ffz[Nb*13+1:Nb*13+nc,1:Nb*13] = Gl
    Ffz[1:Nb*13,Nb*13+1:Nb*13+nc] = mGrt
    Dinv[Nb*12+1:Nb*12+nc,Nb*13+1:Nb*13+nc] = diagm(ones(T,nc))

    A = -Dinv*inv(Ffz)*Fz
    B = -Dinv*inv(Ffz)*Bbody*Bcontrol

    A = A[1:Nb*12,1:Nb*12]
    B = B[1:Nb*12,:]


    return A, B
end
