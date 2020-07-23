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
        if component isa Body
            b[range] = dynamics(mechanism, component)
        else
            b[range] = g(mechanism, component)
        end


        ind1 = ind2+1
    end    
    
    return A, x, b
end

# TODO Gl and Gr probably wrong!
function ∂g∂ʳextension(mechanism::Mechanism{T,N,Nb}) where {T,N,Nb}
    Δt = mechanism.Δt
    eqcs = mechanism.eqconstraints

    nc = 0
    for eqc in eqcs
        isinactive(eqc) && continue
        nc += length(eqc)
    end

    Gl = zeros(T,nc,Nb*12)
    Gr = zeros(T,nc,Nb*13)

    oneindc = 0
    for eqc in eqcs
        isinactive(eqc) && continue
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

                pcol3a12 = offsetrange(parentid,3,12,1)
                pcol3b12 = offsetrange(parentid,3,12,2)
                pcol3c12 = offsetrange(parentid,3,12,3)
                pcol3d12 = offsetrange(parentid,3,12,4)

                pcol3b = offsetrange(parentid,3,13,2)
                pcol3d = offsetrange(parentid,3,13,4)
                pcol3d = first(pcol3d)+1:last(pcol3d)+1

                ccol3a12 = offsetrange(childid,3,12,1)
                ccol3b12 = offsetrange(childid,3,12,2)
                ccol3c12 = offsetrange(childid,3,12,3)
                ccol3d12 = offsetrange(childid,3,12,4)

                ccol3b = offsetrange(childid,3,13,2)
                ccol3d = offsetrange(childid,3,13,4)
                ccol3d = first(ccol3d)+1:last(ccol3d)+1

                pXl, pQl = ∂g∂posa(eqc.constraints[i], posargsnext(pstate, Δt)..., posargsnext(cstate, Δt)...) # x3
                pGr = ∂g∂ʳposa(eqc.constraints[i], posargssol(pstate)..., posargssol(cstate)...) # x2
                pXr, pQr = pGr[:,1:3], pGr[:,4:6]
                cXl, cQl =  ∂g∂posb(eqc.constraints[i], posargsnext(pstate, Δt)..., posargsnext(cstate, Δt)...) # x3
                cGr =  ∂g∂ʳposb(eqc.constraints[i], posargssol(pstate)..., posargssol(cstate)...) # x2
                cXr, cQr = cGr[:,1:3], cGr[:,4:6]

                if eqc.constraints[i] isa Joint{T,2}
                    pXl = eqc.constraints[i].V12 * pXl
                    pQl = eqc.constraints[i].V12 * pQl
                    pXr = eqc.constraints[i].V12 * pXr
                    pQr = eqc.constraints[i].V12 * pQr
                    cXl = eqc.constraints[i].V12 * cXl
                    cQl = eqc.constraints[i].V12 * cQl
                    cXr = eqc.constraints[i].V12 * cXr
                    cQr = eqc.constraints[i].V12 * cQr
                elseif eqc.constraints[i] isa Joint{T,3}
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

                Gl[range,pcol3a12] = pGlx
                Gl[range,pcol3b12] = pGlx*Δt
                Gl[range,pcol3c12] = pGlq*Rmat(ωbar(pstate.ωsol[2],Δt))*LVᵀmat(pstate.qsol[2])
                Gl[range,pcol3d12] = pGlq*Lmat(pstate.qsol[2])*derivωbar(pstate.ωsol[2],Δt)
                Gr[range,pcol3b] = pGrx
                Gr[range,pcol3d] = pGrq

                cGlx = cXl
                cGlq = cQl
                cGrx = cXr
                cGrq = cQr

                Gl[range,ccol3a12] = cGlx
                Gl[range,ccol3b12] = cGlx*Δt
                Gl[range,ccol3c12] = cGlq*Rmat(ωbar(cstate.ωsol[2],Δt))*LVᵀmat(cstate.qsol[2])
                Gl[range,ccol3d12] = cGlq*Lmat(cstate.qsol[2])*derivωbar(cstate.ωsol[2],Δt)
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

                ccol3a12 = offsetrange(childid,3,12,1)
                ccol3b12 = offsetrange(childid,3,12,2)
                ccol3c12 = offsetrange(childid,3,12,3)
                ccol3d12 = offsetrange(childid,3,12,4)

                ccol3b = offsetrange(childid,3,13,2)
                ccol3d = offsetrange(childid,3,13,4)
                ccol3d = first(ccol3d)+1:last(ccol3d)+1


                cXl, cQl =  ∂g∂posb(eqc.constraints[i], posargsnext(cstate, Δt)...) # x3
                cGr =  ∂g∂ʳposb(eqc.constraints[i], posargssol(cstate)...) # x2
                cXr, cQr = cGr[:,1:3], cGr[:,4:6]

                if eqc.constraints[i] isa Joint{T,2}
                    cXl = eqc.constraints[i].V12 * cXl
                    cQl = eqc.constraints[i].V12 * cQl
                    cXr = eqc.constraints[i].V12 * cXr
                    cQr = eqc.constraints[i].V12 * cQr
                elseif eqc.constraints[i] isa Joint{T,3}
                    cXl = convert(SMatrix{3,3,T,9}, cXl)
                    cQl = convert(SMatrix{3,4,T,12}, cQl)
                else
                    @error "not supported"
                end

                cGlx = cXl
                cGlq = cQl
                cGrx = cXr
                cGrq = cQr

                Gl[range,ccol3a12] = cGlx
                Gl[range,ccol3b12] = cGlx*Δt
                Gl[range,ccol3c12] = cGlq*Rmat(ωbar(cstate.ωsol[2],Δt))*LVᵀmat(cstate.qsol[2])
                Gl[range,ccol3d12] = cGlq*Lmat(cstate.qsol[2])*derivωbar(cstate.ωsol[2],Δt)
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

    A, Bu, Bλ, G = lineardynamics(mechanism, eqcids) # TODO check again for Fk , τk

    # restore old state
    for (i,id) in enumerate(bodyids)
        body = getbody(mechanism, id)
        body.state = statesold[i]
    end

    return A, Bu, Bλ, G
end


function lineardynamics(mechanism::Mechanism{T,N,Nb}, eqcids) where {T,N,Nb}
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    graph = mechanism.graph

    nc = 0
    for eqc in eqcs
        isinactive(eqc) && continue
        nc += length(eqc)
    end

    # calculate next state
    discretizestate!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    newton!(mechanism)

    # get state linearization 
    Fz = zeros(T,Nb*13,Nb*12)
    Fu = zeros(T,Nb*13,Nb*6)
    Fλ = zeros(T,Nb*13,nc)
    Ffz = zeros(T,Nb*13,Nb*13)
    invFfzquat = zeros(T,Nb*12,Nb*13)

    Bcontrol = zeros(T,Nb*6,length(eqcids))

    for (id,body) in pairs(bodies)
        col6 = offsetrange(id,6)
        col12 = offsetrange(id,12)
        col13 = offsetrange(id,13)

        Fzi = ∂F∂z(body, Δt)
        Fui = ∂F∂u(body, Δt)
        Ffzi, invFfzquati = ∂F∂fz(body, Δt)

        Fz[col13,col12] = Fzi
        Fu[col13,col6] = Fui
        Ffz[col13,col13] = Ffzi
        invFfzquat[col12,col13] = invFfzquati
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

    # display(round.(Ffz,digits=5))
    # display(round.(invFfzquat * inv(Ffz),digits=5))

    for (i,eqc) in enumerate(mechanism.eqconstraints)
        # display(eqc.λsol[2])

        K = zeros(T,9,9)
        K[1,1] = K[2,4] = K[3,7] = K[4,2] = K[5,5] = K[6,8] = K[7,3] = K[8,6] = K[9,9] = 1
        E = SMatrix{3,3,T,9}(I)

        parentid = eqc.parentid
        if parentid !== nothing
            body1 = getbody(mechanism, parentid)
            statea = body1.state
            pcol13 = offsetrange(parentid,13)

            for (consti, childid) in enumerate(eqc.childids)
                body2 = getbody(mechanism, childid)
                stateb = body2.state
                constraint = eqc.constraints[consti]
                ccol13 = offsetrange(childid,13)

                n1 = 1
                n2 = 0
                for i=1:consti-1
                    n1 += getN(eqc.constraints[i])
                    n2 += getN(eqc.constraints[i])
                end
                n2 += getN(eqc.constraints[consti])
                λ = eqc.λsol[2][n1:n2]

                Aaa = zeros(T,13,13)
                Aab = zeros(T,13,13)
                Aba = zeros(T,13,13)
                Abb = zeros(T,13,13)                

                XX, XQ, QX, QQ = ∂2g∂posaa(constraint, statea.xsol[2], statea.qsol[2], stateb.xsol[2], stateb.qsol[2]).*-1
                Aaa[4:6,1:3] = kron(λ'*reductionmat(constraint),E)*K*XX
                Aaa[4:6,7:10] = kron(λ'*reductionmat(constraint),E)*K*XQ
                Aaa[11:13,1:3] = kron(λ'*reductionmat(constraint),E)*K*QX
                Aaa[11:13,7:10] = kron(λ'*reductionmat(constraint),E)*K*QQ

                XX, XQ, QX, QQ = ∂2g∂posab(constraint, statea.xsol[2], statea.qsol[2], stateb.xsol[2], stateb.qsol[2]).*-1
                Aab[4:6,1:3] = kron(λ'*reductionmat(constraint),E)*K*XX
                Aab[4:6,7:10] = kron(λ'*reductionmat(constraint),E)*K*XQ
                Aab[11:13,1:3] = kron(λ'*reductionmat(constraint),E)*K*QX
                Aab[11:13,7:10] = kron(λ'*reductionmat(constraint),E)*K*QQ

                XX, XQ, QX, QQ = ∂2g∂posba(constraint, statea.xsol[2], statea.qsol[2], stateb.xsol[2], stateb.qsol[2]).*-1
                Aba[4:6,1:3] = kron(λ'*reductionmat(constraint),E)*K*XX
                Aba[4:6,7:10] = kron(λ'*reductionmat(constraint),E)*K*XQ
                Aba[11:13,1:3] = kron(λ'*reductionmat(constraint),E)*K*QX
                Aba[11:13,7:10] = kron(λ'*reductionmat(constraint),E)*K*QQ

                XX, XQ, QX, QQ = ∂2g∂posbb(constraint, statea.xsol[2], statea.qsol[2], stateb.xsol[2], stateb.qsol[2]).*-1
                Abb[4:6,1:3] = kron(λ'*reductionmat(constraint),E)*K*XX
                Abb[4:6,7:10] = kron(λ'*reductionmat(constraint),E)*K*XQ
                Abb[11:13,1:3] = kron(λ'*reductionmat(constraint),E)*K*QX
                Abb[11:13,7:10] = kron(λ'*reductionmat(constraint),E)*K*QQ

                Ffz[pcol13,pcol13] += Aaa
                Ffz[pcol13,ccol13] += Aab
                Ffz[ccol13,pcol13] += Aba
                Ffz[ccol13,ccol13] += Abb
            end
        else
            for (consti, childid) in enumerate(eqc.childids)
                body2 = getbody(mechanism, childid)
                stateb = body2.state
                constraint = eqc.constraints[consti]
                ccol13 = offsetrange(childid,13)

                n1 = 1
                n2 = 0
                for i=1:consti-1
                    n1 += getN(eqc.constraints[i])
                    n2 += getN(eqc.constraints[i])
                end
                n2 += getN(eqc.constraints[consti])
                λ = eqc.λsol[2][n1:n2]

                Abb = zeros(T,13,13)

                XX, XQ, QX, QQ = ∂2g∂posbb(constraint, stateb.xsol[2], stateb.qsol[2]).*-1
                Abb[4:6,1:3] = kron(λ'*reductionmat(constraint),E)*K*XX
                Abb[4:6,7:10] = kron(λ'*reductionmat(constraint),E)*K*XQ
                Abb[11:13,1:3] = kron(λ'*reductionmat(constraint),E)*K*QX
                Abb[11:13,7:10] = kron(λ'*reductionmat(constraint),E)*K*QQ

                Ffz[ccol13,ccol13] += Abb
            end
        end
    end

    G, Fλ = ∂g∂ʳextension(mechanism)

    invFfz = invFfzquat * inv(Ffz)
    A = - invFfz * Fz
    Bu = - invFfz * Fu * Bcontrol
    Bλ = - invFfz * Fλ

    # display(round.(Ffz,digits=5))
    # display(round.(invFfz,digits=5))
    # display(round.(Fz,digits=5))

    return A, Bu, Bλ, G
end


# function linearconstraints(mechanism::Mechanism{T,N,Nb}, eqcids) where {T,N,Nb}

# end