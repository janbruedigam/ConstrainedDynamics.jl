function gc(mechanism::Mechanism{T}) where T
    rangeDict = Dict{Int64,UnitRange}()
    ind1 = 1
    ind2 = 0

    for (i,eqc) in enumerate(mechanism.eqconstraints)
        
        ind2 += length(eqc)
        range = ind1:ind2
        rangeDict[i] = range
        ind1 = ind2+1
    end

    gval = zeros(T,ind2)
    
    for (i,eqc) in enumerate(mechanism.eqconstraints)

        gval[rangeDict[i]] = gc(mechanism, eqc)
        
    end 
    return gval
end   

# Derivatives
function ∂g∂posc(mechanism::Mechanism{T,N,Nb}) where {T,N,Nb}
    rangeDict = Dict{Int64,UnitRange}()
    ind1 = 1
    ind2 = 0

    for (i,eqc) in enumerate(mechanism.eqconstraints)
        ind2 += length(eqc)
        range = ind1:ind2
        rangeDict[i] = range

        ind1 = ind2+1
    end

    G = zeros(ind2,Nb*7)
    
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        for (j,body) in enumerate(mechanism.bodies)
            G[rangeDict[i],offsetrange(j,7)] = ∂g∂posc(mechanism, eqc, body.id)
        end
    end
        
    return G
end


function constraintstep!(mechanism::Mechanism{T,N,Nb}) where {T,N,Nb}
    bodies = mechanism.bodies

    gval=gc(mechanism)
    pinv∂g∂posc=pinv(∂g∂posc(mechanism))
    stepvec = -pinv∂g∂posc*gval

    for (i,body) in enumerate(mechanism.bodies)
        body.state.vsol[1] = stepvec[offsetrange(i, 3, 7, 1)]

        # Limit quaternion step to feasible length
        range = offsetrange(i, 3, 7, 2)
        Δstemp = VLᵀmat(bodies[i].state.qk[1]) * stepvec[first(range):(last(range)+1)]
        if norm(Δstemp) > 1
            Δstemp = Δstemp/norm(Δstemp)
        end
        body.state.ωsol[1] = Δstemp
    end

    return
end

function initializeConstraints!(mechanism::Mechanism{T,N,Nb,Ne}; ε = 1e-5, newtonIter = 100, lineIter = 10) where {T,N,Nb,Ne}
    bodies = mechanism.bodies

    norm0 = norm(gc(mechanism))
    norm1 = norm0
    for n = Base.OneTo(newtonIter)
        for body in bodies  
            body.state.xk[1] = body.state.xc
            body.state.qk[1] = body.state.qc
        end

        constraintstep!(mechanism) 

        for j = Base.OneTo(newtonIter)
            for body in bodies
                body.state.xc = body.state.xk[1] + body.state.vsol[1]/(2^(j-1))
                w = sqrt(1-norm(body.state.ωsol[1]/(2^(j-1)))^2)
                body.state.qc = body.state.qk[1] * UnitQuaternion(w,body.state.ωsol[1]/(2^(j-1))...,false)
            end

            norm1 = norm(gc(mechanism))

            if norm1 < norm0 
                break
            end
        end

        if norm1 < ε
            return
        else
            norm0 = norm1
        end
    end

    display("Constraint initialization did not converge!")
    return 
end
