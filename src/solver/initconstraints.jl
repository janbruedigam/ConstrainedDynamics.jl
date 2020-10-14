

function get_free_bodies(mechanism::Mechanism{T},fixedbodies) where T
    freebodies = Body{T}[]
    for (j,body) in enumerate(mechanism.bodies)
        if !(body.id in fixedbodies)
        push!(freebodies,body)
        end  
    end  
    return UnitDict(freebodies) 
end
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
function ∂g∂posc(mechanism::Mechanism{T,N,Nb},fixedbodies) where {T,N,Nb}
    rangeDict = Dict{Int64,UnitRange}()
    ind1 = 1
    ind2 = 0

    for (i,eqc) in enumerate(mechanism.eqconstraints)
        ind2 += length(eqc)
        range = ind1:ind2
        rangeDict[i] = range

        ind1 = ind2+1
    end
    Nfb=length(fixedbodies)
    G = zeros(ind2,(Nb-Nfb)*7)
    
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        for (j,body) in enumerate(get_free_bodies(mechanism,fixedbodies))
             G[rangeDict[i],offsetrange(j,7)] = ∂g∂posc(mechanism, eqc, body.id)
        end
    end
        
    return G
end


function constraintstep!(mechanism::Mechanism{T,N,Nb},fixedbodies) where {T,N,Nb}
    
    gval=gc(mechanism)
    pinv∂g∂posc=pinv(∂g∂posc(mechanism,fixedbodies))
    stepvec = -pinv∂g∂posc*gval

    for (i,body) in enumerate(get_free_bodies(mechanism,fixedbodies))

            body.state.vsol[1] = stepvec[offsetrange(i, 3, 7, 1)]
            # Limit quaternion step to feasible length
            range = offsetrange(i, 3, 7, 2)
            Δstemp = VLᵀmat(body.state.qk[1]) * stepvec[first(range):(last(range)+1)]
            if norm(Δstemp) > 1
                Δstemp = Δstemp/norm(Δstemp)
            end
            body.state.ωsol[1] = Δstemp
        
    end

    return
end

function initializeConstraints!(mechanism::Mechanism{T,N,Nb,Ne},fixedbodies; ε = 1e-5, newtonIter = 100, lineIter = 10) where {T,N,Nb,Ne}
    
    norm0 = norm(gc(mechanism))
    norm1 = norm0
    for n = Base.OneTo(newtonIter)

        for body in get_free_bodies(mechanism,fixedbodies)  
            body.state.xk[1] = body.state.xc
            body.state.qk[1] = body.state.qc
        end

        constraintstep!(mechanism,fixedbodies) 
    
        for j = Base.OneTo(lineIter)

            for body in get_free_bodies(mechanism,fixedbodies)
                
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
