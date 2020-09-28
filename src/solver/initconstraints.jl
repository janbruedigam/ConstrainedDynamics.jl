using ConstrainedDynamics:offsetrange 

function Liniensuche(xold,qold,sx,sq,j)
    
    sqtemp = sq/(2^(j-1)) 
    sxtemp = sx/(2^(j-1))
    w = sqrt(1-norm(sqtemp)^2) 
    sqexpanded = [w;sqtemp]
    qnew=Lmat(UnitQuaternion(qold...,false))*sqexpanded 
    xnew=xold+sxtemp 
    return xnew,qnew    

end

# Constraint functions
function g(mech)
    G=Float64[]
    for eqc in mech.eqconstraints
        gi=gc(mech,eqc)
        for z in gi
            push!(G,z)
        end
    end 
    return G
end   

# Derivatives
function derivative(mech::Mechanism{T,N,Nb,Ne,Ni},Nc) where {T,N,Nb,Ne,Ni}
    D=zeros(Nc,Nb*7)
    j=1
    jd=0
    for eqc in mech.eqconstraints
        for bd in mech.bodies
            if (bd.id==eqc.parentid)||(bd.id in eqc.childids)
                d=∂g∂posc(mech, eqc, bd.id)
                jd=size(d,1)
                D[j:j+jd-1,offsetrange(bd.id,7)]=d
            end
        end
        j=j+jd
    end
        return D
end



function initializeConstraints!(mech::Mechanism{T,N,Nb,Ne,Ni},epsilon,Nmax) where {T,N,Nb,Ne,Ni}
    
    Xold=zeros(3*Nb,1)
    qold=zeros(4*Nb,1)
    conv=false

    for k=1:Nmax 
        for bd in mech.bodies  

            Xold[offsetrange(bd.id,3)]=bd.state.xc
            qold[offsetrange(bd.id,4)]=params(bd.state.qc)

        end
        #println("[xold,qold] ",[Xold;qold])

        if (norm(g(mech)) <= epsilon)
            conv=true;
            break
        end

        #step calculation 
        M=zeros(6*Nb,7*Nb)
        for bd in mech.bodies 

            Mi=[Matrix{Float64}(I, 3, 3) zeros(3,4); zeros(3,3) VLᵀmat(UnitQuaternion(qold[offsetrange(bd.id,4)]...,false))]
            M[offsetrange(bd.id,6),offsetrange(bd.id,7)]=Mi

        end 

        G=g(mech)
        n=size(G,1)
        invderiv=pinv(derivative(mech,n))
        Δs=-M*invderiv*G 
        for l=1:Nb
            if norm(Δs[6*l-2:6*l])>1
                Δs[6*l-2:6*l] = Δs[6*l-2:6*l]/norm(Δs[6*l-2:6*l])
            end
        end

        for j=1:Nmax
            normold=norm(g(mech))
            for bd in mech.bodies 

                sxi=Δs[6*bd.id-5:6*bd.id-3]
                sqi=Δs[6*bd.id-2:6*bd.id]
                xnewi,qnewi = Liniensuche(Xold[offsetrange(bd.id,3)],qold[offsetrange(bd.id,4)],sxi,sqi,j)
                #println("[xnew,qnew] ",[xnewi;qnewi])
                setPosition!(bd,x=xnewi,q=UnitQuaternion(qnewi...,false))

            end

            if norm(g(mech)) < normold 
                break
            end
        end

    
    end

    println("[xfinal,qfinal] ",[Xold;qold])
    return conv
    
end

    