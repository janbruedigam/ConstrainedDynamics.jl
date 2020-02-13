@inline function setDandŝ!(diagonal::DiagonalEntry,body::Body,mechanism::Mechanism)
    diagonal.D = ∂dyn∂vel(body, mechanism.dt)
    diagonal.ŝ = dynamics(body, mechanism)
    return
end

@inline function extendDandŝ!(diagonal::DiagonalEntry,body::Body,c::InequalityConstraint,mechanism::Mechanism)
    dt = mechanism.dt
    diagonal.D += diagval(c,body,dt) #+ SMatrix{6,6,Float64,36}(1e-5*I)
    diagonal.ŝ += dynineq(c,body,mechanism)
    return
end

@inline function setDandŝ!(d::DiagonalEntry{T,N},c::EqualityConstraint,mechanism::Mechanism) where {T,N}
    # d.D = @SMatrix zeros(T,N,N)
    μ = 1e-5
    d.D = SMatrix{N,N,T,N*N}(μ*I) # TODO Positiv because of weird system? fix generally
    d.ŝ = g(c,mechanism)
    return
end

@inline function setLU!(o::OffDiagonalEntry,bodyid::Int64,c::EqualityConstraint,mechanism)
    o.L = -∂g∂pos(c,bodyid,mechanism)'
    o.U = ∂g∂vel(c,bodyid,mechanism)
    return
end

@inline function setLU!(o::OffDiagonalEntry,c::EqualityConstraint,bodyid::Int64,mechanism)
    o.L = ∂g∂vel(c,bodyid,mechanism)
    o.U = -∂g∂pos(c,bodyid,mechanism)'
    return
end

@inline function setLU!(o::OffDiagonalEntry{T,N1,N2}) where {T,N1,N2}
    o.L = @SMatrix zeros(T,N2,N1)
    o.U = o.L'
    return
end

@inline function updateLU1!(o::OffDiagonalEntry,d::DiagonalEntry,gc::OffDiagonalEntry,cgc::OffDiagonalEntry)
    D = d.D
    o.L -= gc.L*D*cgc.U
    o.U -= cgc.L*D*gc.U
    return
end

@inline function updateLU2!(o::OffDiagonalEntry,d::DiagonalEntry)
    Dinv = d.Dinv
    o.L = o.L*Dinv
    o.U = Dinv*o.U
    return
end

@inline function updateD!(d::DiagonalEntry,c::DiagonalEntry,f::OffDiagonalEntry)
    d.D -= f.L*c.D*f.U
    return
end

function invertD!(d::DiagonalEntry)
    d.Dinv = inv(d.D)
    return
end

@inline function LSol!(d::DiagonalEntry,child::DiagonalEntry,fillin::OffDiagonalEntry)
    d.ŝ -= fillin.L*child.ŝ
    return
end

function DSol!(d::DiagonalEntry)
    d.ŝ = d.Dinv*d.ŝ
    return
end

@inline function USol!(d::DiagonalEntry,parent::DiagonalEntry,fillin::OffDiagonalEntry)
    d.ŝ -= fillin.U*parent.ŝ
    return
end


function factor!(graph::Graph,ldu::SparseLDU)
    for id in graph.dfslist
        sucs = successors(graph,id)
        for cid in sucs
            offdiagonal = getentry(ldu,(id,cid))
            for gcid in sucs
                gcid == cid && break
                if hasdirectchild(graph,cid,gcid)
                    updateLU1!(offdiagonal,getentry(ldu,gcid),getentry(ldu,(id,gcid)),getentry(ldu,(cid,gcid)))
                end
            end
            updateLU2!(offdiagonal,getentry(ldu,cid))
        end

        diagonal = getentry(ldu,id)

        for cid in successors(graph,id)
            updateD!(diagonal,getentry(ldu,cid),getentry(ldu,(id,cid)))
        end
        invertD!(diagonal)
    end
end

function solve!(graph::Graph,ldu::SparseLDU,mechanism)
    dfslist = graph.dfslist

    for id in dfslist
        diagonal = getentry(ldu,id)

        for cid in successors(graph,id)
            LSol!(diagonal,getentry(ldu,cid),getentry(ldu,(id,cid)))
        end
    end

    for id in graph.rdfslist
        diagonal = getentry(ldu,id)

        DSol!(diagonal)

        for pid in predecessors(graph,id)
            USol!(diagonal,getentry(ldu,pid),getentry(ldu,(pid,id)))
        end

        # inequalities
        for cid in ineqchildren(graph,id)
            SLGASol!(getineq(ldu,cid),diagonal,getbody(mechanism,id),getineqconstraint(mechanism,cid),mechanism)
        end
    end
end

@inline update!(component::Component,ldu::SparseLDU) = update!(component,getentry(ldu,component.id))
function update!(component::Component,diagonal::DiagonalEntry)
    component.s1 = component.s0 - diagonal.ŝ
    return
end

@inline update!(eq::EqualityConstraint,ldu::SparseLDU,αγmax) = update!(eq,getentry(ldu,eq.id),αγmax)
function update!(eq::EqualityConstraint,diagonal::DiagonalEntry,αγmax)
    eq.s1 = eq.s0 - αγmax*diagonal.ŝ
    return
end

@inline update!(ineq::InequalityConstraint,ldu::SparseLDU,αsmax,αγmax) = update!(ineq,getineq(ldu,ineq.id),αsmax,αγmax)
function update!(ineq::InequalityConstraint,entry::InequalityEntry,αsmax,αγmax)
    ineq.sl1 = ineq.sl0 - αsmax*entry.sl
    ineq.ga1 = ineq.ga0 - αγmax*entry.ga
    ineq.slf1 = ineq.slf0 - αsmax*entry.slf
    ineq.psi1 = ineq.psi0 - αγmax*entry.psi
    ineq.b1 = ineq.b0 - αsmax*entry.b
    return
end

@inline function s0tos1!(component::Component)
    component.s1 = component.s0
    return
end

@inline function s1tos0!(component::Component)
    component.s0 = component.s1
    return
end

@inline function s0tos1!(ineq::InequalityConstraint)
    ineq.sl1 = ineq.sl0
    ineq.ga1 = ineq.ga0
    ineq.slf1 = ineq.slf0
    ineq.psi1 = ineq.psi0
    ineq.b1 = ineq.b0
    return
end

@inline function s1tos0!(ineq::InequalityConstraint)
    ineq.sl0 = ineq.sl1
    ineq.ga0 = ineq.ga1
    ineq.slf0 = ineq.slf1
    ineq.psi0 = ineq.psi1
    ineq.b0 = ineq.b1
    return
end

@inline function normΔs(component::Component)
    difference = component.s1-component.s0
    return dot(difference,difference)
end

@inline function normΔs(ineq::InequalityConstraint)
    difference = [ineq.sl1-ineq.sl0;ineq.ga1-ineq.ga0]
    return dot(difference,difference)
end


function SLGASol!(ineqentry::InequalityEntry,diagonal::DiagonalEntry,body::Body,ineq::InequalityConstraint,mechanism::Mechanism)
    dt = mechanism.dt
    μ = mechanism.μ



    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = [1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    ga1 = ineq.ga1
    sl1 = ineq.sl1
    slf1 = ineq.slf1
    psi1 = ineq.psi1
    b1 = ineq.b1

    cf = 0.2

    Ax = [-Nx-b1'*D;zeros(6)']
    Av = [Nv;zeros(6)']
    H = [D*s1+2*psi1*b1 2*ga1*b1]
    X = [0 0;2*cf^2*ga1 0]
    B = [zeros(2)';b1']
    Z = [ga1 0;0 psi1]
    S = [sl1 0;0 slf1]

    K = X + 1/2*1/(ga1*psi1)*B*H + Z\S
    φ = body.x[2][3]+dt*body.s1[3]
    ci = [φ;cf^2*ga1^2-b1'*b1]

    Δv = diagonal.ŝ
    temp1 = K\(ci+1/2*1/psi1*B*D*s1+B*b1-μ*inv(Z)*ones(2) - (Av+1/2*1/psi1*B*D)*Δv)
    ineqentry.ga = temp1[1]
    ineqentry.psi = temp1[2]
    temp2 = S*ones(2)-μ*inv(Z)*ones(2)-Z\S*[ineqentry.ga;ineqentry.psi]
    ineqentry.sl = temp2[1]
    ineqentry.slf = temp2[2]
    ineqentry.b = 1/2*1/psi1*(D*s1 + 2*psi1*b1-D*Δv-1/ga1*H*temp1)

    return
end
