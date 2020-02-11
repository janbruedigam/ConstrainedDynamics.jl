@inline function setDandŝ!(diagonal::DiagonalEntry,body::Body,mechanism::Mechanism)
    diagonal.D = ∂dyn∂vel(body, mechanism.dt)
    diagonal.ŝ = dynamics(body, mechanism)
    return
end

@inline function extendDandŝ!(diagonal::DiagonalEntry,body::Body,c::InequalityConstraint,mechanism::Mechanism)
    dt = mechanism.dt
    diagonal.D += diagval(c,body,dt)
    diagonal.ŝ += g(c,body,mechanism)
    return
end

@inline function setDandŝ!(d::DiagonalEntry{T,N},c::EqualityConstraint,mechanism::Mechanism) where {T,N}
    d.D = @SMatrix zeros(T,N,N)
    # μ = 1e-05
    # d.D = SMatrix{N,N,T,N*N}(μ*I)
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

@inline update!(body::Body,ldu::SparseLDU,αsmax) = update!(body,getentry(ldu,body.id),αsmax)
function update!(body::Body,diagonal::DiagonalEntry,αsmax)
    body.s1 = body.s0 - αsmax*diagonal.ŝ
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
    return
end

@inline function s1tos0!(ineq::InequalityConstraint)
    ineq.sl0 = ineq.sl1
    ineq.ga0 = ineq.ga1
    return
end

@inline function normΔs(component::Component)
    difference = component.s1-component.s0
    return dot(difference,difference)
end


function SLGASol!(ineqentry::InequalityEntry,diagonal::DiagonalEntry,body::Body,ineq::InequalityConstraint,mechanism::Mechanism)
    dt = mechanism.dt
    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    γ = ineq.ga1
    s = ineq.sl1
    Σ = γ/s
    Σm = s/γ
    μ = mechanism.μ
    φ = body.x[2][3]+dt*body.s1[3]

    Δv = diagonal.ŝ

    # SHOULD BE SWITCHED ???
    ineqentry.sl = Σm*(γ - ineqentry.ga) - μ/γ
    ineqentry.ga = Σ*(φ - Nv*Δv) - μ/s


    return
end
