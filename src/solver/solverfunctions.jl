@inline function setDandΔs!(diagonal::DiagonalEntry,body::Body,mechanism::Mechanism)
    diagonal.D = ∂dyn∂vel(body, mechanism.dt)
    diagonal.Δs = dynamics(body, mechanism)
    return
end

@inline function extendDandΔs!(diagonal::DiagonalEntry,body::Body,c::InequalityConstraint,mechanism::Mechanism)
    dt = mechanism.dt
    diagonal.D += schurD(c,body,dt) #+ SMatrix{6,6,Float64,36}(1e-5*I)
    diagonal.Δs += schurf(c,body,mechanism)
    return
end

@inline function setDandΔs!(d::DiagonalEntry{T,N},c::EqualityConstraint,mechanism::Mechanism) where {T,N}
    d.D = @SMatrix zeros(T,N,N)
    # μ = 1e-5
    # d.D = SMatrix{N,N,T,N*N}(μ*I) # TODO Positiv because of weird system? fix generally
    d.Δs = g(c,mechanism)
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
    d.Δs -= fillin.L*child.Δs
    return
end

function DSol!(d::DiagonalEntry)
    d.Δs = d.Dinv*d.Δs
    return
end

@inline function USol!(d::DiagonalEntry,parent::DiagonalEntry,fillin::OffDiagonalEntry)
    d.Δs -= fillin.U*parent.Δs
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
            eliminatedSol!(getineq(ldu,cid),diagonal,getbody(mechanism,id),getineqconstraint(mechanism,cid),mechanism)
        end
    end
end

@inline update!(component::Component,ldu::SparseLDU) = update!(component,getentry(ldu,component.id))
function update!(component::Component,diagonal::DiagonalEntry)
    component.s1 = component.s0 - diagonal.Δs
    return
end

@inline update!(eq::EqualityConstraint,ldu::SparseLDU,αγmax) = update!(eq,getentry(ldu,eq.id),αγmax)
function update!(eq::EqualityConstraint,diagonal::DiagonalEntry,αγmax)
    eq.s1 = eq.s0 - αγmax*diagonal.Δs
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
    ineq.s1 = ineq.s0
    ineq.γ1 = ineq.γ0
    return
end

@inline function s1tos0!(ineq::InequalityConstraint)
    ineq.s0 = ineq.s1
    ineq.γ0 = ineq.γ1
    return
end

@inline function normΔs(component::Component)
    difference = component.s1-component.s0
    return dot(difference,difference)
end

@inline function normΔs(ineq::InequalityConstraint)
    difference = [ineq.s1-ineq.s0;ineq.γ1-ineq.γ0]
    return dot(difference,difference)
end


# function eliminatedSol!(ineqentry::InequalityEntry,diagonal::DiagonalEntry,body::Body,ineq::InequalityConstraint,mechanism::Mechanism)
#     dt = mechanism.dt
#     μ = mechanism.μ
#     No = 2
#
#     φ = body.x[No][3]+dt*body.s1[3]
#     cf = ineq.constraints.cf
#
#     Nx = SVector{6,Float64}(0,0,1,0,0,0)'
#     Nv = dt*Nx
#     D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]
#
#     s1 = body.s1
#     γ1 = ineq.γ1[1]
#     sl1 = ineq.s1[1]
#
#     Δv = diagonal.Δs
#     ineqentry.Δγ = [γ1/sl1*φ - μ/sl1 - γ1/sl1*Nv*Δv]
#     ineqentry.Δs = [sl1 - μ/γ1 - sl1/γ1*ineqentry.Δγ[1]]
#
#     return
# end

# function eliminatedSol!(ineqentry::InequalityEntry,diagonal::DiagonalEntry,body::Body,ineq::InequalityConstraint,mechanism::Mechanism)
#     dt = mechanism.dt
#     μ = mechanism.μ
#     No = 2
#
#     impact = ineq.constraints[1]
#
#     φ = g(impact,body,dt,No)
#
#     Nx = impact.Nx
#     Nv = dt*Nx
#
#     γ1 = ineq.γ1[1]
#     s1 = ineq.s1[1]
#
#     Δv = diagonal.Δs
#     ineqentry.Δγ = [γ1/s1*φ - μ/s1 - γ1/s1*Nv*Δv]
#     ineqentry.Δs = [s1 - μ/γ1 - s1/γ1*ineqentry.Δγ[1]]
#
#     return
# end

function eliminatedSol!(ineqentry::InequalityEntry,diagonal::DiagonalEntry,body::Body,ineqc::InequalityConstraint,mechanism::Mechanism)
    dt = mechanism.dt
    μ = mechanism.μ
    No = 2

    φ = g(ineqc,mechanism)

    Nx = ∂g∂pos(ineqc,body,mechanism)
    Nv = ∂g∂vel(ineqc,body,mechanism)

    γ1 = ineqc.γ1
    s1 = ineqc.s1

    Δv = diagonal.Δs
    ineqentry.Δγ = γ1./s1.*φ - μ./s1 - γ1./s1.*(Nv*Δv)
    ineqentry.Δs = s1 .- μ./γ1 - s1./γ1.*ineqentry.Δγ

    return
end

# function eliminateds(ineqentry::InequalityEntry,ineqc::InequalityConstraint,φ,μ,γ,Nv,Δγ)
#     γ1 = ineqc.γ1
#     s1 = ineqc.s1
#     ineqentry.Δs = s1 .- μ./γ1 - s1./γ1.*Δγ
#     return
# end

# function eliminatedγ(ineqentry::InequalityEntry,ineqc::InequalityConstraint,φ,μ,γ,Nv,Δv)
#     γ1 = ineqc.γ1
#     s1 = ineqc.s1
#     ineqentry.Δγ = γ1./s1.*φ - μ./s1 - γ1./s1.*(Nv*Δv)
#     return
# end