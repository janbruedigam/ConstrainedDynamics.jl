@inline function setDandΔs!(mechanism::Mechanism, diagonal::DiagonalEntry, body::Body)
    diagonal.D = ∂dyn∂vel(body, mechanism.Δt)
    diagonal.Δs = dynamics(mechanism, body)
    return
end

@inline function extendDandΔs!(mechanism::Mechanism, diagonal::DiagonalEntry, body::Body, ineqc::InequalityConstraint)
    diagonal.D += schurD(ineqc, body, mechanism.Δt) # + SMatrix{6,6,Float64,36}(1e-5*I)
    diagonal.Δs += schurf(mechanism, ineqc, body)
    return
end

@inline function setDandΔs!(mechanism::Mechanism, diagonal::DiagonalEntry{T,N}, eqc::EqualityConstraint) where {T,N}
    diagonal.D = @SMatrix zeros(T, N, N)
    # μ = 1e-5
    # diagonal.D = SMatrix{N,N,T,N*N}(μ*I) # TODO Positiv because of weird system? fix generally
    diagonal.Δs = g(mechanism, eqc)
    return
end

@inline function setLU!(mechanism::Mechanism, offdiagonal::OffDiagonalEntry, bodyid::Integer, eqc::EqualityConstraint)
    offdiagonal.L = ∂g∂pos(mechanism, eqc, bodyid)'
    offdiagonal.U = ∂g∂vel(mechanism, eqc, bodyid)
    return
end

@inline function setLU!(mechanism::Mechanism, offdiagonal::OffDiagonalEntry, eqc::EqualityConstraint, bodyid::Integer)
    offdiagonal.L = ∂g∂vel(mechanism, eqc, bodyid)
    offdiagonal.U = ∂g∂pos(mechanism, eqc, bodyid)'
    return
end

@inline function setLU!(offdiagonal::OffDiagonalEntry{T,N1,N2}) where {T,N1,N2}
    offdiagonal.L = @SMatrix zeros(T, N2, N1)
    offdiagonal.U = offdiagonal.L'
    return
end

@inline function updateLU1!(offdiagonal::OffDiagonalEntry, diagonal::DiagonalEntry, maintograndchild::OffDiagonalEntry, childtograndchild::OffDiagonalEntry)
    D = diagonal.D
    offdiagonal.L -= maintograndchild.L * D * childtograndchild.U
    offdiagonal.U -= childtograndchild.L * D * maintograndchild.U
    return
end

@inline function updateLU2!(offdiagonal::OffDiagonalEntry, diagonal::DiagonalEntry)
    Dinv = diagonal.Dinv
    offdiagonal.L = offdiagonal.L * Dinv
    offdiagonal.U = Dinv * offdiagonal.U
    return
end

@inline function updateD!(diagonal::DiagonalEntry, childdiagonal::DiagonalEntry, fillin::OffDiagonalEntry)
    diagonal.D -= fillin.L * childdiagonal.D * fillin.U
    return
end

function invertD!(diagonal::DiagonalEntry)
    diagonal.Dinv = inv(diagonal.D)
    return
end

@inline function LSol!(diagonal::DiagonalEntry, child::DiagonalEntry, fillin::OffDiagonalEntry)
    diagonal.Δs -= fillin.L * child.Δs
    return
end

function DSol!(diagonal::DiagonalEntry)
    diagonal.Δs = diagonal.Dinv * diagonal.Δs
    return
end

@inline function USol!(diagonal::DiagonalEntry, parent::DiagonalEntry, fillin::OffDiagonalEntry)
    diagonal.Δs -= fillin.U * parent.Δs
    return
end


function feasibilityStepLength!(mechanism::Mechanism)
    ldu = mechanism.ldu

    τ = 0.995
    mechanism.α = 1.

    for ineqc in mechanism.ineqconstraints
        feasibilityStepLength!(mechanism, ineqc, getineqentry(ldu, ineqc.id), τ)
    end

    return
end

function feasibilityStepLength!(mechanism, ineqc::InequalityConstraint{T,N}, ineqentry::InequalityEntry, τ) where {T,N}
    s1 = ineqc.ssol[2]
    γ1 = ineqc.γsol[2]
    Δs = ineqentry.Δs
    Δγ = ineqentry.Δγ

    for i = 1:N
        αmax = τ * s1[i] / Δs[i]
        (αmax > 0) && (αmax < mechanism.α) && (mechanism.α = αmax)
        αmax = τ * γ1[i] / Δγ[i]
        (αmax > 0) && (αmax < mechanism.α) && (mechanism.α = αmax)
    end
    return
end


function setentries!(mechanism::Mechanism)
    graph = mechanism.graph
    ldu = mechanism.ldu

    for (id, body) in pairs(mechanism.bodies)
        isinactive(graph, id) && continue
        
        for cid in directchildren(graph, id)
            setLU!(mechanism, getentry(ldu, (id, cid)), id, geteqconstraint(mechanism, cid))
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(mechanism, diagonal, body)
        for cid in ineqchildren(graph, id)
            ineqc = getineqconstraint(mechanism, cid)
            calcFrictionForce!(mechanism, ineqc)
            extendDandΔs!(mechanism, diagonal, body, ineqc)
        end
    end

    for eqc in mechanism.eqconstraints
        id = eqc.id
        isinactive(graph, id) && continue
        
        for cid in directchildren(graph, id)
            setLU!(mechanism, getentry(ldu, (id, cid)), eqc, cid)
        end

        for cid in loopchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)))
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(mechanism, diagonal, eqc)
    end

    return 
end

function factor!(graph::Graph, ldu::SparseLDU)
    for id in graph.dfslist
        isinactive(graph, id) && continue

        sucs = successors(graph, id)
        for childid in sucs
            isinactive(graph, childid) && continue
            
            offdiagonal = getentry(ldu, (id, childid))
            for grandchildid in sucs
                grandchildid == childid && break
                isinactive(graph, grandchildid) && continue

                if hasdirectchild(graph, childid, grandchildid)
                    updateLU1!(offdiagonal, getentry(ldu, grandchildid), getentry(ldu, (id, grandchildid)), getentry(ldu, (childid, grandchildid)))
                end
            end
            updateLU2!(offdiagonal, getentry(ldu, childid))
        end

        diagonal = getentry(ldu, id)

        for childid in successors(graph, id)
            isinactive(graph, childid) && continue
            updateD!(diagonal, getentry(ldu, childid), getentry(ldu, (id, childid)))
        end
        invertD!(diagonal)
    end

    return 
end

function solve!(mechanism)
    ldu = mechanism.ldu
    graph = mechanism.graph
    dfslist = graph.dfslist

    for id in dfslist
        isinactive(graph, id) && continue

        diagonal = getentry(ldu, id)

        for childid in successors(graph, id)
            isinactive(graph, childid) && continue
            LSol!(diagonal, getentry(ldu, childid), getentry(ldu, (id, childid)))
        end
    end

    for id in graph.rdfslist
        isinactive(graph, id) && continue

        diagonal = getentry(ldu, id)
        DSol!(diagonal)

        for parentid in predecessors(graph, id)
            isinactive(graph, parentid) && continue
            USol!(diagonal, getentry(ldu, parentid), getentry(ldu, (parentid, id)))
        end

        for childid in ineqchildren(graph, id)
            isinactive(graph, childid) && continue
            eliminatedsolve!(mechanism, getineqentry(ldu, childid), diagonal, getbody(mechanism, id), getineqconstraint(mechanism, childid))
        end
    end

    return 
end

function eliminatedsolve!(mechanism::Mechanism, ineqentry::InequalityEntry, diagonal::DiagonalEntry, body::Body, ineqc::InequalityConstraint)
    Δt = mechanism.Δt
    μ = mechanism.μ

    φ = g(mechanism, ineqc)

    Nx = ∂g∂pos(mechanism, ineqc, body)
    Nv = ∂g∂vel(mechanism, ineqc, body)

    γ1 = ineqc.γsol[2]
    s1 = ineqc.ssol[2]

    Δv = diagonal.Δs
    ineqentry.Δγ = -γ1 ./ s1 .* φ + μ ./ s1 + γ1 ./ s1 .* (Nv * Δv)
    ineqentry.Δs = s1 .- μ ./ γ1 + s1 ./ γ1 .* ineqentry.Δγ

    return
end


@inline function updatesolution!(body::Body)
    body.state.vsol[1] = body.state.vsol[2]
    body.state.ωsol[1] = body.state.ωsol[2]
    return
end

@inline function updatesolution!(eqc::EqualityConstraint)
    eqc.λsol[1] = eqc.λsol[2]
    return
end

@inline function updatesolution!(ineqc::InequalityConstraint)
    ineqc.ssol[1] = ineqc.ssol[2]
    ineqc.γsol[1] = ineqc.γsol[2]
    return
end

@inline function normΔs(body::Body)
    d1 = body.state.vsol[2] - body.state.vsol[1]
    d2 = body.state.ωsol[2] - body.state.ωsol[1]
    return dot(d1, d1) + dot(d2, d2)
end

@inline function normΔs(eqc::EqualityConstraint)
    d = eqc.λsol[2] - eqc.λsol[1]
    return dot(d, d)
end

@inline function normΔs(ineqc::InequalityConstraint)
    d1 = ineqc.ssol[2] - ineqc.ssol[1]
    d2 = ineqc.γsol[2] - ineqc.γsol[1]
    return dot(d1, d1) + dot(d2, d2)
end


