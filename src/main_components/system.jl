function create_system(origin::Origin{T}, bodies::Vector{<:Body},
        eqconstraints::Vector{<:EqualityConstraint}, ineqconstraints::Vector{<:InequalityConstraint}
    ) where T

    adjacency, ids, dims = adjacencyMatrix(eqconstraints, bodies)
    oid = origin.id
    graphoid = findfirst(x->x==oid, ids)
    adjacency = deleteat(adjacency, graphoid)
    dims = deleteat!(dims, graphoid)
    ids = deleteat!(ids, graphoid)

    system = System{T}(adjacency, dims; ids = ids)

    for eqc in eqconstraints
        eqc.parentid == oid && (eqc.parentid = nothing)
    end

    return system
end

function adjacencyMatrix(eqconstraints::Vector{<:EqualityConstraint}, bodies::Vector{<:Body})
    A = zeros(Bool, 0, 0)
    ids = Int64[]
    dims = zeros(Int64, 0)
    n = 0

    for constraint in eqconstraints
        constraintid = constraint.id
        A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
        n += 1
        append!(ids, constraintid)
        append!(dims,length(constraint))
        for bodyid in unique([constraint.parentid;constraint.childids])
            if bodyid ∉ ids
                A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
                n += 1
                append!(ids, bodyid)
                append!(dims,6)
            end
            A[findfirst(x->x==constraintid, ids),findfirst(x->x==bodyid, ids)] = true
        end
    end
    for body in bodies # add unconnected bodies
        if body.id ∉ ids
            A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
            n += 1
            append!(ids, body.id)
            append!(dims,6)
        end
    end

    A = A .| A'
    return convert(Matrix{Int64}, A), ids, dims
end

@inline getentry(system, id1, id2) = system.matrix_entries[(id1, id2)]
@inline getentry(system, id) = system.vector_entries[id]

function recursivedirectchildren!(system, id::Integer)
    dirs = copy(children(system, id))
    dirslocal = copy(dirs)
    for childid in dirslocal
        append!(dirs, recursivedirectchildren!(system, childid))
    end
    return dirs
end


function densesystem(mechanism::Mechanism{T,Nn,Nb}) where {T,Nn,Nb}
    eqcs = mechanism.eqconstraints
    system = mechanism.system
    system = mechanism.system

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

    for id in system.dfs_list
        component = getcomponent(mechanism, id)
        ind2 += length(component)
        range = ind1:ind2
        rangeDict[id] = range


        # A
        diagonal = getentry(system,id,id)
        A[range,range] = diagonal.value

        for childid in system.acyclic_children[id]
            offdiagonal_L = getentry(system, id, childid)
            offdiagonal_U = getentry(system, childid, id)
            nc1 = first(rangeDict[childid])
            nc2 = last(rangeDict[childid])

            A[range,nc1:nc2] = offdiagonal_L.value
            A[nc1:nc2,range] = offdiagonal_U.value
        end

        for cyclic_children in system.cycles[id]
            for childid in cyclic_children
                offdiagonal_L = getentry(system, id, childid)
                offdiagonal_U = getentry(system, childid, id)
                nc1 = first(rangeDict[childid])
                nc2 = last(rangeDict[childid])

                A[range,nc1:nc2] = offdiagonal_L.value
                A[nc1:nc2,range] = offdiagonal_U.value
            end
        end

        # x
        sol = getentry(system,id)
        x[range] = sol.value

        # b
        if component isa Body
            b[range] = -dynamics(mechanism, component)
        else
            b[range] = -g(mechanism, component)
        end


        ind1 = ind2+1
    end    
    
    return A, x, b
end