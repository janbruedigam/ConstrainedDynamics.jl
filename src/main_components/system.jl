function create_system(origin::Origin{T}, bodies::Vector{<:Body},
        eqconstraints::Vector{<:EqualityConstraint}, ineqconstraints::Vector{<:InequalityConstraint}
    ) where T

    adjacency, dims = adjacencyMatrix(eqconstraints, bodies, ineqconstraints)
    system = System{T}(adjacency, dims)

    for eqc in eqconstraints
        eqc.parentid == origin.id && (eqc.parentid = nothing)
    end

    return system
end

function adjacencyMatrix(eqcs::Vector{<:EqualityConstraint}, bodies::Vector{<:Body}, ineqcs::Vector{<:InequalityConstraint})
    A = zeros(Bool, 0, 0)
    ids = Int64[]
    dims = zeros(Int64, 0)
    n = 0

    for eqc in eqcs
        eqcid = eqc.id
        A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
        n += 1
        append!(ids, eqcid)
        append!(dims,length(eqc))

        parentid = eqc.parentid
        if parentid ∉ ids
            A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
            n += 1
            append!(ids, parentid)
            append!(dims,6)
        end
        A[findfirst(x->x==eqcid, ids),findfirst(x->x==parentid, ids)] = true

        for bodyid in unique(eqc.childids)
            if bodyid ∉ ids
                A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
                n += 1
                append!(ids, bodyid)
                append!(dims,6)
            end
            A[findfirst(x->x==eqcid, ids),findfirst(x->x==bodyid, ids)] = true

            if eqc.isdamper
                A[findfirst(x->x==parentid, ids),findfirst(x->x==bodyid, ids)] = true
            end
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

    for ineqc in ineqcs 
        if typeof(ineqc.constraint) <: Impact
            ineqcid = ineqc.id
            A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
            n += 1
            append!(ids, ineqcid)
            append!(dims,length(ineqc)*2)

            parentid = ineqc.parentid
            A[findfirst(x->x==ineqcid, ids),findfirst(x->x==parentid, ids)] = true
        end
    end

    p = [findfirst(x->x==i,ids) for i=1:length(ids)][1:end-1]
    A = convert(Matrix{Int64}, A .| A')
    A = A[p,p]
    dims = dims[p]

    return A, dims
end

@inline getentry(system, id1, id2) = system.matrix_entries[id1, id2]
@inline getentry(system, id) = system.vector_entries[id]

function recursivedirectchildren!(system, id::Integer)
    dirs = copy(children(system, id))
    dirslocal = copy(dirs)
    for childid in dirslocal
        append!(dirs, recursivedirectchildren!(system, childid))
    end
    return dirs
end


# TODO does not include ineqcs yet
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

        for childid in system.cyclic_children[id]
            offdiagonal_L = getentry(system, id, childid)
            offdiagonal_U = getentry(system, childid, id)
            nc1 = first(rangeDict[childid])
            nc2 = last(rangeDict[childid])

            A[range,nc1:nc2] = offdiagonal_L.value
            A[nc1:nc2,range] = offdiagonal_U.value
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