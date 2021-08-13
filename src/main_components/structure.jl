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