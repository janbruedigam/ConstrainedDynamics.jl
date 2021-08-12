struct Structure{N}
    system::System{N}
    dict::Dict{Int64,Int64} # maps ids to the graph-interal numbering (not dfs order). Necessary, e.g., if ids don't start at 1. 
    rdict::Dict{Int64,Int64} # reverse mapping

    function Structure(origin::Origin{T}, bodies::Vector{<:Body},
            eqconstraints::Vector{<:EqualityConstraint}, ineqconstraints::Vector{<:InequalityConstraint}
        ) where T

        adjacency, dict, dims = adjacencyMatrix2(eqconstraints, bodies)
        oid = origin.id
        graphoid = dict[oid]
        adjacency = deleteat(adjacency, graphoid)
        dims = deleteat!(dims, graphoid)
        for (id, ind) in dict
            ind > graphoid && (dict[id] = ind - 1)
        end
        pop!(dict, oid)
        rdict = Dict(ind => id for (id, ind) in dict)

        system = System{T}(adjacency, dims)

        new{length(dims)}(system, dict, rdict)
    end
end

function adjacencyMatrix2(eqconstraints::Vector{<:EqualityConstraint}, bodies::Vector{<:Body})
    A = zeros(Bool, 0, 0)
    dict = Dict{Int64,Int64}()
    dims = zeros(Int64, 0)
    n = 0

    for constraint in eqconstraints
        childid = constraint.id
        A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
        dict[childid] = n += 1
        append!(dims,length(constraint))
        for bodyid in unique([constraint.parentid;constraint.childids])
            if !haskey(dict, bodyid)
                A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
                dict[bodyid] = n += 1
                append!(dims,6)
            end
            A[dict[childid],dict[bodyid]] = true
        end
    end
    for body in bodies # add unconnected bodies
        if !haskey(dict, body.id)
            A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
            dict[body.id] = n += 1
            append!(dims,6)
        end
    end

    A = A .| A'
    return convert(Matrix{Int64}, A), dict, dims
end

@inline function getentry(structure, id1, id2)
    dict = structure.dict
    structure.system.matrix_entries[(dict[id1],dict[id2])]
end
@inline function getentry(structure, id)
    structure.system.vector_entries[structure.dict[id]]
end