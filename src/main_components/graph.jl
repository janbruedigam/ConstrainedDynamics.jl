struct Graph{N}
    directchildren::Vector{Vector{Int64}} # direct child nodes
    ineqchildren::Vector{Vector{Int64}} # direct child nodes for inequality constraints (contact)
    dampergrandchildren::Vector{Vector{Int64}} # direct grandchild nodes for bodies with damped children 
    loopchildren::Vector{Vector{Int64}} # successor nodes in a loop excluding direct children
    successors::Vector{Vector{Int64}} # direct and loop successors
    predecessors::Vector{Vector{Int64}} # direct parent and loop-opening predecessor(s?) (for numerics?)
    dampergrandparent::Vector{Vector{Int64}} # direct grandparent node for bodies with damped grandparents. Should only be one
    connections::Vector{Vector{Int64}} # direct connections
    springconnections::Vector{Vector{Int64}} # direct connections from bodies to eqconstraints with springs  
    damperconnections::Vector{Vector{Int64}} # direct connections for eqconstraints with damping  

    dfslist::SVector{N,Int64} # depth-first-seach list (dfslist[end] = root)
    rdfslist::SVector{N,Int64} # reverse dfslist

    dict::UnitDict{Base.OneTo{Int64},Int64} # maps ids to the graph-interal numbering (not dfs order). Necessary, e.g., if ids don't start at 1. 
    rdict::UnitDict{Base.OneTo{Int64},Int64} # reverse mapping
    activedict::UnitDict{Base.OneTo{Int64},Bool}

    function Graph(origin::Origin,bodies::Vector{<:Body},
            eqconstraints::Vector{<:EqualityConstraint},ineqconstraints::Vector{<:InequalityConstraint}
        )

        oid = origin.id
        adjacency, dict = adjacencyMatrix(eqconstraints, bodies)
        dfsgraph, dfslist, loops = dfs(adjacency, dict, oid)
        pat = pattern(dfsgraph, dict, loops)
        fil, originals = fillins(dfsgraph, pat, dict, loops)

        adjacency = deleteat(adjacency, dict[oid])
        dfsgraph = deleteat(dfsgraph, dict[oid])
        pat = deleteat(pat, dict[oid])
        fil = deleteat(fil, dict[oid])
        originals = deleteat(originals, dict[oid])
        dfslist = StaticArrays.deleteat(dfslist, length(dfslist))
        rdfslist = reverse(dfslist)

        for (id, ind) in dict
            ind > dict[oid] && (dict[id] = ind - 1)
        end
        pop!(dict, oid)
        rdict = Dict(ind => id for (id, ind) in dict)
        activedict = Dict(id => true for (id, ind) in dict)

        for eqc in eqconstraints
            eqc.parentid == oid && (eqc.parentid = nothing)
            activedict[eqc.id] = isactive(eqc)
        end

        for body in bodies
            activedict[body.id] = isactive(body)
        end

        for ineqc in ineqconstraints
            activedict[ineqc.id] = isactive(ineqc)
        end

        N = length(dict)

        adjacency = convert(Vector{SVector{N,Bool}}, adjacency)
        dfsgraph = convert(Vector{SVector{N,Bool}}, dfsgraph)
        pat = convert(Vector{SVector{N,Bool}}, pat)
        fil = convert(Vector{SVector{N,Bool}}, fil)
        originals = convert(Vector{SVector{N,Bool}}, originals)

        dirs = directchildren(dfslist, originals, dict)
        ineqs = ineqchildren(dfslist, bodies, ineqconstraints, dict)
        damps = dampergrandchildren(dfslist, eqconstraints, dict)
        loops = loopchildren(dfslist, fil, dict)
        succs = successors(dfslist, pat, dict)
        preds = predecessors(rdfslist, pat, dict)
        dampgrand = dampergrandparent(rdfslist, damps, dict, rdict)
        cons = connections(dfslist, adjacency, dict)
        springcons = springconnections(dfslist, bodies, eqconstraints, dict)
        dampcons = damperconnections(dfslist, bodies, eqconstraints, dict)

        dict = UnitDict(dict)
        rdict = UnitDict(rdict)
        activedict = UnitDict(activedict)

        new{N}(dirs, ineqs, damps, loops, succs, preds, dampgrand, cons, springcons, dampcons, dfslist, rdfslist, dict, rdict, activedict)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, graph::Graph{N}) where {N}
    summary(io, graph)
end


### Construction functions
function adjacencyMatrix(eqconstraints::Vector{<:EqualityConstraint}, bodies::Vector{<:Body})
    A = zeros(Bool, 0, 0)
    dict = Dict{Int64,Int64}()
    n = 0

    for constraint in eqconstraints
        childid = constraint.id
        A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
        dict[childid] = n += 1
        for bodyid in unique([constraint.parentid;constraint.childids])
            if !haskey(dict, bodyid)
                A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
                dict[bodyid] = n += 1
            end
            A[dict[childid],dict[bodyid]] = true
        end
    end
    for body in bodies # add unconnected bodies
        if !haskey(dict, body.id)
            A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
            dict[body.id] = n += 1
        end
    end

    A = A .| A'
    return convert(Matrix{Bool}, A), dict
end

function dfs(adjacency::Matrix, dict::Dict, originid::Integer)
    N = size(adjacency)[1]
    dfsgraph = zeros(Bool, N, N)
    dfslist = zeros(Int64, N)
    visited = zeros(Bool, N)
    loops = Vector{Int64}[]
    index = N

    dfslist[index] = originid
    visited[dict[originid]] = true
    dfs!(adjacency, dfsgraph, dict, dfslist, visited, loops, N, originid, -1)
    loops = loops[sortperm(sort.(loops))[1:2:length(loops)]] # removes double entries of loop connections and keeps the first found pair

    return dfsgraph, convert(SVector{N}, dfslist), loops
end

function dfs!(A::Matrix, Adfs::Matrix, dict::Dict, list::Vector, visited::Vector, loops::Vector{<:Vector{<:Integer}},
        index::Integer, currentid::Integer, parentid::Integer
    )
    
    i = dict[currentid]
    for (childid, j) in dict
        if A[i,j] && parentid != childid # connection from i to j in adjacency && not a direct connection back to the parent
            if visited[j]
                push!(loops, [childid,currentid]) # childid is actually a predecessor of currentid since it's a loop
            else
                index -= 1
                list[index] = childid
                visited[j] = true
                Adfs[i,j] = true
                index = dfs!(A, Adfs, dict, list, visited, loops, index, childid, currentid)
            end
        end
    end
    return index
end

function pattern(dfsgraph::Matrix, dict::Dict, loops::Vector{<:Vector{<:Integer}})
    pat = deepcopy(dfsgraph)

    for loop in loops # loop = [index of starting constraint, index of last body in loop]
        startid = loop[1]
        endid = loop[2]
        currentid = endid
        nextid = parent(dfsgraph, dict, currentid)
        while nextid != startid
            pat[dict[startid],dict[currentid]] = true
            currentid = nextid
            nextid = parent(dfsgraph, dict, nextid)
        end
    end

    return pat
end

function fillins(dfsgraph::Matrix, pattern::Matrix, dict::Dict, loops::Vector{<:Vector{<:Integer}})
    fil = deepcopy(dfsgraph .⊻ pattern) # xor so only fillins remain (+ to be removed loop closure since this is a child)
    originals = deepcopy(dfsgraph)

    for loop in loops # loop = [index of starting constraint, index of last body in loop]
        startid = loop[1]
        endid = loop[2]
        fil[dict[startid],dict[endid]] = false
        originals[dict[startid],dict[endid]] = true
    end

    return convert(Matrix{Bool}, fil), convert(Matrix{Bool}, originals)
end

function parent(dfsgraph::Matrix, dict::Dict, childid::Integer)
    j = dict[childid]
    for (parentid, i) in dict
        dfsgraph[i,j] && (return parentid)
    end
    return -1
end


### Graph functions

# this is done in order!
function directchildren(dfslist, dfsgraph, dict::Dict)
    N = length(dfslist)
    dirs = [Int64[] for i = 1:N]
    for i = 1:N
        for childid in dfslist
            dfsgraph[i][dict[childid]] && push!(dirs[i], childid)
        end
    end

    return dirs
end

function ineqchildren(dfslist, bodies, ineqconstraints, dict::Dict)
    N = length(dfslist)
    ineqs = [Int64[] for i = 1:N]
    for body in bodies
        for ineqc in ineqconstraints
            ineqc.parentid == body.id && push!(ineqs[dict[body.id]], ineqc.id)
        end
    end

    return ineqs
end

function dampergrandchildren(dfslist, eqconstraints, dict::Dict)
    N = length(dfslist)
    damps = [Int64[] for i = 1:N]
    for eqc in eqconstraints
        (!eqc.isdamper || eqc.parentid === nothing) && continue

        for id in unique(eqc.childids)
            push!(damps[dict[eqc.parentid]], id)
        end
    end

    return damps
end

# this is done in order!
function loopchildren(dfslist, fillins, dict::Dict)
    N = length(dfslist)
    loos = [Int64[] for i = 1:N]
    for i = 1:N
        for childid in dfslist
            fillins[i][dict[childid]] && push!(loos[i], childid)
        end
    end

    return loos
end

# this is done in order!
function successors(dfslist, pattern, dict::Dict)
    N = length(dfslist)
    sucs = [Int64[] for i = 1:N]
    for i = 1:N
        for childid in dfslist
            pattern[i][dict[childid]] && push!(sucs[i], childid)
        end
    end

    return sucs
end

# this is done in reverse order (but this is not really important for predecessors)
function predecessors(rdfslist, pattern, dict::Dict)
    N = length(rdfslist)
    preds = [Int64[] for i = 1:N]
    for i = 1:N
        for childid in rdfslist
            pattern[dict[childid]][i] && push!(preds[i], childid)
        end
    end

    return preds
end

# this is done in reverse order (but this is not really important for dampergrandparent, should only be one)
function dampergrandparent(rdfslist, dampergrandchildren, dict::Dict, rdict::Dict)
    N = length(rdfslist)
    dampgrand = [Int64[] for i = 1:N]
    for i = 1:N
        for grandparentid in rdfslist
            rdict[i] ∈ dampergrandchildren[dict[grandparentid]] && push!(dampgrand[i], grandparentid)
        end
    end

    return dampgrand
end

# this is done in order (but this is not really important for connections)
function connections(dfslist, adjacency, dict::Dict)
    N = length(dfslist)
    cons = [Int64[] for i = 1:N]
    for i = 1:N
        for childid in dfslist
            adjacency[i][dict[childid]] && push!(cons[i], childid)
        end
    end

    return cons
end

# this is done in order so the spring parent is always the first entry (but this is not really important for spring connections)
function springconnections(dfslist, bodies, eqconstraints, dict::Dict)
    N = length(dfslist)
    springs = [Int64[] for i = 1:N]
    for bodyid in dfslist
        for eqc in eqconstraints
            !eqc.isspring && continue

            eqc.parentid == bodyid && push!(springs[dict[bodyid]], eqc.id)

            for childid in unique(eqc.childids)
                childid == bodyid && push!(springs[dict[bodyid]], eqc.id)
            end
        end
    end

    return springs
end

# this is done in order so the damper parent is always the first entry 
function damperconnections(dfslist, bodies, eqconstraints, dict::Dict)
    N = length(dfslist)
    damps = [Int64[] for i = 1:N]
    for bodyid in dfslist
        for eqc in eqconstraints
            !eqc.isdamper && continue

            eqc.parentid == bodyid && push!(damps[dict[bodyid]], eqc.id)

            for childid in unique(eqc.childids)
                childid == bodyid && push!(damps[dict[bodyid]], eqc.id)
            end
        end
    end

    return damps
end

function recursivedirectchildren!(graph, id::Integer)
    dirs = copy(directchildren(graph, id))
    dirslocal = copy(dirs)
    for childid in dirslocal
        append!(dirs, recursivedirectchildren!(graph, childid))
    end
    return dirs
end


@inline directchildren(graph, id::Integer) = graph.directchildren[graph.dict[id]]
@inline dampergrandchildren(graph, id::Integer) = graph.dampergrandchildren[graph.dict[id]]
@inline ineqchildren(graph, id::Integer) = graph.ineqchildren[graph.dict[id]]
@inline loopchildren(graph, id::Integer) = graph.loopchildren[graph.dict[id]]
@inline successors(graph, id::Integer) = graph.successors[graph.dict[id]]
@inline predecessors(graph, id::Integer) = graph.predecessors[graph.dict[id]]
@inline dampergrandparent(graph, id::Integer) = graph.dampergrandparent[graph.dict[id]]
@inline connections(graph, id::Integer) = graph.connections[graph.dict[id]]
@inline springconnections(graph, id::Integer) = graph.springconnections[graph.dict[id]]
@inline damperconnections(graph, id::Integer) = graph.damperconnections[graph.dict[id]]
@inline isactive(graph, id::Integer) = graph.activedict[id]
@inline isinactive(graph, id::Integer) = !isactive(graph, id)

@inline function hasdirectchild(graph::Graph, id, childid)
    for val in graph.directchildren[graph.dict[id]]
        val == childid && (return true)
    end
    return false
end

@inline function hassuccessor(graph::Graph, id, childid)
    for val in graph.successors[graph.dict[id]]
        val == childid && (return true)
    end
    return false
end

@inline function haspredecessor(graph::Graph, id, parentid)
    for val in graph.predecessors[graph.dict[id]]
        val == parentid && (return true)
    end
    return false
end

@inline function activate!(graph::Graph, id::Integer)
    graph.activedict[id] = true
    return
end

@inline function deactivate!(graph::Graph, id::Integer)
    graph.activedict[id] = false
    return
end