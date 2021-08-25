"""
$(TYPEDEF)

An `EqualityConstraint` is a component of a [`Mechanism`](@ref) and is used to describe the kinematic relation between two or more [`Body`](@ref)s. 
Typically, an `EqualityConstraint` should not be created directly. Use the joint prototypes instead, for example: 
```julia
EqualityConstraint(Revolute(body1, body2, rotation_axis)).
```
# Important attributes
* `id`:       The unique ID of a constraint. Assigned when added to a `Mechanism`.  
* `name`:     The name of a constraint. The name is taken from a URDF or can be assigned by the user.
* `parentid`: The ID of the parent body.  
* `childids`: The IDs of the child bodies.  
"""
mutable struct EqualityConstraint{T,N,Nc,Cs} <: AbstractConstraint{T,N}
    id::Int64
    name::String
    isspring::Bool
    isdamper::Bool

    constraints::Cs
    parentid::Union{Int64,Nothing}
    childids::SVector{Nc,Int64}
    inds::SVector{Nc,SVector{2,Int64}} # indices for minimal coordinates, assumes joints

    λsol::Vector{SVector{N,T}}

    function EqualityConstraint(data...; name::String="")
        if typeof(data[1][1]) <: Friction
            data = data[1]
            global CURRENTID
            data[4].parentid = data[3].parentid = CURRENTID 
    
            T = Float64
    
            isspring = false
            isdamper = false
            parentid = data[2]
            childids = SA[data[3].id;data[4].id]
            constraints = (data[1],)
            inds = [[1;2] for i=1:2]
            N = 4
            Nc = 2
            
    
            λsol = [zeros(T, N) for i=1:2]
    
            return new{T,N,Nc,typeof(constraints)}(getGlobalID(), name, isspring, isdamper, constraints, parentid, childids, inds, λsol)

        else
            jointdata = Tuple{AbstractJoint,Int64,Int64}[]
            for info in data
                if info[1] isa AbstractJoint
                    push!(jointdata, info)
                else
                    for subinfo in info
                        push!(jointdata, subinfo)
                    end
                end
            end
    
            T = getT(jointdata[1][1])# .T
    
            isspring = false
            isdamper = false
            parentid = jointdata[1][2]
            childids = Int64[]
            constraints = AbstractJoint{T}[]
            inds = Vector{Int64}[]
            N = 0
            for set in jointdata
                set[1].spring != 0 && (isspring = true)
                set[1].damper != 0 && (isdamper = true)
                
                push!(constraints, set[1])
                @assert set[2] == parentid
                push!(childids, set[3])
    
                Nset = length(set[1])
                if isempty(inds)
                    push!(inds, [1;3-Nset])
                else
                    push!(inds, [last(inds)[2]+1;last(inds)[2]+3-Nset])
                end
                N += Nset
            end
            constraints = Tuple(constraints)
            Nc = length(constraints)
            
    
            λsol = [zeros(T, N) for i=1:2]
    
            return new{T,N,Nc,typeof(constraints)}(getGlobalID(), name, isspring, isdamper, constraints, parentid, childids, inds, λsol)
        end


        
    end
end


"""
    setPosition!(mechanism, eqconstraint, xθ)

Sets the minimal coordinates (vector) of joint `eqconstraint`. 

Revolute joint example:
    setPosition!(mechanism, geteqconstraint(mechanism, "joint_name"), [pi/2])
"""
function setPosition!(mechanism, eqc::EqualityConstraint, xθ; iter::Bool = true)
    if !iter
        _setPosition!(mechanism, eqc, xθ)
    else
        currentvals = minimalCoordinates(mechanism)
        _setPosition!(mechanism, eqc, xθ)
        for id in recursivedirectchildren!(mechanism.system, eqc.id)
            component = getcomponent(mechanism, id)
            if component isa EqualityConstraint
                _setPosition!(mechanism, component, currentvals[id])
            end
        end
    end

    return
end

# TODO make zero alloc
# TODO currently assumed constraints are in order and only joints which is the case unless very low level constraint setting
function _setPosition!(mechanism, eqc::EqualityConstraint{T,N,Nc}, xθ) where {T,N,Nc}
    @assert length(xθ)==3*Nc-N
    n = Int64(Nc/2)
    body1 = getbody(mechanism, eqc.parentid)
    for i = 1:n
        body2 = getbody(mechanism, eqc.childids[i])
        Δx = getPositionDelta(eqc.constraints[i], body1, body2, xθ[SUnitRange(eqc.inds[i][1],eqc.inds[i][2])]) 
        Δq = getPositionDelta(eqc.constraints[i+1], body1, body2, xθ[SUnitRange(eqc.inds[i+1][1],eqc.inds[i+1][2])])
        
        p1, p2 = eqc.constraints[i].vertices
        setPosition!(body1, body2; p1 = p1, p2 = p2, Δx = Δx, Δq = Δq)
    end
    return
end

# TODO make zero alloc
# TODO currently assumed constraints are in order and only joints which is the case unless very low level constraint setting
"""
    setVelocity!(mechanism, eqconstraint, vω)

Sets the minimal coordinate velocities (vector) of joint `eqconstraint`. Note that currently this function sets the velocity of the directly connected body.

Planar joint example:
    setVelocity!(mechanism, geteqconstraint(mechanism, jointid), [0.5;2.0])
"""
function setVelocity!(mechanism, eqc::EqualityConstraint{T,N,Nc}, vω) where {T,N,Nc}
    @assert length(vω)==3*Nc-N
    n = Int64(Nc/2)
    body1 = getbody(mechanism, eqc.parentid)
    for i = 1:n
        body2 = getbody(mechanism, eqc.childids[i])
        Δv = getVelocityDelta(eqc.constraints[i], body1, body2, vω[SUnitRange(eqc.inds[i][1],eqc.inds[i][2])])
        Δω = getVelocityDelta(eqc.constraints[i+1], body1, body2, vω[SUnitRange(eqc.inds[i+1][1],eqc.inds[i+1][2])])
        
        p1, p2 = eqc.constraints[i].vertices
        setVelocity!(body1, body2; p1 = p1, p2 = p2, Δv = Δv, Δω = Δω)
    end
    return
end

# TODO make zero alloc
"""
    setForce!(mechanism, eqconstraint, Fτ)

Sets the minimal coordinate forces (vector) of joint `eqconstraint`.

Prismatic joint example:
    setVelocity!(mechanism, geteqconstraint(mechanism, jointid), [-1.0])
"""
function setForce!(mechanism, eqc::EqualityConstraint{T,N,Nc}, Fτ::AbstractVector) where {T,N,Nc}
    @assert length(Fτ)==3*Nc-N
    for i = 1:Nc
        setForce!(eqc.constraints[i], Fτ[SUnitRange(eqc.inds[i][1], eqc.inds[i][2])])
    end
    return
end

function addForce!(mechanism, eqc::EqualityConstraint{T,N,Nc}, Fτ::AbstractVector) where {T,N,Nc}
    @assert length(Fτ)==3*Nc-N
    for i = 1:Nc
        addForce!(eqc.constraints[i], Fτ[SUnitRange(eqc.inds[i][1], eqc.inds[i][2])])
    end
    return
end

"""
    minimalCoordinates(mechanism, eqconstraint)

Gets the minimal coordinates of joint `eqconstraint`.
"""
@generated function minimalCoordinates(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(minimalCoordinates(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, eqc.childids[$i]))) for i = 1:Nc]
    return :(svcat($(vec...)))
end

"""
    minimalVelocities(mechanism, eqconstraint)

Gets the minimal coordinate velocities of joint `eqconstraint`.
"""
@generated function minimalVelocities(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(minimalVelocities(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, eqc.childids[$i]))) for i = 1:Nc]
    return :(svcat($(vec...)))
end

@inline function GtλTof!(mechanism, body::Body, eqc::EqualityConstraint)
    body.state.d -= zerodimstaticadjoint(∂g∂ʳpos(mechanism, eqc, body.id)) * eqc.λsol[2]
    return
end

@inline function springTof!(mechanism, body::Body, eqc::EqualityConstraint)
    body.state.d -= springforce(mechanism, eqc, body.id)
    return
end

@inline function damperTof!(mechanism, body::Body, eqc::EqualityConstraint)
    body.state.d -= damperforce(mechanism, eqc, body.id)
    return
end

@inline function damperToD!(mechanism, body::Body, eqc::EqualityConstraint)
    body.state.D -= diagonal∂damper∂ʳvel(mechanism, eqc, body.id)
    return
end

@generated function g(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(g(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, eqc.childids[$i]), mechanism.Δt)) for i = 1:Nc]
    return :(svcat($(vec...)))
end
@generated function gc(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(g(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, eqc.childids[$i]))) for i = 1:Nc]
    return :(svcat($(vec...)))
end

@inline function ∂damper∂ʳvel(mechanism, eqc::EqualityConstraint, ida::Integer, idb::Integer)
    ida == eqc.parentid ? (return ∂damper∂ʳvela(mechanism, eqc, ida, idb)) : (return ∂damper∂ʳvelb(mechanism, eqc, ida, idb))
end

@inline function springforce(mechanism, eqc::EqualityConstraint, id::Integer)
    id == eqc.parentid ? (return springforcea(mechanism, eqc, id)) : (return springforceb(mechanism, eqc, id))
end
@inline function damperforce(mechanism, eqc::EqualityConstraint, id::Integer)
    id == eqc.parentid ? (return damperforcea(mechanism, eqc, id)) : (return damperforceb(mechanism, eqc, id))
end

@generated function ∂g∂ʳposa(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [:(∂g∂ʳposa(eqc.constraints[$i], getbody(mechanism, id), getbody(mechanism, eqc.childids[$i]), eqc.childids[$i])) for i = 1:Nc]
    return :(vcat($(vec...)))
end
@generated function ∂g∂ʳposb(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [:(∂g∂ʳposb(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, id), eqc.childids[$i])) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂ʳvela(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [:(∂g∂ʳvela(eqc.constraints[$i], getbody(mechanism, id), getbody(mechanism, eqc.childids[$i]), eqc.childids[$i], mechanism.Δt)) for i = 1:Nc]
    return :(vcat($(vec...)))
end
@generated function ∂g∂ʳvelb(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [:(∂g∂ʳvelb(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, id), eqc.childids[$i], mechanism.Δt)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

# Currently assumes no coupling between translational and rotational velocities
@inline function diagonal∂damper∂ʳvel(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    D = szeros(T, 6, 6)
    for i=1:Nc
        if id == eqc.parentid || id == eqc.childids[i]
            D += diagonal∂damper∂ʳvel(eqc.constraints[i])
        end
    end
    return D
end
@inline function offdiagonal∂damper∂ʳvel(mechanism, eqc::EqualityConstraint{T,N,Nc}, id1::Integer, id2::Integer) where {T,N,Nc}
    D = szeros(T, 6, 6)
    body1 = getbody(mechanism, id1)
    body2 = getbody(mechanism, id2)
    for i=1:Nc
        D += offdiagonal∂damper∂ʳvel(eqc.constraints[i], body1, body2, eqc.childids[i])
    end
    return D
end

@inline function springforcea(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += springforcea(eqc.constraints[i], getbody(mechanism, id), getbody(mechanism, eqc.childids[i]), eqc.childids[i])
    end
    return vec
end
@inline function springforceb(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += springforceb(eqc.constraints[i], getbody(mechanism, eqc.parentid), getbody(mechanism, id), eqc.childids[i])
    end
    return vec
end

@inline function damperforcea(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += damperforcea(eqc.constraints[i], getbody(mechanism, id), getbody(mechanism, eqc.childids[i]), eqc.childids[i])
    end
    return vec
end
@inline function damperforceb(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += damperforceb(eqc.constraints[i], getbody(mechanism, eqc.parentid), getbody(mechanism, id), eqc.childids[i])
    end
    return vec
end

@generated function ∂Fτ∂ua(mechanism, eqc::EqualityConstraint{T,N,Nc}, id) where {T,N,Nc}
    vec = [:(∂Fτ∂ua(eqc.constraints[$i], getbody(mechanism, id), getbody(mechanism, eqc.childids[$i]), eqc.childids[$i])) for i = 1:Nc]
    return :(hcat($(vec...)))
end
@generated function ∂Fτ∂ub(mechanism, eqc::EqualityConstraint{T,N,Nc}, id) where {T,N,Nc}
    vec = [:(∂Fτ∂ub(eqc.constraints[$i], getbody(mechanism, eqc.parentid), getbody(mechanism, id), eqc.childids[$i])) for i = 1:Nc]
    return :(hcat($(vec...)))
end

@inline function applyFτ!(eqc::EqualityConstraint{T,N,Nc}, mechanism, clear::Bool=true) where {T,N,Nc}
    for i=1:Nc
        applyFτ!(eqc.constraints[i], getbody(mechanism, eqc.parentid), getbody(mechanism, eqc.childids[i]), clear)
    end
    return
end

# Derivatives NOT accounting for quaternion specialness

function ∂g∂posc(mechanism, eqc::EqualityConstraint, id::Integer)
    id == eqc.parentid ? (return ∂g∂posac(mechanism, eqc, id)) : (return ∂g∂posbc(mechanism, eqc, id))
end

function ∂g∂posac(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [hcat(∂g∂posac(eqc.constraints[i], getbody(mechanism, id), getbody(mechanism, eqc.childids[i]), eqc.childids[i])) for i = 1:Nc]
    return vcat(vec...)
end
function ∂g∂posbc(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [hcat(∂g∂posbc(eqc.constraints[i], getbody(mechanism, eqc.parentid), getbody(mechanism, id), eqc.childids[i])) for i = 1:Nc]
    return vcat(vec...)
end

function Base.cat(eqc1::EqualityConstraint{T,N1,Nc1}, eqc2::EqualityConstraint{T,N2,Nc2}) where {T,N1,N2,Nc1,Nc2}
    @assert eqc1.parentid == eqc2.parentid "Can only concatenate constraints with the same parentid"
    parentid = eqc1.parentid
    if parentid === nothing 
        parentid = -1
        nothingflag = true
    else
        nothingflag = false
    end

    constraints = [[eqc1.constraints[i] for i=1:Nc1];[eqc2.constraints[i] for i=1:Nc2]]
    childids = [[eqc1.childids[i] for i=1:Nc1];[eqc2.childids[i] for i=1:Nc2]]

    eqc = EqualityConstraint([(constraints[i],parentid,childids[i]) for i=1:Nc1+Nc2]..., name="combined_"*eqc1.name*"_and_"*eqc2.name)
    nothingflag && (eqc.parentid = nothing)

    return eqc
end