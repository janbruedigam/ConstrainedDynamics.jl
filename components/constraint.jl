mutable struct Constraint{T,N,Nc,Cs} <: Component{T}
    id::Int64

    constraints::Cs
    pid::Union{Int64,Nothing}
    linkids::SVector{Nc,Int64}

    s0::SVector{N,T}
    s1::SVector{N,T}

    function Constraint(jointdata...)
        T = getT(jointdata[1][1])#.T

        pid = jointdata[1][2]
        linkids = Vector{Int64}(undef,0)
        constraints = Vector{Joint{T}}(undef,0)
        N = 0
        for set in jointdata
            push!(constraints,set[1])
            @assert set[2] == pid
            push!(linkids,set[3])
            N += getNc(set[1])#.Nc
        end
        constraints = Tuple(constraints)
        Nc = length(constraints)

        s0 = @SVector zeros(T,N)
        s1 = @SVector zeros(T,N)

        new{T,N,Nc,typeof(constraints)}(getGlobalID(),constraints,pid,linkids,s0,s1)
    end
end

Base.length(c::Constraint{T,N}) where {T,N} = N

@generated function g(c::Constraint{T,N,Nc},robot) where {T,N,Nc}
    vec = [:(g(c.constraints[$i],getlink(robot,c.pid),getlink(robot,c.linkids[$i]),robot.dt,robot.No)) for i=1:Nc]
    :(vcat($(vec...)))
end

@inline function ∂g∂pos(c::Constraint,id::Int64,robot)
    id == c.pid ? ∂g∂posa(c,id,robot) : ∂g∂posb(c,id,robot)
end

@inline function ∂g∂vel(c::Constraint,id::Int64,robot)
    id == c.pid ? ∂g∂vela(c,id,robot) : ∂g∂velb(c,id,robot)
end

@generated function ∂g∂posa(c::Constraint{T,N,Nc},id::Int64,robot) where {T,N,Nc}
    vec = [:(∂g∂posa(c.constraints[$i],getlink(robot,id),getlink(robot,c.linkids[$i]),robot.No)) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂posb(c::Constraint{T,N,Nc},id::Int64,robot) where {T,N,Nc}
    vec = [:(∂g∂posb(c.constraints[$i],getlink(robot,c.pid),getlink(robot,id),robot.No)) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂vela(c::Constraint{T,N,Nc},id::Int64,robot) where {T,N,Nc}
    vec = [:(∂g∂vela(c.constraints[$i],getlink(robot,id),getlink(robot,c.linkids[$i]),robot.dt,robot.No)) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂velb(c::Constraint{T,N,Nc},id::Int64,robot) where {T,N,Nc}
    vec = [:(∂g∂velb(c.constraints[$i],getlink(robot,c.pid),getlink(robot,id),robot.dt,robot.No)) for i=1:Nc]
    return :(vcat($(vec...)))
end
