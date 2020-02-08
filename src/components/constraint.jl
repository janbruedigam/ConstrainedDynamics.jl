mutable struct Constraint{T,N,Nc,Cs} <: Component{T}
    id::Int64

    constraints::Cs
    pid::Union{Int64,Nothing}
    bodyids::SVector{Nc,Int64}

    s0::SVector{N,T}
    s1::SVector{N,T}

    function Constraint(data...)
        jointdata = Vector{Tuple{Joint,Int64,Int64}}(undef,0)
        for info in data
            if typeof(info[1]) <: Joint
                push!(jointdata,info)
            else
                for subinfo in info
                    push!(jointdata,subinfo)
                end
            end
        end

        T = getT(jointdata[1][1])#.T

        pid = jointdata[1][2]
        bodyids = Vector{Int64}(undef,0)
        constraints = Vector{Joint{T}}(undef,0)
        N = 0
        for set in jointdata
            push!(constraints,set[1])
            @assert set[2] == pid
            push!(bodyids,set[3])
            N += getNc(set[1])#.Nc
        end
        constraints = Tuple(constraints)
        Nc = length(constraints)

        s0 = @SVector zeros(T,N)
        s1 = @SVector zeros(T,N)

        new{T,N,Nc,typeof(constraints)}(getGlobalID(),constraints,pid,bodyids,s0,s1)
    end

    # Constraint(joint) = Constraint(joint...)
end

Base.length(c::Constraint{T,N}) where {T,N} = N

@generated function g(c::Constraint{T,N,Nc},robot) where {T,N,Nc}
    vec = [:(g(c.constraints[$i],getbody(robot,c.pid),getbody(robot,c.bodyids[$i]),robot.dt,robot.No)) for i=1:Nc]
    :(vcat($(vec...)))
end

@inline function ∂g∂pos(c::Constraint,id::Int64,robot)
    id == c.pid ? ∂g∂posa(c,id,robot) : ∂g∂posb(c,id,robot)
end

@inline function ∂g∂vel(c::Constraint,id::Int64,robot)
    id == c.pid ? ∂g∂vela(c,id,robot) : ∂g∂velb(c,id,robot)
end

@generated function ∂g∂posa(c::Constraint{T,N,Nc},id::Int64,robot) where {T,N,Nc}
    vec = [:(∂g∂posa(c.constraints[$i],getbody(robot,id),getbody(robot,c.bodyids[$i]),robot.No)) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂posb(c::Constraint{T,N,Nc},id::Int64,robot) where {T,N,Nc}
    vec = [:(∂g∂posb(c.constraints[$i],getbody(robot,c.pid),getbody(robot,id),robot.No)) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂vela(c::Constraint{T,N,Nc},id::Int64,robot) where {T,N,Nc}
    vec = [:(∂g∂vela(c.constraints[$i],getbody(robot,id),getbody(robot,c.bodyids[$i]),robot.dt,robot.No)) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂velb(c::Constraint{T,N,Nc},id::Int64,robot) where {T,N,Nc}
    vec = [:(∂g∂velb(c.constraints[$i],getbody(robot,c.pid),getbody(robot,id),robot.dt,robot.No)) for i=1:Nc]
    return :(vcat($(vec...)))
end
