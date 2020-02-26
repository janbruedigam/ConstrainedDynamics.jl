mutable struct InequalityConstraint{T,N,Cs} <: AbstractConstraint{T}
    id::Int64

    constraints::Cs
    pid::Int64
    # bodyid::Int64

    s0::SVector{N,T}
    s1::SVector{N,T}
    γ0::SVector{N,T}
    γ1::SVector{N,T}

    # function InequalityConstraint(input1,input2)
    #     impact1,body = input1
    #     impact2,body = input2
    #
    #     T = Float64
    #     N = 2
    #     pid = body.id
    #     constraints = (impact1,impact2)
    #
    #     s0 = ones(T,N)
    #     s1 = ones(T,N)
    #     γ0 = ones(T,N)
    #     γ1 = ones(T,N)
    #
    #     new{T,N,typeof(constraints)}(getGlobalID(),constraints,pid,s0,s1,γ0,γ1)
    # end

    function InequalityConstraint(data...)
        contactdata = Tuple{Contact,Int64}[]
        for info in data
            if typeof(info[1]) <: Contact
                push!(contactdata,info)
            else
                for subinfo in info
                    push!(contactdata,subinfo)
                end
            end
        end

        T = getT(contactdata[1][1])#.T

        pid = contactdata[1][2]
        bodyids = Int64[]
        constraints = Contact{T}[]
        N = 0
        for set in contactdata
            push!(constraints,set[1])
            @assert set[2] == pid
            N += 1 # getNc(set[1])
        end
        constraints = Tuple(constraints)
        # Nc = length(constraints)

        s0 = ones(T,N)
        s1 = ones(T,N)
        γ0 = ones(T,N)
        γ1 = ones(T,N)

        new{T,N,typeof(constraints)}(getGlobalID(),constraints,pid,s0,s1,γ0,γ1)
    end
end


Base.length(::InequalityConstraint{T,N}) where {T,N} = N

function setOne(ineqc::InequalityConstraint{T,N}) where {T,N}
    ineqc.s0 = @SVector ones(T,N)
    ineqc.s1 = @SVector ones(T,N)
    ineqc.γ0 = @SVector ones(T,N)
    ineqc.γ1 = @SVector ones(T,N)
    return 
end

@generated function gs(ineqc::InequalityConstraint{T,N},mechanism) where {T,N}
    vec = [:(g(ineqc.constraints[$i],getbody(mechanism,ineqc.pid),mechanism.dt,mechanism.No) - ineqc.s1[$i]) for i=1:N]
    # :(vcat($(vec...)))
    :(SVector{N,T}($(vec...)))
end

function gs(ineqc::InequalityConstraint{T,1},mechanism) where {T}
    g(ineqc.constraints[1],getbody(mechanism,ineqc.pid),mechanism.dt,mechanism.No) - ineqc.s1[1]
end

@generated function g(ineqc::InequalityConstraint{T,N},mechanism) where {T,N}
    vec = [:(g(ineqc.constraints[$i],getbody(mechanism,ineqc.pid),mechanism.dt,mechanism.No)) for i=1:N]
    # :(vcat($(vec...)))
    :(SVector{N,T}($(vec...)))
end

function g(ineqc::InequalityConstraint{T,1},mechanism) where {T}
    g(ineqc.constraints[1],getbody(mechanism,ineqc.pid),mechanism.dt,mechanism.No)
end

function hμ(ineqc::InequalityConstraint{T},μ::T) where T
    ineqc.s1.*ineqc.γ1 .- μ
end

h(ineqc::InequalityConstraint) = hμ(ineqc,0.)

function schurf(ineqc::InequalityConstraint{T,N},body,mechanism) where {T,N}
    val = @SVector zeros(T,6)
    for i=1:N
        val += schurf(ineqc,ineqc.constraints[i],i,body,mechanism)
    end
    return val
end

function schurD(ineqc::InequalityConstraint{T,N},body,mechanism) where {T,N}
    val = @SMatrix zeros(T,6,6)
    for i=1:N
        val += schurD(ineqc,ineqc.constraints[i],i,body,mechanism)
    end
    return val
end

@generated function ∂g∂pos(ineqc::InequalityConstraint{T,N},body,mechanism) where {T,N}
    vec = [:(∂g∂pos(ineqc.constraints[$i])) for i=1:N]
    :(vcat($(vec...)))
end

@generated function ∂g∂vel(ineqc::InequalityConstraint{T,N},body,mechanism) where {T,N}
    vec = [:(∂g∂vel(ineqc.constraints[$i],mechanism.dt)) for i=1:N]
    :(vcat($(vec...)))
end

function extrafriction!(ineqc::InequalityConstraint{T,N},body,mechanism) where {T,N}
    for i=1:N
        extrafriction!(ineqc,ineqc.constraints[i],i,body,mechanism)
    end
    return
end
