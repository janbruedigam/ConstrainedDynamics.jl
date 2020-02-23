mutable struct InequalityConstraint{T,N,Cs} <: AbstractConstraint{T}
    id::Int64

    constraints::Cs
    pid::Int64
    # bodyid::Int64

    s0::SVector{N,T}
    s1::SVector{N,T}
    γ0::SVector{N,T}
    γ1::SVector{N,T}

    function InequalityConstraint(input)
        impact,body = input

        T = Float64
        N = 1
        pid = body.id
        constraints = impact

        s0 = ones(T,N)
        s1 = ones(T,N)
        γ0 = ones(T,N)
        γ1 = ones(T,N)

        new{T,N,typeof(constraints)}(getGlobalID(),constraints,pid,s0,s1,γ0,γ1)
    end
end


Base.length(::InequalityConstraint{T,N}) where {T,N} = N

function g(ineqc::InequalityConstraint,mechanism)
    val = g(ineqc.constraints,getbody(mechanism,ineqc.pid),mechanism.dt,mechanism.No)
    val - ineqc.s1[1]
end

function hμ(ineqc::InequalityConstraint,mechanism)
    ineqc.s1[1]*ineqc.γ1[1] - mechanism.μ
end

function h(ineqc::InequalityConstraint,mechanism)
    ineqc.s1[1]*ineqc.γ1[1]
end

function dynineq(ineqc::InequalityConstraint,body::Body,mechanism)
    dynineq(ineqc,ineqc.constraints,body,mechanism.dt,mechanism.No,mechanism.μ)
end

function diagval(ineqc::InequalityConstraint,body::Body,dt)
    diagval(ineqc,ineqc.constraints,body,dt)
end

function ∂g∂pos(ineqc::InequalityConstraint,id,mechanism)
    ineqc.constraints.Nx
end
