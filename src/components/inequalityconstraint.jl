mutable struct InequalityConstraint{T,N,Cs} <: AbstractConstraint{T}
    id::Int64

    constraints::Cs
    pid::Int64
    # bodyid::Int64

    sl0::Float64
    sl1::Float64
    ga0::Float64
    ga1::Float64

    slf0::Float64
    slf1::Float64
    psi0::Float64
    psi1::Float64
    b0::Vector{Float64}
    b1::Vector{Float64}

    function InequalityConstraint(body::Body{T}) where T
        N = 1
        pid = body.id
        constraints = Impact(body)

        s0 = 1.
        s1 = 1.
        γ0 = 1.
        γ1 = 1.

        new{T,N,typeof(constraints)}(getGlobalID(),constraints,pid,s0,s1,γ0,γ1,1,1,1,1,zeros(2),zeros(2))
    end
end


Base.length(c::InequalityConstraint{T,N}) where {T,N} = N

function g(c::InequalityConstraint,mechanism)
    g(c,c.constraints,getbody(mechanism,c.pid),mechanism.dt,mechanism.No)
end

function hμ(c::InequalityConstraint,mechanism)
    [c.sl1*c.ga1 - mechanism.μ;c.slf1*c.psi1 - mechanism.μ]
end

function h(c::InequalityConstraint,mechanism)
    [c.sl1*c.ga1;c.slf1*c.psi1]
end

function dynineq(ineq::InequalityConstraint,body::Body,mechanism)
    dynineq(ineq,ineq.constraints,body,mechanism.dt,mechanism.No,mechanism.μ)
end

function diagval(ineq::InequalityConstraint,body::Body,dt)
    diagval(ineq,ineq.constraints,body,dt)
end
