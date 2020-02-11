mutable struct InequalityConstraint{T,N,Cs} <: AbstractConstraint{T}
    id::Int64

    constraints::Cs
    pid::Int64
    # bodyid::Int64

    sl0::Float64
    sl1::Float64
    ga0::Float64
    ga1::Float64

    function InequalityConstraint(body::Body{T}) where T
        N = 1
        pid = body.id
        constraints = Impact(body)

        s0 = 0
        s1 = 0
        γ0 = 0
        γ1 = 0

        new{T,N,typeof(constraints)}(getGlobalID(),constraints,pid,s0,s1,γ0,γ1)
    end
end


Base.length(c::InequalityConstraint{T,N}) where {T,N} = N

function g(ineq::InequalityConstraint,body::Body,mechanism)
    g(ineq,ineq.constraints,body,mechanism.dt,mechanism.No,mechanism.μ)
end

function diagval(ineq::InequalityConstraint,body::Body,dt)
    diagval(ineq,ineq.constraints,body,dt)
end
