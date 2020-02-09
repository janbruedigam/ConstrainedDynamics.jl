mutable struct InequalityConstraint{T,Cs} <: AbstractConstraint{T}
    id::Int64

    constraints::Cs
    pid::Int64
    # bodyid::Int64

    sl0::Float64
    sl1::Float64
    ga0::Float64
    ga1::Float64

    function InequalityConstraint(body::Body{T}) where T
        pid = body.id
        constraints = Impact(body)

        s0 = 0
        s1 = 0
        γ0 = 0
        γ1 = 0

        new{T,typeof(constraints)}(getGlobalID(),constraints,pid,s0,s1,γ0,γ1)
    end
end


function g(ineq::InequalityConstraint,body::Body,dt,No,μ)
    g(ineq.constraints,body,dt,No,μ)
end

function diagval(ineq::InequalityConstraint,body::Body,dt)
    diagval(ineq.constraints,body,dt)
end
