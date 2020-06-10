mutable struct State{T}
    order::Integer

    # Continuous states
    xc::SVector{3,T}
    qc::UnitQuaternion{T}
    vc::SVector{3,T}
    ωc::SVector{3,T}

    # Knot points
    xk::Vector{SVector{3,T}}
    qk::Vector{UnitQuaternion{T}}
    Fk::Vector{SVector{3,T}}
    τk::Vector{SVector{3,T}}

    # Current solution estimate [before step;after step] (xsol and qsol are not set in code since they are trivially x2 and q2)
    xsol::Vector{SVector{3,T}}
    qsol::Vector{UnitQuaternion{T}}
    vsol::Vector{SVector{3,T}}
    ωsol::Vector{SVector{3,T}}

    function State{T}() where T
        xc = zeros(T, 3)
        qc = one(UnitQuaternion{T})
        vc = zeros(T, 3)
        ωc = zeros(T, 3)

        xk = [zeros(T, 3)]
        qk = [one(UnitQuaternion{T})]
        Fk = [zeros(T, 3)]
        τk = [zeros(T, 3)]

        xsol = [zeros(T, 3) for i=1:2]
        qsol = [one(UnitQuaternion{T}) for i=1:2]
        vsol = [zeros(T, 3) for i=1:2]
        ωsol = [zeros(T, 3) for i=1:2]
        new{T}(0, xc, qc, vc, ωc, xk, qk, Fk, τk, xsol, qsol, vsol, ωsol)
    end
end

function initknotpoints!(state::State, order)
    state.order = order

    state.xk = [state.xk[1] for i = 1:order]
    state.qk = [state.qk[1] for i = 1:order]
    state.Fk = [state.Fk[1] for i = 1:order]
    state.τk = [state.τk[1] for i = 1:order]

    return
end