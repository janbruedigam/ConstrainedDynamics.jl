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

    # Current solution estimate [before step;after step]
    xsol::Vector{SVector{3,T}}
    qsol::Vector{UnitQuaternion{T}}
    vsol::Vector{SVector{3,T}}
    ωsol::Vector{SVector{3,T}}
    
    # Current equations of motion evaluation
    d::SVector{6,T}

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

        d = zeros(T, 6)

        new{T}(0, xc, qc, vc, ωc, xk, qk, Fk, τk, xsol, qsol, vsol, ωsol, d)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, state::State{T}) where {T}
    summary(io, state)
    println(io,"")
    println(io,"xc:   "*string(state.xc))
    println(io,"qc:   "*string(state.qc))
    println(io,"ωc:   "*string(state.ωc))
    println(io,"xk:   "*string(state.xk))
    println(io,"qk:   "*string(state.qk))
    println(io,"Fk:   "*string(state.Fk))
    println(io,"τk:   "*string(state.τk))
    println(io,"xsol: "*string(state.xsol))
    println(io,"qsol: "*string(state.qsol))
    println(io,"vsol: "*string(state.vsol))
    println(io,"ωsol: "*string(state.ωsol))
    println(io,"d:    "*string(state.d))
end


function initknotpoints!(state::State, order)
    state.order = order

    state.xk = [state.xk[1] for i = 1:order]
    state.qk = [state.qk[1] for i = 1:order]
    state.Fk = [state.Fk[1] for i = 1:order]
    state.τk = [state.τk[1] for i = 1:order]

    return
end