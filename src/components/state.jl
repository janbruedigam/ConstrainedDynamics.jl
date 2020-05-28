mutable struct State{T}
    # Continuous states
    xc::SVector{3,T}
    qc::Quaternion{T}
    vc::SVector{3,T}
    ωc::SVector{3,T}

    # Knot points
    xd::Vector{SVector{3,T}}
    qd::Vector{Quaternion{T}}

    # Current solution estimate [before step;after step] (xsol and qsol are not set in code since they are trivially x2 and q2)
    xsol::Vector{SVector{3,T}}
    qsol::Vector{Quaternion{T}}
    vsol::Vector{SVector{3,T}}
    ωsol::Vector{SVector{3,T}}

    function State{T}() where T
        xc = zeros(T, 3)
        qc = Quaternion{T}()
        vc = zeros(T, 3)
        ωc = zeros(T, 3)

        xd = [zeros(T, 3)]
        qd = [Quaternion{T}()]

        xsol = [zeros(T, 3) for i=1:2]
        qsol = [Quaternion{T}() for i=1:2]
        vsol = [zeros(T, 3) for i=1:2]
        ωsol = [zeros(T, 3) for i=1:2]
        new{T}(xc, qc, vc, ωc, xd, qd, xsol, qsol, vsol, ωsol)
    end
end


