mutable struct State{T}
    # Continuous
    xc::Vector{SVector{3,T}}
    qc::Vector{Quaternion{T}}
    vc::Vector{SVector{3,T}}
    ωc::Vector{SVector{3,T}}

    # Discrete
    xd::Vector{SVector{3,T}}
    qd::Vector{Quaternion{T}}

    function State{T}() where T
        xc = [zeros(T, 3)]
        qc = [Quaternion{T}()]
        vc = [zeros(T, 3)]
        ωc = [zeros(T, 3)]

        xd = [zeros(T, 3)]
        qd = [Quaternion{T}()]
        new{T}(xc, qc, vc, ωc, xd, qd)
    end
end

function discretizex(xc, vc)
    xd1 = xc
    xd2 = xc
    return xd1, xd2
end

function discretizeq(qc, ωc)
    qd1 = qc
    qd2 = qc
    return qd1, qd2
end
