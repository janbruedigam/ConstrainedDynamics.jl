mutable struct Body{T} <: AbstractBody{T}
    id::Int64
    name::String
    active::Bool

    m::T
    J::SMatrix{3,3,T,9}

    state::State{T}


    function Body(m::T, J::AbstractArray{T,2}; name::String="") where T
        new{T}(getGlobalID(), name, true, m, J, State{T}())
    end

    function Body(shape::Shape; name::String="")
        body = Body(shape.m, shape.J; name=name)
        push!(shape.bodyids, body.id)
        return body
    end

    function Body(::Type{T},contents...) where T
        new{T}(getGlobalID(), contents...)
    end
end

mutable struct Origin{T} <: AbstractBody{T}
    id::Int64
    name::String

    Origin{T}(; name::String="") where T = new{T}(getGlobalID(), name)

    function Origin(shape::Shape{T}; name::String="") where T
        origin = Origin{T}(name=name)
        push!(shape.bodyids, origin.id)
        return origin
    end
end

function Base.deepcopy(b::Body{T}) where T
    contents = []
    for i = 2:getfieldnumber(b)
        push!(contents, deepcopy(getfield(b, i)))
    end

    return Body(T, contents...)
end

function Base.deepcopy(b::Body{T}, shape::Shape) where T
    contents = []
    
    for i = 2:getfieldnumber(b)
        push!(contents, deepcopy(getfield(b, i)))
    end

    body = Body(T, contents...)
    @assert (body.m == shape.m && body.J == shape.J)
    push!(shape.bodyids, body.id)
    return body
end

Base.length(::Body) = 6

@inline getM(body::Body{T}) where T = [[I*body.m;szeros(T,3,3)] [szeros(T,3,3);body.J]]


@inline function ∂dyn∂vel(body::Body{T}, Δt) where T
    J = body.J
    ω2 = body.state.ωsol[2]
    sq = sqrt(4 / Δt^2 - ω2' * ω2)

    dynT = I * body.m / Δt
    dynR = skewplusdiag(ω2, sq) * J - J * ω2 * (ω2' / sq) - skew(J * ω2)

    Z = szeros(T, 3, 3)

    return [[dynT; Z] [Z; dynR]]
end

@inline function ∂zp1∂z(mechanism, body::Body{T}, Δt) where T
    state = body.state
    Z = szeros(T,3,3)
    E = SMatrix{3,3,T,9}(I)
    Z2 = [Z Z]
    Z2t = [Z;Z]

    # Position
    AposT = [I E*Δt] # TODO is there a UniformScaling way for this instead of E?
    BposT = Z

        # This calculates the ϵ for q⊗Δq = q⊗(1 ϵᵀ)ᵀ
    AposR = VLᵀmat(state.qsol[2]) * [Rmat(ωbar(state.ωc, Δt))*LVᵀmat(state.qc) Lmat(state.qc)*derivωbar(state.ωc, Δt)]
    BposR = Z

    # Velocity
    AvelT = [Z I]
    BvelT = I*Δt/body.m

    J = body.J
    ω1 = state.ωc
    ω2 = state.ωsol[2]
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)

    ω1func = skewplusdiag(-ω1, sq1) * J - J * ω1 * (ω1' / sq1) + skew(J * ω1)
    ω2func = skewplusdiag(ω2, sq2) * J - J * ω2 * (ω2' / sq2) - skew(J * ω2)
    invω2func = inv(ω2func)
    AvelR = [Z invω2func*ω1func]
    BvelR = 2*invω2func


    AT = [[AposT;AvelT] [Z2;Z2]]
    AR = [[Z2;Z2] [AposR;AvelR]]
    BT = [[BposT;BvelT] Z2t]
    BR = [Z2t [BposR;BvelR]]

    return [AT;AR], [BT;BR]
end

@inline function ∂F∂z(mechanism, body::Body{T}, Δt) where T
    state = body.state
    Z3 = szeros(T,3,3)
    Z43 = szeros(T,4,3)
    Z6 = szeros(T,6,6)
    Z63 = szeros(T,6,3)
    Z73 = szeros(T,7,3)
    Z76 = szeros(T,7,6)
    E3 = SMatrix{3,3,T,9}(I)
    

    # Position
    AposT = [-I -E3*Δt] # TODO is there a UniformScaling way for this instead of E3?
    BposT = Z3

    AposR = [-Rmat(ωbar(state.ωc, Δt))*LVᵀmat(state.qc) -Lmat(state.qc)*derivωbar(state.ωc, Δt)]
    BposR = Z43

    # Velocity
    AvelT = [Z3 -I*body.m/Δt]
    BvelT = -I

    J = body.J
    ω1 = state.ωc
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)

    ω1func = skewplusdiag(-ω1, sq1) * J - J * ω1 * (ω1' / sq1) + skew(J * ω1)
    AvelR = [Z3 ω1func]
    BvelR = -2.0*I


    AT = [[AposT;AvelT] Z6]
    AR = [Z76 [AposR;AvelR]]
    BT = [[BposT;BvelT] Z63]
    BR = [Z73 [BposR;BvelR]]

    return [AT;AR], [BT;BR]
end

@inline function ∂F∂fz(mechanism, body::Body{T}, Δt) where T
    state = body.state
    Z3 = szeros(T,3,3)
    Z4 = szeros(T,4,4)
    Z34 = szeros(T,3,4)
    Z43 = szeros(T,4,3)
    Z67 = szeros(T,6,7)
    Z76 = szeros(T,7,6)
    E3 = SMatrix{3,3,T,9}(I)
    E4 = SMatrix{3,3,T,9}(I)

    # Position
    AposT = [I Z3] # TODO is there a UniformScaling way for this instead of E?

        # This calculates the ϵ for q⊗Δq = q⊗(1 ϵᵀ)ᵀ
    AposR = [I Z43]

    # Velocity
    AvelT = [Z3 I*body.m/Δt]

    J = body.J
    ω2 = state.ωsol[2]
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)

    ω2func = skewplusdiag(ω2, sq2) * J - J * ω2 * (ω2' / sq2) - skew(J * ω2)
    AvelR = [Z34 ω2func]


    AT = [[AposT;AvelT] Z67]
    AR = [Z76 [AposR;AvelR]]

    Dinv = [
        E3 Z3 Z34                   Z3
        Z3 E3 Z34                   Z3
        Z3 Z3 VLᵀmat(state.qsol[2]) Z3
        Z3 Z3 Z34                   E3
    ]

    return [AT;AR], Dinv
end