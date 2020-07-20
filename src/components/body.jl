mutable struct Body{T} <: AbstractBody{T}
    id::Int64
    name::String
    active::Bool

    m::T
    J::SMatrix{3,3,T,9}

    state::State{T}


    function Body(m::Real, J::AbstractArray; name::String="")
        T = promote_type(eltype.((m, J))...)
        new{T}(getGlobalID(), name, true, m, J, State{T}())
    end

    function Body(shape::Shape; name::String="")
        body = Body(shape.m, shape.J; name=name)
        push!(shape.bodyids, body.id)
        return body
    end

    function Body{T}(contents...) where T
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

    return Body{T}(contents...)
end

function Base.deepcopy(b::Body{T}, shape::Shape) where T
    contents = []
    
    for i = 2:getfieldnumber(b)
        push!(contents, deepcopy(getfield(b, i)))
    end

    body = Body{T}(contents...)
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


function ∂F∂z(body::Body{T}, Δt) where T
    state = body.state
    Z3 = szeros(T,3,3)
    Z76 = szeros(T,7,6)
    Z6 = szeros(T,6,6)
    E3 = SMatrix{3,3,T,9}(I)


    AposT = [-I -E3*Δt] # TODO is there a UniformScaling way for this instead of E3?

    AvelT = [Z3 -I*body.m/Δt]

    AposR = [-Rmat(ωbar(state.ωc, Δt))*LVᵀmat(state.qc) -Lmat(state.qc)*derivωbar(state.ωc, Δt)]

    J = body.J
    ω1 = state.ωc
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    ω1func = skewplusdiag(ω1, -sq1) * J + J * ω1 * (ω1' / sq1) - skew(J * ω1)
    AvelR = [Z3 ω1func]

    
    return [[AposT;AvelT] Z6;Z76 [AposR;AvelR]]
end

function ∂F∂u(body::Body{T}, Δt) where T
    Z3 = szeros(T,3,3)
    Z43 = szeros(T,4,3)


    BposT = [Z3 Z3] # TODO is there a UniformScaling way for this instead of E3?

    BvelT = [-I Z3]

    BposR = [Z43 Z43]

    BvelR = [Z3 -2.0*I]

    
    return [BposT;BvelT;BposR;BvelR]
end

function ∂F∂g(body::Body{T}, Δt) where T
    Z3 = szeros(T,3,3)
    Z43 = szeros(T,4,3)


    BposT = [0;0;0] # TODO is there a UniformScaling way for this instead of E3?

    BvelT = [0;0;body.m]

    BposR = [0;0;0;0]

    BvelR = [0;0;0]

    
    return [BposT;BvelT;BposR;BvelR]
end

# function ∂F∂λ()

# end

function ∂F∂fz(body::Body{T}, Δt) where T
    state = body.state
    Z3 = szeros(T,3,3)
    Z34 = szeros(T,3,4)
    Z43 = szeros(T,4,3)
    Z67 = szeros(T,6,7)
    Z76 = szeros(T,7,6)


    AposT = [I Z3]

    AvelT = [Z3 I*body.m/Δt]

    AposR = [I Z43]

    J = body.J
    ω2 = state.ωsol[2]
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
    ω2func = skewplusdiag(ω2, sq2) * J - J * ω2 * (ω2' / sq2) - skew(J * ω2)
    AvelR = [Z34 inv(ω2func)]


    return [[AposT;AvelT] Z67;Z76 [AposR;AvelR]], [I Z3 Z34 Z3;Z3 I Z34 Z3;Z3 Z3 VLᵀmat(state.qsol[2]) Z3;Z3 Z3 Z34 I]
end
