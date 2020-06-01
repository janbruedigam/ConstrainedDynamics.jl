mutable struct Body{T} <: AbstractBody{T}
    id::Int64
    name::String

    m::T
    J::SMatrix{3,3,T,9}

    state::State{T}

    # Solver variables
    f::SVector{6,T}


    function Body(m::T, J::AbstractArray{T,2}; name::String="") where T
        state = State{T}()

        f = zeros(T, 6)

        new{T}(getGlobalID(), name, m, J, state, f)
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

    Body(T, contents...)
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

@inline getM(body::Body{T}) where T = [[SMatrix{3,3,T,9}(I*body.m);@SMatrix zeros(T,3,3)] [@SMatrix zeros(T,3,3);body.J]]

@inline function ∂zp1∂z(mechanism, body::Body{T}, xc, vc, Fk, qc, ωc, τk, Δt) where T
    stateold, fold = settempvars!(body, xc, vc, Fk, qc, ωc, τk, zeros(T,6))

    discretizestate!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    newton!(mechanism)

    A, B = ∂zp1∂z(mechanism, body::Body{T}, Δt)
    
    resettempvars!(body, stateold, fold)

    return A, B
end
