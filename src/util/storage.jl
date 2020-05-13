struct Storage{T,N}
    x::Vector{Vector{SVector{3,T}}}
    q::Vector{Vector{Quaternion{T}}}
    v::Vector{Vector{SVector{3,T}}}
    ω::Vector{Vector{SVector{3,T}}}

    function Storage{T}(steps, nbodies) where T
        x = [[(@SVector zeros(3)) for i = steps] for j = 1:nbodies]
        q = [[Quaternion{T}() for i = steps] for j = 1:nbodies]
        v = [[(@SVector zeros(3)) for i = steps] for j = 1:nbodies]
        ω = [[(@SVector zeros(3)) for i = steps] for j = 1:nbodies]
        new{T,length(steps)}(x, q, v, ω)
    end

    Storage{T}() where T = Storage{T}(Base.OneTo(0),0)
end
