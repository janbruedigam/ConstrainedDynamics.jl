struct Storage{T}
    x::Vector{Vector{SVector{3,T}}}
    q::Vector{Vector{Quaternion{T}}}
    λ::Vector{Vector{SVector}}

    function Storage{T}(steps,nbodies,nconstraints) where T
        x = [[(@SVector zeros(3)) for i=1:steps] for j=1:nbodies]
        q = [[Quaternion{T}() for i=1:steps] for j=1:nbodies]
        λ = [[(@SVector zeros(1)) for i=1:steps] for j=1:nconstraints]
        new{T}(x,q,λ)
    end
end
