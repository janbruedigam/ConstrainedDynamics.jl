struct Storage{T}
    x::Vector{Vector{SVector{3,T}}}
    q::Vector{Vector{Quaternion{T}}}

    function Storage{T}(steps,nlinks) where T
        x = [[(@SVector zeros(3)) for i=1:steps] for j=1:nlinks]
        q = [[Quaternion{T}() for i=1:steps] for j=1:nlinks]
        new{T}(x,q)
    end
end
