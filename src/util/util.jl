function deleteat(M::Array{T,2}, i1::Integer, i2::Integer) where T
    return [M[1:i1 - 1,1:i2 - 1] M[1:i1 - 1,i2 + 1:end];M[i1 + 1:end,1:i2 - 1] M[i1 + 1:end,i2 + 1:end]]
end

deleteat(M::Array{T,2},i::Integer) where T = deleteat(M, i, i)

@inline function skewplusdiag(v::AbstractVector{T},w::T) where T
    SA[
         w    -v[3]  v[2]
         v[3]  w    -v[1]
        -v[2]  v[1]  w
    ]
end

function getfieldnumber(obj)
    i = 1
    while true
        !isdefined(obj, i) ? break : (i+=1)
    end
    return i-1
end
