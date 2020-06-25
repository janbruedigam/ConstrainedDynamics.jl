function deleteat(M::Array, i1::Integer, i2::Integer)
    return [M[1:i1 - 1,1:i2 - 1] M[1:i1 - 1,i2 + 1:end];M[i1 + 1:end,1:i2 - 1] M[i1 + 1:end,i2 + 1:end]]
end

deleteat(M::Array,i::Integer) = deleteat(M, i, i)

function orthogonalrows(axis::AbstractVector)
    if norm(axis) > 0 
        axis = normalize(axis)
    end
    A = svd(skew(axis)).Vt
    inds = SA[1; 2; 3]
    V1 = A[1,inds]'
    V2 = A[2,inds]'
    V3 = axis' # instead of A[3,:] for correct sign: abs(axis) = abs(A[3,:])

    return V1, V2, V3
end

function orthogonalcols(axis::AbstractVector)
    V1, V2, V3 = orthogonalrows(axis)
    return V1', V2', V3'
end

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

@inline offsetrange(offset, length) = (offset-1)*length+1:offset*length
@inline offsetrange(offset, length, totallength, inneroffset) = (offset-1)*totallength+(inneroffset-1)*length+1:(offset-1)*totallength+inneroffset*length