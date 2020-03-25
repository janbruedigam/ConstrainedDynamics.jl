using StaticArrays
using BenchmarkTools
using LinearAlgebra

@generated function submatrix(A::SMatrix{N1,N2,T,N1N2}, r, c) where {N1,N2,T,N1N2}
    inds = [:(($i < r ? $i : $i + 1) + (($j < c ? $j : $j + 1) - 1) * N1) for j in 1:(N2 - 1) for i in 1:(N1 - 1)]
    :(SMatrix{N1 - 1,N2 - 1,T,(N1 - 1) * (N2 - 1)}(A[SVector($(inds...))]))
end

@inline function StaticArrays._det(::Size{(3, 3)}, A::StaticMatrix{N,N,T}) where {N,T}
    cof11 =  det(submatrix(A, 1, 1))
    cof21 = -det(submatrix(A, 2, 1))
    cof31 =  det(submatrix(A, 3, 1))
    return A[1] * cof11 + A[2] * cof21 + A[3] * cof31
end

@inline function StaticArrays._det(::Size{(5, 5)}, A::StaticMatrix{N,N,T}) where {N,T}
    cof11 =  det(submatrix(A, 1, 1))
    cof21 = -det(submatrix(A, 2, 1))
    cof31 =  det(submatrix(A, 3, 1))
    cof41 = -det(submatrix(A, 4, 1))
    cof51 =  det(submatrix(A, 5, 1))
    return A[1] * cof11 + A[2] * cof21 + A[3] * cof31 + A[4] * cof41 + A[5] * cof51
end

@inline function StaticArrays._det(::Size{(6, 6)}, A::StaticMatrix{N,N,T}) where {N,T}
    cof11 =  det(submatrix(A, 1, 1))
    cof21 = -det(submatrix(A, 2, 1))
    cof31 =  det(submatrix(A, 3, 1))
    cof41 = -det(submatrix(A, 4, 1))
    cof51 =  det(submatrix(A, 5, 1))
    cof61 =  -det(submatrix(A, 6, 1))
    return A[1] * cof11 + A[2] * cof21 + A[3] * cof31 + A[4] * cof41 + A[5] * cof51 + A[6] * cof61
end

@inline function StaticArrays._inv(::Size{(3, 3)}, A)
    adj11 =  det(submatrix(A, 1, 1))
    adj12 = -det(submatrix(A, 2, 1))
    adj13 =  det(submatrix(A, 3, 1))
    adj21 = -det(submatrix(A, 1, 2))
    adj22 =  det(submatrix(A, 2, 2))
    adj23 = -det(submatrix(A, 3, 2))
    adj31 =  det(submatrix(A, 1, 3))
    adj32 = -det(submatrix(A, 2, 3))
    adj33 =  det(submatrix(A, 3, 3))

    invdet = 1 / (A[1] * adj11 + A[2] * adj12 + A[3] * adj13)

    B =  @SMatrix [
        adj11 * invdet
    	adj21 * invdet
        adj31 * invdet
        adj12 * invdet
        adj22 * invdet
        adj32 * invdet
        adj13 * invdet
        adj23 * invdet
        adj33 * invdet]
    return similar_type(A)(B)
end

@inline function StaticArrays._inv(::Size{(5, 5)}, A)
    adj11 =  det(submatrix(A, 1, 1))
    adj12 = -det(submatrix(A, 2, 1))
    adj13 =  det(submatrix(A, 3, 1))
    adj14 = -det(submatrix(A, 4, 1))
    adj15 =  det(submatrix(A, 5, 1))

    adj21 = -det(submatrix(A, 1, 2))
    adj22 =  det(submatrix(A, 2, 2))
    adj23 = -det(submatrix(A, 3, 2))
    adj24 =  det(submatrix(A, 4, 2))
    adj25 = -det(submatrix(A, 5, 2))

    adj31 =  det(submatrix(A, 1, 3))
    adj32 = -det(submatrix(A, 2, 3))
    adj33 =  det(submatrix(A, 3, 3))
    adj34 = -det(submatrix(A, 4, 3))
    adj35 =  det(submatrix(A, 5, 3))

    adj41 = -det(submatrix(A, 1, 4))
    adj42 =  det(submatrix(A, 2, 4))
    adj43 = -det(submatrix(A, 3, 4))
    adj44 =  det(submatrix(A, 4, 4))
    adj45 = -det(submatrix(A, 5, 4))

    adj51 =  det(submatrix(A, 1, 5))
    adj52 = -det(submatrix(A, 2, 5))
    adj53 =  det(submatrix(A, 3, 5))
    adj54 = -det(submatrix(A, 4, 5))
    adj55 =  det(submatrix(A, 5, 5))

    invdet = 1 / (A[1] * adj11 + A[2] * adj12 + A[3] * adj13 + A[4] * adj14 + A[5] * adj15)

    B =  @SMatrix [
        adj11 * invdet
    	adj21 * invdet
        adj31 * invdet
        adj41 * invdet
        adj51 * invdet
        adj12 * invdet
        adj22 * invdet
        adj32 * invdet
        adj42 * invdet
        adj52 * invdet
        adj13 * invdet
        adj23 * invdet
        adj33 * invdet
        adj43 * invdet
        adj53 * invdet
        adj14 * invdet
    	adj24 * invdet
        adj34 * invdet
        adj44 * invdet
        adj54 * invdet
        adj15 * invdet
    	adj25 * invdet
        adj35 * invdet
        adj45 * invdet
        adj55 * invdet]
    return similar_type(A)(B)
end

function StaticArrays._inv(::Size{(6, 6)}, A)
    adj11 =  det(submatrix(A, 1, 1))
    adj12 = -det(submatrix(A, 2, 1))
    adj13 =  det(submatrix(A, 3, 1))
    adj14 = -det(submatrix(A, 4, 1))
    adj15 =  det(submatrix(A, 5, 1))
    adj16 = -det(submatrix(A, 6, 1))

    adj21 =  det(submatrix(A, 1, 2))
    adj22 = -det(submatrix(A, 2, 2))
    adj23 =  det(submatrix(A, 3, 2))
    adj24 = -det(submatrix(A, 4, 2))
    adj25 =  det(submatrix(A, 5, 2))
    adj26 = -det(submatrix(A, 6, 2))

    adj31 =  det(submatrix(A, 1, 3))
    adj32 = -det(submatrix(A, 2, 3))
    adj33 =  det(submatrix(A, 3, 3))
    adj34 = -det(submatrix(A, 4, 3))
    adj35 =  det(submatrix(A, 5, 3))
    adj36 = -det(submatrix(A, 6, 3))

    adj41 =  det(submatrix(A, 1, 4))
    adj42 = -det(submatrix(A, 2, 4))
    adj43 =  det(submatrix(A, 3, 4))
    adj44 = -det(submatrix(A, 4, 4))
    adj45 =  det(submatrix(A, 5, 4))
    adj46 = -det(submatrix(A, 6, 4))

    adj51 =  det(submatrix(A, 1, 5))
    adj52 = -det(submatrix(A, 2, 5))
    adj53 =  det(submatrix(A, 3, 5))
    adj54 = -det(submatrix(A, 4, 5))
    adj55 =  det(submatrix(A, 5, 5))
    adj56 = -det(submatrix(A, 6, 5))

    adj61 =  det(submatrix(A, 1, 6))
    adj62 = -det(submatrix(A, 2, 6))
    adj63 =  det(submatrix(A, 3, 6))
    adj64 = -det(submatrix(A, 4, 6))
    adj65 =  det(submatrix(A, 5, 6))
    adj66 = -det(submatrix(A, 6, 6))

    invdet = 1 / (A[1] * adj11 + A[2] * adj12 + A[3] * adj13 + A[4] * adj14 + A[5] * adj15  + A[6] * adj16)

    B =  @SMatrix [
        adj11 * invdet
    	adj21 * invdet
        adj31 * invdet
        adj41 * invdet
        adj51 * invdet
        adj61 * invdet
        adj12 * invdet
        adj22 * invdet
        adj32 * invdet
        adj42 * invdet
        adj52 * invdet
        adj62 * invdet
        adj13 * invdet
        adj23 * invdet
        adj33 * invdet
        adj43 * invdet
        adj53 * invdet
        adj63 * invdet
        adj14 * invdet
    	adj24 * invdet
        adj34 * invdet
        adj44 * invdet
        adj54 * invdet
        adj64 * invdet
        adj15 * invdet
    	adj25 * invdet
        adj35 * invdet
        adj45 * invdet
        adj55 * invdet
        adj65 * invdet
        adj16 * invdet
    	adj26 * invdet
        adj36 * invdet
        adj46 * invdet
        adj56 * invdet
        adj66 * invdet]
    return similar_type(A)(B)
end
