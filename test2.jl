using StaticArrays
using BenchmarkTools

function testa(sγ1::SVector{N,T},Δsγ::SVector{N,T},τ,αmax) where {N,T}
    for i=1:N
        temp = τ*sγ1[i]/Δsγ[i]
        (temp > 0) && (temp < αmax) && (αmax = temp)
    end

    return αmax
end

a = @SVector rand(3)
b = @SVector rand(3)
c = .995
d = 1.