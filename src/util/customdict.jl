struct UnitDict{R,T}
    keys::R
    values::Vector{T}

    UnitDict(unitrange, values::Vector{T}) where T = new{typeof(unitrange),T}(unitrange, values)
    UnitDict(values::Vector) = UnitDict(Base.OneTo(length(values)), values)
    function UnitDict(dict::Dict{T1,T2}) where {T1,T2}
        if isempty(dict)
            UnitDict(zero(T1):zero(T1), T2[])
        else
            a = minimum(keys(dict))
            b = maximum(keys(dict))

            a == 1 ? (range = Base.OneTo(b)) : (range = a:b)
            values = T2[]
            for k in range
                push!(values, dict[k])
            end

            UnitDict(range, values)
        end
    end
end


@inline Base.length(dict::UnitDict) = length(dict.values)
@inline Base.firstindex(dict::UnitDict) = first(dict.keys)
@inline Base.lastindex(dict::UnitDict) = last(dict.keys)

@inline Base.getindex(dict::UnitDict{Base.OneTo{K},T}, key::K) where {T,K} = dict.values[key]
@inline Base.getindex(dict::UnitDict{UnitRange{K},T}, key::K) where {T,K} = dict.values[key - first(dict.keys) + 1]
@inline Base.getindex(dict::UnitDict{Base.OneTo{K},T}, key::AbstractVector{K}) where {T,K} = dict.values[key]
@inline Base.getindex(dict::UnitDict{UnitRange{K},T}, key::AbstractVector{K}) where {T,K} = dict.values[key .- (first(dict.keys) + 1)]
@inline Base.getindex(dict::UnitDict{Base.OneTo{K},T}, key::UnitRange{K}) where {T,K} = dict.values[key]
@inline Base.getindex(dict::UnitDict{UnitRange{K},T}, key::UnitRange{K}) where {T,K} = dict.values[key .- (first(dict.keys) + 1)]
@inline Base.setindex!(dict::UnitDict{Base.OneTo{K},T}, value::T, key::K) where {T,K} = setindex!(dict.values, value, key)
@inline Base.setindex!(dict::UnitDict{UnitRange{K},T}, value::T, key::K) where {T,K} = setindex!(dict.values, value, key - first(dict.keys) + 1)

@inline Base.iterate(d::UnitDict, i = 1) = i > length(d) ? nothing : (d.values[i], i + 1)
@inline Base.pairs(d::UnitDict, i = 1) = Base.Generator(=>, d.keys, d.values)

@inline Base.foreach(f, itr::UnitDict, arg...) = (for x in itr; f(x, arg...); end; return)
@inline function foreachactive(f, itr::UnitDict, arg...)
    for x in itr
        isinactive(x) && continue
        f(x, arg...)
    end
    return
end

@inline Base.haskey(d::UnitDict, key) = Base.ht_keyindex(d, key)

function Base.ht_keyindex(h::UnitDict, key)
    for el in h.keys
        if el == key
            return true
        end
    end
    return false
end
