abstract type Joint{T,Nc} end

mutable struct XMLJoint{T}
    # From joint tag
    name::String
    type::String
    x::Vector{T}
    q::Quaternion{T}
    axis::Vector{T}
    parent::String
    child::String
    parentid::Int64
    childid::Int64
    pids::Tuple{Int64,Int64}

    dof::Int64
end


function Base.getproperty(C::Joint{T,Nc},x::Symbol) where {T,Nc}
    if x == :T
        return T
    elseif x == :Nc
        return Nc
    else
        return getfield(C,x)
    end
end

@inline g(C::Joint{T,Nc}) where {T,Nc} = @SVector zeros(T,Nc)

@inline zeroBlock(C::Joint{T,Nc}) where {T,Nc} = @SMatrix zeros(T,Nc,6)
@inline ∂g∂posa(C::Joint) = zeroBlock(C)
@inline ∂g∂posb(C::Joint) = zeroBlock(C)
@inline ∂g∂vela(C::Joint) = zeroBlock(C)
@inline ∂g∂velb(C::Joint) = zeroBlock(C)
