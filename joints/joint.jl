abstract type Joint{T,Nc} end

function Base.getproperty(J::Joint{T,Nc},x::Symbol) where {T,Nc}
    if x == :T
        return T
    elseif x == :Nc
        return Nc
    else
        return getfield(J,x)
    end
end

@inline g(J::Joint{T,Nc}) where {T,Nc} = @SVector zeros(T,Nc)

@inline zeroBlock(J::Joint{T,Nc}) where {T,Nc} = @SMatrix zeros(T,Nc,6)
@inline ∂g∂posa(J::Joint) = zeroBlock(J)
@inline ∂g∂posb(J::Joint) = zeroBlock(J)
@inline ∂g∂vela(J::Joint) = zeroBlock(J)
@inline ∂g∂velb(J::Joint) = zeroBlock(J)
