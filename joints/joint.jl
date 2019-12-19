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

g(J::Joint{T,Nc}) where {T,Nc} = @SVector zeros(T,Nc)

zeroBlock(J::Joint{T,Nc}) where {T,Nc} = @SMatrix zeros(T,Nc,6)
∂g∂posa(J::Joint) = zeroBlock(J)
∂g∂posb(J::Joint) = zeroBlock(J)
∂g∂vela(J::Joint) = zeroBlock(J)
∂g∂velb(J::Joint) = zeroBlock(J)
