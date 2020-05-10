abstract type Joint{T,N} end

Base.show(io::IO, joint::Joint) = summary(io, joint)

getT(joint::Joint{T}) where T = T
getN(joint::Joint{T,N}) where {T,N} = N

@inline setForce!(joint::Joint) = return
@inline setForce!(joint::Joint, body1::AbstractBody, body2::AbstractBody, Fτ::SVector{0,T}, No) where T = setForce!(joint)
@inline function clearForce!(joint::Joint, body2::Body, No)
    body2.F[No] -= joint.F[2]
    body2.τ[No] -= joint.τ[2]
    return
end
@inline function clearForce!(joint::Joint, body1::Body, body2::Body, No)
    body1.F[No] -= joint.F[1]
    body1.τ[No] -= joint.τ[1]
    body2.F[No] -= joint.F[2]
    body2.τ[No] -= joint.τ[2]
    return
end
@inline function updateForce!(joint::Joint, body2::Body, F2, τ2, No)
    body2.F[No] += F2
    body2.τ[No] += τ2

    joint.F[2] = F2
    joint.τ[2] = τ2
    return
end
@inline function updateForce!(joint::Joint, body1::Body, body2::Body, F1, τ1, F2, τ2, No)
    body1.F[No] += F1
    body1.τ[No] += τ1
    body2.F[No] += F2
    body2.τ[No] += τ2

    joint.F[1] = F1
    joint.τ[1] = τ1
    joint.F[2] = F2
    joint.τ[2] = τ2
    return
end

@inline minimalCoordinates(joint::Joint{T,N}) where {T,N} = @SVector zeros(T, 3 - N)
@inline g(joint::Joint{T,N}) where {T,N} = @SVector zeros(T, N)

@inline ∂g∂posa(joint::Joint{T,N}) where {T,N} = @SMatrix zeros(T, N, 6)
@inline ∂g∂posb(joint::Joint{T,N}) where {T,N} = @SMatrix zeros(T, N, 6)
@inline ∂g∂vela(joint::Joint{T,N}) where {T,N} = @SMatrix zeros(T, N, 6)
@inline ∂g∂velb(joint::Joint{T,N}) where {T,N} = @SMatrix zeros(T, N, 6)
@inline ∂g∂con(joint::Joint{T,N}) where {T,N} = @SMatrix zeros(T, N, 6)
