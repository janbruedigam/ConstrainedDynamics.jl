mutable struct Body{T} <: AbstractBody{T}
    id::Int64
    name::String

    m::T
    J::SMatrix{3,3,T,9}

    state::State{T}

    F::Vector{SVector{3,T}}
    τ::Vector{SVector{3,T}}

    # Solver variables
    f::SVector{6,T}
    solv::SVector{3,T}
    solω::SVector{3,T}


    function Body(m::T, J::AbstractArray{T,2}; name::String="") where T
        state = State{T}()

        F = [zeros(T, 3)]
        τ = [zeros(T, 3)]

        f = zeros(T, 6)
        solv = zeros(T, 3)
        solω = zeros(T, 3)

        new{T}(getGlobalID(), name, m, J, state, F, τ, f, solv, solω)
    end

    function Body(shape::Shape; name::String="")
        body = Body(shape.m, shape.J; name=name)
        push!(shape.bodyids, body.id)
        return body
    end

    function Body(::Type{T},contents...) where T
        new{T}(getGlobalID(), contents...)
    end
end

mutable struct Origin{T} <: AbstractBody{T}
    id::Int64
    name::String

    Origin{T}(; name::String="") where T = new{T}(getGlobalID(), name)

    function Origin(shape::Shape{T}; name::String="") where T
        origin = Origin{T}(name=name)
        push!(shape.bodyids, origin.id)
        return origin
    end
end

function Base.deepcopy(b::Body{T}) where T
    contents = []
    for i = 2:getfieldnumber(b)
        push!(contents, deepcopy(getfield(b, i)))
    end

    Body(T, contents...)
end

function Base.deepcopy(b::Body{T}, shape::Shape) where T
    contents = []
    
    for i = 2:getfieldnumber(b)
        push!(contents, deepcopy(getfield(b, i)))
    end

    body = Body(T, contents...)
    @assert (body.m == shape.m && body.J == shape.J)
    push!(shape.bodyids, body.id)
    return body
end

Base.length(::Body) = 6

@inline getM(body::Body{T}) where T = [[SMatrix{3,3,T,9}(I*body.m);@SMatrix zeros(T,3,3)] [@SMatrix zeros(T,3,3);body.J]]

# @inline function settempvars(body::Body{T}, x, v, F, q, ω, τ, s0, s1, f, Δt) where T
#     xold = body.x[2]
#     vold = getv1(body,Δt)
#     Fold = body.F[2]
#     qold = body.q[2]
#     ωold = getω1(body,Δt)
#     τold = body.τ[2]
#     s0old = body.s0
#     s1old = body.s1
#     fold = body.f
    

#     body.x[2] = x
#     body.x[1] = x - v*Δt
#     body.F[2] = F
#     body.q[2] = q
#     body.q[1] = Δt/2 * (q * Quaternion(sqrt(4 / Δt^2 - dot(ω, ω)), -SVector{3,T}(ω))) # this accounts for the non-unit-quaternion inverse (q2 * conj(ωbar))
#     body.τ[2] = τ
#     body.s0 = s0
#     body.s1 = s1
#     body.f = f
#     return xold, vold, Fold, qold, ωold, τold, s0old, s1old, fold
# end

# @inline function ∂zp1∂z(mechanism, body::Body{T}, xd, vd, Fd, qd, ωd, τd, Δt) where T
#     xold, vold, Fold, qold, ωold, τold, s0old, s1old, fold = settempvars(body, xd, vd, Fd, qd, ωd, τd, [vd;ωd], [vd;ωd], zeros(T,6), Δt)

#     newton!(mechanism)

#     # Velocity
#     Z = @SMatrix zeros(T,3,3)
#     Z43 = @SMatrix zeros(T,4,3)
#     Z2 = [Z Z]
#     Z2t = Z2'
#     E = SMatrix{3,3,T,9}(I)

#     AvelT = [Z E]
#     BvelT = E*Δt/body.m

#     J = body.J
#     ω1 = getω1(body, Δt)
#     ωnew = getωnew(body)
#     sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
#     sq2 = sqrt(4 / Δt^2 - ωnew' * ωnew)

#     ω1func = skewplusdiag(-ω1, sq1) * J - J * ω1 * (ω1' / sq1) + skew(J * ω1)
#     ω2func = skewplusdiag(ωnew, sq2) * J - J * ωnew * (ωnew' / sq2) - skew(J * ωnew)
#     AvelR = [Z ω2func\ω1func]
#     BvelR = 2*E

#     # Position
#     AposT = [E Z] + AvelT*Δt
#     BposT = BvelT*Δt

#     # This calculates the ϵ for q⊗Δq = q⊗(1 ϵᵀ)ᵀ
#     AposR = VLmat(getq3(body,Δt)) * ([Rmat(ωbar(body, Δt))*LVᵀmat(qd) Z43] + Lmat(qd)*derivωbar(body, Δt)*AvelR)
#     BposR = VLmat(getq3(body,Δt)) * Lmat(qd)*derivωbar(body, Δt)*BvelR

#     AT = [[AposT;AvelT] [Z2;Z2]]
#     AR = [[Z2;Z2] [AposR;AvelR]]
#     BT = [[BposT;BvelT] Z2t]
#     BR = [Z2t [BposR;BvelR]]

#     settempvars(body, xold, vold, Fold, qold, ωold, τold, s0old, s1old, fold, Δt)
#     return [AT;AR], [BT;BR]
# end
