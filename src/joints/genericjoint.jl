# User-defined joint
mutable struct GenericJoint{T,N} <: AbstractJoint{T,N}
    spring::T
    damper::T

    F::SVector{3,T}
    τ::SVector{3,T}

    childid::Int64
    
    g::Function  
        
    
    function GenericJoint{N}(body1::AbstractBody{T}, body2::AbstractBody, gfunc; spring = zero(T), damper = zero(T)) where {T,N}
        F = zeros(T,3)
        τ = zeros(T,3)
        childid = body2.id

        new{T,N}(spring, damper, F, τ, childid, gfunc), body1.id, body2.id
    end

end 


### General functions
@inline constraintmat(::GenericJoint) = I


### Constraints and derivatives
## Discrete-time position wrappers (for dynamics)
g(joint::GenericJoint, statea::State, stateb::State, Δt) = joint.g(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...)
g(joint::GenericJoint, stateb::State, Δt) = joint.g(joint, posargsnext(stateb, Δt)...)
g(joint::GenericJoint, statea::State, stateb::State) = joint.g(joint, posargsc(statea)..., posargsc(stateb)...)
g(joint::GenericJoint, stateb::State) = joint.g(joint, posargsc(stateb)...)

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::GenericJoint, xa::AbstractVector, qa::QuatRotation, xb::AbstractVector, qb::QuatRotation)
    z = [xa;params(qa);xb;params(qb)]

    G = ForwardDiff.jacobian(x -> joint.g(joint,x[1:3],QuatRotation(x[4:7]...,false),x[8:10],QuatRotation(x[11:14]...,false)), z)
    X = G[:,1:3]
    Q = G[:,4:7]

    return X, Q
end
@inline function ∂g∂posb(joint::GenericJoint, xa::AbstractVector, qa::QuatRotation, xb::AbstractVector, qb::QuatRotation)
    z = [xa;params(qa);xb;params(qb)]
    
    G = ForwardDiff.jacobian(x -> joint.g(joint,x[1:3],QuatRotation(x[4:7]...,false),x[8:10],QuatRotation(x[11:14]...,false)), z)
    X = G[:,8:10]
    Q = G[:,11:14]

    return X, Q
end
@inline function ∂g∂posb(joint::GenericJoint, xb::AbstractVector, qb::QuatRotation)
    z = [xb;params(qb)]
    
    G = ForwardDiff.jacobian(x -> joint.g(joint,x[1:3],QuatRotation(x[4:7]...,false)), z)
    X = G[:,1:3]
    Q = G[:,4:7]

    return X, Q
end
