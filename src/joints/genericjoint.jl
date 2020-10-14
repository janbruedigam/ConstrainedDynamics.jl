mutable struct GenericJoint{T,N} <: AbstractJoint{T,N}
    F::SVector{3,T}
    τ::SVector{3,T}

    childid::Int64
    
    g::Function  
        
    
    function GenericJoint{N}(body1::AbstractBody{T}, body2::AbstractBody,gfunc) where {T,N}
        F = zeros(T,3)
        τ = zeros(T,3)
        childid = body2.id

        new{T,N}(F, τ, childid, gfunc), body1.id, body2.id
    end

end 

@inline constraintmat(::GenericJoint) = I

g(joint::GenericJoint, statea::State, stateb::State, Δt) = joint.g(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...)
g(joint::GenericJoint, stateb::State, Δt) = joint.g(joint, posargsnext(stateb, Δt)...)
g(joint::GenericJoint, statea::State, stateb::State) = joint.g(joint, posargsc(statea)..., posargsc(stateb)...)
g(joint::GenericJoint, stateb::State) = joint.g(joint, posargsc(stateb)...)


@inline function g2(x)
    xa=x[1:3]
    qa=UnitQuaternion(x[4:7]...,false)
    xb=x[8:10]
    qb=UnitQuaternion(x[11:14]...,false)
    return _joint.g(_joint,xa,qa,xb,qb) ####
end
@inline function g1(x)
    xb=x[1:3]
    qb=UnitQuaternion(x[4:7]...,false)
    return _joint.g(_joint,xb,qb)       #### 
end

@inline function ∂g∂posa(joint::GenericJoint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    global _joint=joint
    x=[xa;params(qa);xb;params(qb)]
    X= ForwardDiff.jacobian(g2,x)[:,1:3]
    Q= ForwardDiff.jacobian(g2,x)[:,4:7]
    return X,Q
end

@inline function ∂g∂posb(joint::GenericJoint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    global _joint=joint
    x=[xa;params(qa);xb;params(qb)]
    X= ForwardDiff.jacobian(g2,x)[:,8:10]
    Q= ForwardDiff.jacobian(g2,x)[:,11:14]
    return X,Q
end

@inline function ∂g∂posb(joint::GenericJoint, xb::AbstractVector, qb::UnitQuaternion)
    global _joint=joint
    x=[xb;params(qb)]
    X= ForwardDiff.jacobian(g1,x)[:,1:3]
    Q= ForwardDiff.jacobian(g1,x)[:,4:7]
    return X,Q 
end