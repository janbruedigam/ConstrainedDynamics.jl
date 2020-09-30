using StaticArrays
#using ForwardDiff 
using ConstrainedDynamics:GenericJoint

mutable struct my_constraint{T,N} <: GenericJoint{T,N}

    F::SVector{3,T}
    τ::SVector{3,T}
    childid::Int64

    function my_constraint{T,N}(body1::AbstractBody, body2::AbstractBody) where {T,N}

        F = zeros(T,3)
        τ = zeros(T,3)
        childid = body2.id

        new{T,N}(F, τ, childid), body1.id, body2.id

    end

end 
@inline function my_Joint(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T})) where T  
    global vertices = (p1, p2)
    global V1, V2, V3 = orthogonalrows(axis)
    global V12 = [V1;V2]
return my_constraint{T,2}(body1, body2)
end 
@inline function g(joint::my_constraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)


    G= vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
    if length(G)==1 
        return [G] 
    else 
        return G 
    end 
    
end

@inline function g(joint::my_constraint, xb::AbstractVector, qb::UnitQuaternion)
    
    
    G= xb + vrotate(vertices[2], qb) - vertices[1]
    if length(G)==1 
        return [G] 
    else 
        return G 
    end 
end

@inline function g2(x)
    xa=x[1:3]
    qa=UnitQuaternion(x[4:7]...,false)
    xb=x[8:10]
    qb=UnitQuaternion(x[11:14]...,false)
    return g(_joint,xa,qa,xb,qb)
end
@inline function g1(x)
    xb=x[1:3]
    qb=UnitQuaternion(x[4:7]...,false)
    return g(_joint,xb,qb)
end

@inline function ∂g∂posa(joint::my_constraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    global _joint=joint
    x=[xa,params(qa),xb,params(qb)]
    D= ForwardDiff.jacobian(g2,x)[:,1:7]#*[Matrix{Float64}(I, 3, 3) zeros(3,3); zeros(4,3) LVᵀmat(qa) ] for rotationalderivative
    return D
end
@inline function ∂g∂posb(joint::my_constraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    global _joint=joint
    x=[xa,params(qa),xb,params(qb)]
    D= ForwardDiff.jacobian(g2,x)[:,8:14]#*[Matrix{Float64}(I, 3, 3) zeros(3,3); zeros(4,3) LVᵀmat(qb) ] for rotationalderivative
    return D
end
@inline function ∂g∂posb(joint::my_constraint, xb::AbstractVector, qb::UnitQuaternion)
    global _joint=joint
    x=[xb,params(qb)]
    D= ForwardDiff.jacobian(g1,x)#*[Matrix{Float64}(I, 3, 3) zeros(3,3); zeros(4,3) LVᵀmat(qb) ] for rotationalderivative
    return D  
end
