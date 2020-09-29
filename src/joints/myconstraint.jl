using StaticArrays
#using ForwardDiff 
#using ConstrainedDynamics: AbstractBody,UnitQuaternion
mutable struct myconstraint{T} #<: Joint{T,N}

    F::SVector{3,T}
    τ::SVector{3,T}
    childid::Int64

    function myconstraint{T}(body1::AbstractBody, body2::AbstractBody) where T

        F = zeros(T,3)
        τ = zeros(T,3)
        childid = body2.id

        new{T}(F, τ, childid), body1.id, body2.id

    end

end 

@inline function g(joint::myconstraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return #=write the suitable constraint function and return a vector in case the function g:n=>1:
if length(g)==1 return [g] else return g=#
    
end

@inline function g(joint::myconstraint, xb::AbstractVector, qb::UnitQuaternion)
    return #write the suitable constraint function and return a vector in case the function g:n=>1 
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

@inline function ∂g∂posa(joint::myconstraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    global _joint=joint
    x=[xa,params(qa),xb,params(qb)]
    D= ForwardDiff.jacobian(g2,x)[:,1:7]#*[Matrix{Float64}(I, 3, 3) zeros(3,3); zeros(4,3) LVᵀmat(qa) ] for rotationalderivative
    return D
end
@inline function ∂g∂posb(joint::myconstraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    global _joint=joint
    x=[xa,params(qa),xb,params(qb)]
    D= ForwardDiff.jacobian(g2,x)[:,8:14]#*[Matrix{Float64}(I, 3, 3) zeros(3,3); zeros(4,3) LVᵀmat(qb) ] for rotationalderivative
    return D
end
@inline function ∂g∂posb(joint::myconstraint, xb::AbstractVector, qb::UnitQuaternion)
    global _joint=joint
    x=[xb,params(qb)]
    D= ForwardDiff.jacobian(g1,x)#*[Matrix{Float64}(I, 3, 3) zeros(3,3); zeros(4,3) LVᵀmat(qb) ] for rotationalderivative
    return D  
end
