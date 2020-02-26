# mutable struct Friction{T}
#     cf::Float64

#     function Friction(body::Body{T},cf) where T
#         new{T}(cf)
#     end
# end

mutable struct Friction{T} <: Contact{T}
    Nx::Adjoint{T,SVector{6,T}}
    D::SMatrix{2,6,T,12}
    cf::T
    offset::SVector{6,T}


    function Friction(body::Body{T},normal::AbstractVector{T},cf::T;offset::AbstractVector{T}=zeros(3)) where T
        normal = normal/norm(normal)

        A = Array(svd(skew(normal)).V)
        A[:,3] = normal # to ensure correct sign
        Ainv = inv(A)
        ainv3 = Ainv[3,:]
        Nx = [ainv3;0;0;0]'
        D = [A[:,1:2];zeros(3,2)]'
        # D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]

        new{T}(Nx,D,cf,[offset;0;0;0]), body.id
    end
end


function g(friction::Friction,body::Body,dt,No)
    friction.Nx[SVector(1,2,3)]'*(getx3(body,dt)-friction.offset[SVector(1,2,3)])
end

∂g∂pos(friction::Friction) = friction.Nx
∂g∂vel(friction::Friction,dt) = friction.Nx*dt
function extrafriction!(ineq,friction::Friction{T},i,body,mechanism) where T
    dt = mechanism.dt
    cf = friction.cf
    γ1 = ineq.γ1[i]

    v = body.s1
    D = friction.D

    f = body.f
    body.s1 = @SVector zeros(T,6)
    b = D*dynamics0(body,mechanism)
    body.s1 = v
    body.f = f


    if norm(b)>0
        b = b/norm(b)*minimum([norm(b);cf*γ1])
    end
    
    body.f -= D'*b
    return
end

function schurf(ineq,friction::Friction,i,body::Body,mechanism)
    dt = mechanism.dt
    μ = mechanism.μ
    No = mechanism.No
    φ = g(friction,body,dt,No)
    cf = friction.cf

    v = body.s1
    γ1 = ineq.γ1[i]
    s1 = ineq.s1[i]

    Nx = friction.Nx
    D = friction.D
    Dv = friction.D*v

    f = body.f
    body.s1 = @SVector zeros(6)
    b = D*dynamics0(body,mechanism)
    body.s1 = v
    body.f = f


    if norm(b)>0
        b = b/norm(b)*minimum([norm(b);cf*γ1])
    end

    if norm(b) < cf*γ1
        return Nx'*(γ1/s1*φ - μ/s1)
    else
        X = D*(body.m*[getv1(body,dt);getω1(body,dt)]/dt - [body.F[2];body.τ[2]])
        if norm(X) == 0
            throw("should not happen")
        end
        return (Nx'- D'*X/norm(X)*cf)*(γ1/s1*φ - μ/s1)
    end
end

function schurD(ineq,friction::Friction,i,body::Body,mechanism)
    dt = mechanism.dt
    cf = friction.cf

    v = body.s1
    γ1 = ineq.γ1[i]
    s1 = ineq.s1[i]

    Nx = friction.Nx
    Nv = dt*Nx
    D = friction.D
    Dv = friction.D*v
    
    f = body.f
    body.s1 = @SVector zeros(6)
    b = D*dynamics0(body,mechanism)
    body.s1 = v
    body.f = f

    if norm(b)>0
        b = b/norm(b)*minimum([norm(b);cf*γ1])
    end

    if norm(b) < cf*γ1
        return Nx'*γ1/s1*Nv
    else
        X = D*(body.m*[getv1(body,dt);getω1(body,dt)]/dt - [body.F[2];body.τ[2]])
        if norm(X) == 0
            throw("should not happen")
        end
        return (Nx'-D'*X/norm(X)*cf)*γ1/s1*Nv
    end
end

