mutable struct Friction{T}
    cf::Float64

    function Friction(body::Body{T},cf) where T
        new{T}(cf)
    end
end


function g(ineqc,friction::Friction,body::Body,dt,No)
    body.x[No][3]+dt*body.s1[3]
end

function dynineq(ineqc,friction::Friction,body::Body,dt,No,μ)
    φ = g(ineqc,friction,body,dt,No)
    cf = friction.cf

    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    γ1 = ineqc.γ1[1]
    sl1 = ineqc.s1[1]

    Dv = D*s1

    ezg = SVector{3,Float64}(0,0,-mechanism.g)
    b = D[:,1:3]*(body.m*(( - getv1(body,dt))/dt + ezg) - body.F[2])

    if norm(b)>0
        b = b/norm(b)*minimum([norm(b);cf*γ1])
    end

    if norm(b) < cf*γ1
        return Nx'*(γ1/sl1*φ - μ/sl1)
    else
        X = D*(body.m*[getv1(body,dt);getω1(body,dt)]/dt - [body.F[2];body.τ[2]])
        if norm(X) == 0
            throw("should not happen")
        end
        return (Nx'- D'*X/norm(X)*cf)*(γ1/sl1*φ - μ/sl1)
    end


    # return Nx'*(γ1/sl1*φ - μ/sl1)
    # if norm(Dv) < γ1/(dt*g)
    #     return Nx'*(γ1/sl1*φ - μ/sl1)
    # else
    #     return Nx'*(γ1/sl1*φ - μ/sl1) + D'*D/norm(Dv)*s1*cf*(-γ1/sl1*φ+μ/sl1)
    # end
end

function diagval(ineqc,friction::Friction,body::Body,dt)
    No = 2
    g = 9.81
    cf = friction.cf

    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    γ1 = ineqc.γ1[1]
    sl1 = ineqc.s1[1]

    Dv = D*s1

    ezg = SVector{3,Float64}(0,0,g)
    b = D[:,1:3]*(body.m*(( - getv1(body,dt))/dt + ezg) - body.F[2])

    if norm(b)>0
        b = b/norm(b)*minimum([norm(b);cf*γ1])
    end

    if norm(b) < cf*γ1
        return Nx'*γ1/sl1*Nv
    else
        X = D*(body.m*[getv1(body,dt);getω1(body,dt)]/dt - [body.F[2];body.τ[2]])
        if norm(X) == 0
            throw("should not happen")
        end
        return (Nx'-D'*X/norm(X)*cf)*γ1/sl1*Nv
    end



    # return Nx'*γ1/sl1*Nv
    # if norm(Dv) < γ1/(dt*g)
    #     return Nx'*γ1/sl1*Nv + dt*cf*g*D'*D
    # else
    #     return Nx'*γ1/sl1*Nv + cf*γ1/(norm(Dv))*(D'*D - s1*s1'*D'*D/(s1'*D'*D*s1) - D'*D*s1*Nv/sl1)
    # end
end
