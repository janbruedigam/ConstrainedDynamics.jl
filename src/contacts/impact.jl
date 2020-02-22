mutable struct Impact{T}
    cf::Float64

    function Impact(body::Body{T},cf) where T
        new{T}(cf)
    end
end


function g(ineq,impact::Impact,body::Body,dt,No)
    φ = body.x[No][3]+dt*body.s1[3]
    φ - ineq.sl1
end

function dynineq(ineq,impact::Impact,body::Body,dt,No,μ)
    φ = body.x[No][3]+dt*body.s1[3]
    g = 9.81
    cf = impact.cf

    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    ga1 = ineq.ga1
    sl1 = ineq.sl1

    Dv = D*s1

    ezg = SVector{3,Float64}(0,0,g)
    b = D[:,1:3]*(body.m*(( - getv1(body,dt))/dt + ezg) - body.F[2])

    if norm(b)>0
        b = b/norm(b)*minimum([norm(b);cf*ga1])
    end

    if norm(b) < cf*ga1
        return Nx'*(ga1/sl1*φ - μ/sl1)
    else
        X = D*(body.m*[getv1(body,dt);getω1(body,dt)]/dt - [body.F[2];body.τ[2]])
        if norm(X) == 0
            throw("should not happen")
        end
        return (Nx'- D'*X/norm(X)*cf)*(ga1/sl1*φ - μ/sl1)
    end


    # return Nx'*(ga1/sl1*φ - μ/sl1)
    # if norm(Dv) < ga1/(dt*g)
    #     return Nx'*(ga1/sl1*φ - μ/sl1)
    # else
    #     return Nx'*(ga1/sl1*φ - μ/sl1) + D'*D/norm(Dv)*s1*cf*(-ga1/sl1*φ+μ/sl1)
    # end
end

function diagval(ineq,impact::Impact,body::Body,dt)
    No = 2
    g = 9.81
    cf = impact.cf

    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    ga1 = ineq.ga1
    sl1 = ineq.sl1

    Dv = D*s1

    ezg = SVector{3,Float64}(0,0,g)
    b = D[:,1:3]*(body.m*(( - getv1(body,dt))/dt + ezg) - body.F[2])

    if norm(b)>0
        b = b/norm(b)*minimum([norm(b);cf*ga1])
    end

    if norm(b) < cf*ga1
        return Nx'*ga1/sl1*Nv
    else
        X = D*(body.m*[getv1(body,dt);getω1(body,dt)]/dt - [body.F[2];body.τ[2]])
        if norm(X) == 0
            throw("should not happen")
        end
        return (Nx'-D'*X/norm(X)*cf)*ga1/sl1*Nv
    end



    # return Nx'*ga1/sl1*Nv
    # if norm(Dv) < ga1/(dt*g)
    #     return Nx'*ga1/sl1*Nv + dt*cf*g*D'*D
    # else
    #     return Nx'*ga1/sl1*Nv + cf*ga1/(norm(Dv))*(D'*D - s1*s1'*D'*D/(s1'*D'*D*s1) - D'*D*s1*Nv/sl1)
    # end
end
