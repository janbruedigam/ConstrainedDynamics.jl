mutable struct Impact{T}
    cf::Float64

    function Impact(body::Body{T},cf) where T
        new{T}(cf)
    end
end


function g(ineq,impact::Impact,body::Body,dt,No)
    φ = body.x[No][3]+dt*body.s1[3]
    f1 = φ - ineq.sl1

    cf = impact.cf
    D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]
    Dv = D*body.s1
    f2 = 0.
    if Dv == zeros(2)
        f2 = ineq.b1
    else
        f2 = cf*ineq.ga1*Dv + ineq.b1*norm(Dv)
    end

    [f1;f2]
end

function dynineq(ineq,impact::Impact,body::Body,dt,No,μ)
    φ = body.x[No][3]+dt*body.s1[3]
    cf = impact.cf

    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    ga1 = ineq.ga1
    sl1 = ineq.sl1

    Dv = D*s1
    if Dv == zeros(2)
        return Nx'*(ga1/sl1*φ - μ/sl1)
    else
        return Nx'*(ga1/sl1*φ - μ/sl1) + cf*D'*D/norm(Dv)*s1*(ga1-ga1/sl1*φ+μ/sl1)
    end
end

function diagval(ineq,impact::Impact,body::Body,dt)
    No = 2
    cf = impact.cf

    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    ga1 = ineq.ga1
    sl1 = ineq.sl1
    b1 = ineq.b1

    Dv = D*s1
    if Dv == zeros(2)
        return Nx'*ga1/sl1*Nv
    else
        X = cf*ga1*D + b1*s1'*D'*D/norm(Dv)
        return D'*X/norm(Dv) - cf*D'*D/norm(Dv)*s1*ga1/sl1*Nv + Nx'*ga1/sl1*Nv
    end
end
