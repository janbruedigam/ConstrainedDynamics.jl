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
    cf = impact.cf

    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = Float64[1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    ga1 = ineq.ga1
    sl1 = ineq.sl1

    Dv = D*s1
    if norm(Dv) < 1e-10
        return Nx'*(ga1/sl1*φ - μ/sl1)
    else
        return Nx'*(ga1/sl1*φ - μ/sl1) + D'*D/norm(Dv)*s1*cf*(-ga1/sl1*φ+μ/sl1)
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

    Dv = D*s1
    if norm(Dv) < 1e-10
        return Nx'*ga1/sl1*Nv
    else
        return Nx'*ga1/sl1*Nv + cf*ga1/(norm(Dv))*(D'*D - s1*s1'*D'*D/(s1'*D'*D*s1) - D'*D*s1*Nv/sl1)
    end
end
