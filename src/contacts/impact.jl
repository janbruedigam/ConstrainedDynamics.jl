mutable struct Impact{T}

    function Impact(body::Body{T}) where T
        new{T}()
    end
end


function g(ineq,impact::Impact,body::Body,dt,No)
    φ = body.x[No][3]+dt*body.s1[3]
    f1 = φ - ineq.sl1

    D = [1 0 0 0 0 0;0 1 0 0 0 0]
    cf = 0.2

    f2 = cf^2*ineq.ga1^2 - ineq.b1'*ineq.b1
    f3 = D*body.s1*ineq.ga1 + 2*ineq.psi1*ineq.ga1*ineq.b1

    [f1;f2;f3]
end

function dynineq(ineq,impact::Impact,body::Body,dt,No,μ)
    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = [1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    ga1 = ineq.ga1
    sl1 = ineq.sl1
    slf1 = ineq.slf1
    psi1 = ineq.psi1
    b1 = ineq.b1

    cf = 0.2

    Ax = [-Nx-b1'*D;zeros(6)']
    Av = [Nv;zeros(6)']
    H = [D*s1+2*psi1*b1 2*ga1*b1]
    X = [0 0;2*cf^2*ga1 0]
    B = [zeros(2)';b1']
    Z = [ga1 0;0 psi1]
    S = [sl1 0;0 slf1]

    K = X + 1/2*1/(ga1*psi1)*B*H + Z\S
    φ = body.x[No][3]+dt*body.s1[3]
    ci = [φ;cf^2*ga1^2-b1'*b1]

    1/2*ga1*D'*(1/psi1*D*s1+2*b1)-(1/2*D'*1/psi1*H-Ax')/K*(ci+1/2*1/psi1*B*D*s1+B*b1-μ*inv(Z)*ones(2))
end

function diagval(ineq,impact::Impact,body::Body,dt)
    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    D = [1 0 0 0 0 0;0 1 0 0 0 0]

    s1 = body.s1
    ga1 = ineq.ga1
    sl1 = ineq.sl1
    slf1 = ineq.slf1
    psi1 = ineq.psi1
    b1 = ineq.b1

    cf = 0.2

    Ax = [-Nx-b1'*D;zeros(6)']
    Av = [Nv;zeros(6)']
    H = [D*s1+2*psi1*b1 2*ga1*b1]
    X = [0 0;2*cf^2*ga1 0]
    B = [zeros(2)';b1']
    Z = [ga1 0;0 psi1]
    S = [sl1 0;0 slf1]

    K = X + 1/2*1/(ga1*psi1)*B*H + Z\S

    1/2*ga1/psi1*D'*D + (Ax' - 1/2*D'*1/psi1*H)/K*(Av+1/2*1/psi1*B*D)
end
