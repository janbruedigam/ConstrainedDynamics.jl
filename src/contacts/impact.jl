mutable struct Impact{T}

    function Impact(body::Body{T}) where T
        new{T}()
    end
end


function g(ineq,impact::Impact,body::Body,dt,No,μ)
    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    ga1 = ineq.ga1
    sl1 = ineq.sl1
    Σ = ga1/sl1
    φ = body.x[No][3]+dt*body.s1[3]

    Nx'*(Σ*φ - ga1 - μ/sl1)
end

function diagval(ineq,impact::Impact,body::Body,dt)
    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    Nv = dt*Nx
    ga1 = ineq.ga1
    sl1 = ineq.sl1
    Σ = ga1/sl1

    Nx'*Σ*Nv
end
