mutable struct Impact{T}
    Nx::Adjoint{Float64,SVector{6,Float64}}
    # D::SMatrix{2,6,T,12}
    offset::SVector{6,T}


    function Impact(body::Body{T},normal::AbstractVector{T};offset::AbstractVector{T}=zeros(3)) where T
        normal = normal/norm(normal)

        A = Array(svd(skew(normal)).V)
        A[:,3] = normal # to ensure correct sign
        Ainv = inv(A)
        ainv3 = Ainv[3,:]
        Nx = [ainv3;0;0;0]'
        # D = [A[:,1:2];zeros(3,2)]'

        new{T}(Nx,[offset;0;0;0]), body
    end
end


function g(impact::Impact,body::Body,dt,No)
    impact.Nx[SVector(1,2,3)]'*(getx3(body,dt)-impact.offset[SVector(1,2,3)])
end

function dynineq(ineq,impact::Impact,body::Body,dt,No,μ)
    φ = g(impact,body,dt,No)

    γ1 = ineq.γ1[1]
    s1 = ineq.s1[1]

    return impact.Nx'*(γ1/s1*φ - μ/s1)
end

function diagval(ineq,impact::Impact,body::Body,dt)
    Nx = impact.Nx
    Nv = dt*Nx

    γ1 = ineq.γ1[1]
    s1 = ineq.s1[1]

    return Nx'*γ1/s1*Nv
end
