using ForwardDiff
using LinearAlgebra

function f(x)
    x[1]^2*x[2] + x[2]^2 + x[3]^2 + x[4]^2
end

function ce(x)
    [
        x[1] - 2.
        x[2] - 3.
    ]
end

function ci(x)
    [
        x[3]-2.
        x[4]-3.
    ]
end


function solve(f,ce,ci,x)
    nx = length(x)
    ns = length(ci(x))
    ny = length(ce(x))
    nz = ns

    y = zeros(ny)
    z = zeros(nz)
    s = rand(ns)
    In = Matrix{Float64}(I,ns,ns)
    e = ones(ns)
    μ = 1.
    τ = 0.995
    σ = 0.2

    for iter=1:100
        L(x) = f(x) - y'*ce(x) - z'*(ci(x)-s)

        Ae = ForwardDiff.jacobian(ce,x)
        Ai = ForwardDiff.jacobian(ci,x)
        ∇f = ForwardDiff.gradient(f,x)
        Z = diagm(z)
        S = diagm(s)
        ∇²f = ForwardDiff.hessian(L,x)

        Mat = [
            ∇²f          zeros(nx,ns) -Ae'        -Ai'
            zeros(ns,nx) S\Z          zeros(ns,ny) In
            Ae           zeros(ny,ns) zeros(ny,ny) zeros(ny,nz)
            Ai           -In          zeros(nz,ny) zeros(nz,nz)
        ]

        vec = -[
            ∇f - Ae'*y - Ai'*z
            z-inv(S)*μ
            ce(x)
            ci(x) - s
        ]

        sol = Mat\vec
        px = sol[1:nx]
        ps = sol[nx+1:nx+ns]
        py = sol[nx+ns+1:nx+ns+ny]
        pz = sol[nx+ns+ny+1:nx+ns+ny+nz]

        αs = ones(ns)
        αsmax = 1.
        αz = ones(nz)
        αzmax = 1.

        for i=1:ns
            if ps[i] < 0
                αs[i] = minimum([1.;-τ*s[i]/ps[i]])
            end
            αsmax = minimum(αs)

            if pz[i] < 0
                αz[i] = minimum([1.;-τ*z[i]/pz[i]])
            end
            αzmax = minimum(αz)
        end

        x += αsmax.*px
        s += αsmax.*ps
        y += αzmax.*py
        z += αzmax.*pz

        norm1 = vec[1:nx]
        norm2 = vec[nx+1:nx+ns]
        norm3 = vec[nx+ns+1:nx+ns+ny]
        norm4 = vec[nx+ns+ny+1:nx+ns+ny+nz]

        E = maximum([norm1;norm2;norm3;norm4])
        μ = σ*μ

        if E<μ
            display(iter)
            break
        end
    end

    return x

end
