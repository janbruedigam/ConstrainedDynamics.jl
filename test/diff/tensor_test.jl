using ConstrainedDynamics
using ConstrainedDynamics: LVᵀmat, Rotational, Translational, ∂2g∂posaa, ∂2g∂posab, ∂2g∂posba, ∂2g∂posbb, g
using ForwardDiff
using Rotations: params, QuatRotation
using LinearAlgebra


function transtest()
    xa = srand(3)
    qa = rand(QuatRotation).q
    xb = srand(3)
    qb = rand(QuatRotation).q
    z = [xa;params(qa);xb;params(qb)]

    axis = srand(3)
    vert = (srand(3),srand(3))

    body1 = Box(1., 1., 1., 1.)
    body2 = Box(1., 1., 1., 1.)
    trans = Translational{Float64,2}(body1, body2; axis = axis, p1 = vert[1], p2 = vert[2])[1]

    XXaa, XQaa, QXaa, QQaa = ∂2g∂posaa(trans, xa, qa, xb, qb)
    XXab, XQab, QXab, QQab = ∂2g∂posab(trans, xa, qa, xb, qb)
    XXba, XQba, QXba, QQba = ∂2g∂posba(trans, xa, qa, xb, qb)
    XXbb, XQbb, QXbb, QQbb = ∂2g∂posbb(trans, xa, qa, xb, qb)

    G = [
        XXaa XQaa XXab XQab
        QXaa QQaa QXab QQab
        XXba XQba XXbb XQbb
        QXba QQba QXbb QQbb
    ]

    function firstderivative(z)
        E = diagm(ones(3))
        Z34 = zeros(3,4)
        Z43 = zeros(4,3)
        postmultmat = [
            E 0*I 0*I 0*I
            Z43 LVᵀmat(Quaternion(z[4:7]...)) Z43 Z43
            0*I 0*I E 0*I
            Z43 Z43 Z43 LVᵀmat(Quaternion(z[11:14]...))
        ]
        ForwardDiff.jacobian(x->g(trans,x[1:3],Quaternion(x[4:7]...),x[8:10],Quaternion(x[11:14]...)),z)*postmultmat
    end
    function secondderivative(z)
        ForwardDiff.jacobian(x->vec(firstderivative(x)),z)
    end
    @test isapprox(norm(G-secondderivative(z)), 0.0; atol = 1e-8)
end

function rottest()
    xa = srand(3)
    qa = rand(QuatRotation).q
    xb = srand(3)
    qb = rand(QuatRotation).q
    z = [xa;params(qa);xb;params(qb)]

    axis = srand(3)
    qoffset = rand(QuatRotation).q

    body1 = Box(1., 1., 1., 1.)
    body2 = Box(1., 1., 1., 1.)
    rot = Rotational{Float64,2}(body1, body2; axis = axis, qoffset = qoffset)[1]

    XXaa, XQaa, QXaa, QQaa = ∂2g∂posaa(rot, xa, qa, xb, qb)
    XXab, XQab, QXab, QQab = ∂2g∂posab(rot, xa, qa, xb, qb)
    XXba, XQba, QXba, QQba = ∂2g∂posba(rot, xa, qa, xb, qb)
    XXbb, XQbb, QXbb, QQbb = ∂2g∂posbb(rot, xa, qa, xb, qb)

    G = [
        XXaa XQaa XXab XQab
        QXaa QQaa QXab QQab
        XXba XQba XXbb XQbb
        QXba QQba QXbb QQbb
    ]

    function firstderivative(z)
        E = diagm(ones(3))
        Z34 = zeros(3,4)
        Z43 = zeros(4,3)
        postmultmat = [
            E 0*I 0*I 0*I
            Z43 LVᵀmat(Quaternion(z[4:7]...)) Z43 Z43
            0*I 0*I E 0*I
            Z43 Z43 Z43 LVᵀmat(Quaternion(z[11:14]...))
        ]
        ForwardDiff.jacobian(x->g(rot,x[1:3],Quaternion(x[4:7]...),x[8:10],Quaternion(x[11:14]...)),z)*postmultmat
    end
    function secondderivative(z)
        ForwardDiff.jacobian(x->vec(firstderivative(x)),z)
    end
    @test isapprox(norm(G-secondderivative(z)), 0.0; atol = 1e-8)
end

for i=1:10
    transtest()
    rottest()
end
