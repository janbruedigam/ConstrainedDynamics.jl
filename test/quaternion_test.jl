using ConstrainedDynamics
using ConstrainedDynamics: params, Lmat, Lᵀmat, Rmat, Rᵀmat, Tmat, Tᵀmat, Vmat, Vᵀmat, VLmat, VLᵀmat, VRmat, VRᵀmat, LVᵀmat, LᵀVᵀmat, RVᵀmat, RᵀVᵀmat, slerp
using StaticArrays
using LinearAlgebra

p1 = normalize(rand(4))
w1 = p1[1]
v1 = p1[2:4]
p2 = normalize(rand(4))
w2 = p2[1]
v2 = p2[2:4]
qvref1 = Quaternion(v1)
qref1 = Quaternion(w1,v1...)
qref2 = Quaternion(w2,v2...)
q1 = Quaternion(w1,v1...)
@test q1 == qref1
q1 = Quaternion(w1,v1...)
@test q1 == qref1
qv1 = Quaternion(v1)
@test qv1 == qvref1
qv1 = Quaternion(v1)
@test qv1 == qvref1


@test isapprox(norm(params(qref1*qref2) - Lmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref1\qref2) - Lᵀmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2*qref1) - Rmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2/qref1) - Rᵀmat(q1)*params(qref2)), 0.0; atol = 1e-10)

@test isapprox(norm(params(inv(qref1)) - Tmat()*params(q1)), 0.0; atol = 1e-10)
@test isapprox(norm(params(inv(qref1)) - Tᵀmat()*params(q1)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qvref1)[2:4] - Vmat()*params(qv1)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qvref1) - Vᵀmat()*v1), 0.0; atol = 1e-10)

@test isapprox(norm(params(qref1*qref2)[2:4] - VLmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref1\qref2)[2:4] - VLᵀmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2*qref1)[2:4] - VRmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2/qref1)[2:4] - VRᵀmat(q1)*params(qref2)), 0.0; atol = 1e-10)

@test isapprox(norm(params(qref1*qvref1) - LVᵀmat(q1)*v1), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref1\qvref1) - LᵀVᵀmat(q1)*v1), 0.0; atol = 1e-10)
@test isapprox(norm(params(qvref1*qref1) - RVᵀmat(q1)*v1), 0.0; atol = 1e-10)
@test isapprox(norm(params(qvref1/qref1) - RᵀVᵀmat(q1)*v1), 0.0; atol = 1e-10)

@test isapprox(norm(params(qref1) - params(slerp(qref1,qref2,0.0))), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2) - params(slerp(qref1,qref2,1.0))), 0.0; atol = 1e-10)