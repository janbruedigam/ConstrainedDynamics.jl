using LinearAlgebra
using RigidBodyDynamics
using StaticArrays
using NLsolve


l1 = 1. # length 1
l2 = sqrt(2)/2 # length 2
m1 = l1 # mass 1
m2 = l2 # mass 2
I1 = 0.0841667 # Inertia 1
I2 = 0.030052 # Inertia 2

# Set up mechanism

joint_axis = SVector(1., 0., 0.);

origin = RigidBody{Float64}("origin")
threebar = Mechanism(origin; gravity = SVector(0., 0., -9.81))

joint1 = Joint("joint1", Revolute(joint_axis))
inertia1 = SpatialInertia(frame_after(joint1),com=SVector(0, 0, -l1/2),moment_about_com=I1*joint_axis*transpose(joint_axis),mass=m1)
link1 = RigidBody(inertia1)
before_joint1_to_origin = one(Transform3D,frame_before(joint1), default_frame(origin))
attach!(threebar, origin, link1, joint1,joint_pose = before_joint1_to_origin)

joint2 = Joint("joint2", Revolute(joint_axis))
inertia2 = SpatialInertia(frame_after(joint2),com=SVector(0, 0, -l2/2),moment_about_com=I2*joint_axis*transpose(joint_axis),mass=m2)
link2 = RigidBody(inertia2)
before_joint2_to_after_joint1 = Transform3D(frame_before(joint2), frame_after(joint1), SVector(0, 0., -l1))
attach!(threebar, link1, link2, joint2,joint_pose = before_joint2_to_after_joint1)

joint3 = Joint("joint3", Revolute(joint_axis))
inertia3 = SpatialInertia(frame_after(joint3),com=SVector(0, 0, -l1/2),moment_about_com=I1*joint_axis*transpose(joint_axis),mass=m1)
link3 = RigidBody(inertia3)
before_joint3_to_origin = Transform3D(frame_before(joint3), default_frame(origin), SVector(0, l1,0))
attach!(threebar, origin, link3, joint3,joint_pose = before_joint3_to_origin)

joint4 = Joint("joint4", Revolute(joint_axis))
before_joint4_to_joint2 = Transform3D(frame_before(joint4), frame_after(joint2), SVector(0, 0., -l2))
joint3_to_after_joint4 = Transform3D(frame_after(joint3), frame_after(joint4), SVector(0, 0., l1))
attach!(threebar, link2, link3, joint4,joint_pose = before_joint4_to_joint2, successor_pose = joint3_to_after_joint4)

state = MechanismState(threebar)

# Find joint coordinates for loop closure
a3 = pi
function f!(F, x)
    F[1] = l1*sin(x[1])+l2*sin(x[1]+x[2])-l1-sin(a3)*l1
    F[2] = -l1*cos(x[1])-l2*cos(x[1]+x[2])+cos(a3)*l1
end
sol = nlsolve(f!, [a3;sign(cos(a3))*pi/2]).zero
set_configuration!(state, joint1, sol[1])
set_configuration!(state, joint2, sol[2])
set_configuration!(state, joint3, a3)


setdirty!(state)

ts, qs, vs = simulate(state, 600., Δt = 0.01, stabilization_gains=nothing);

drift = zeros(length(ts))
for (i,el) in enumerate(qs)
    p1x = l1*sin(el[1])+l2*sin(el[1]+el[2])
    p1y = -l1*cos(el[1])-l2*cos(el[1]+el[2])
    p2x = l1+sin(el[3])*l1
    p2y = -cos(el[3])*l1
    drift[i] = norm([p1x;p1y]-[p2x;p2y])
end
