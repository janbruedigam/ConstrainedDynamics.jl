using LinearAlgebra
using RigidBodyDynamics
using StaticArrays
using Blink

g = -9.81

l1 = 1.
l2 = sqrt(2)/2
m1 = l1
m2 = l2
I1 = 0.0841667
I2 = 0.030052

axis = SVector(1., 0., 0.);

world = RigidBody{Float64}("world")
fourbar = Mechanism(world; gravity = SVector(0., 0., g))

joint1 = Joint("joint1", Revolute(axis))
inertia1 = SpatialInertia(frame_after(joint1),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
link1 = RigidBody(inertia1)
before_joint1_to_world = one(Transform3D,frame_before(joint1), default_frame(world))
attach!(fourbar, world, link1, joint1,joint_pose = before_joint1_to_world)

joint2 = Joint("joint2", Revolute(axis))
inertia2 = SpatialInertia(frame_after(joint2),com=SVector(0, 0, -l2/2),moment_about_com=I2*axis*transpose(axis),mass=m2)
link2 = RigidBody(inertia2)
before_joint2_to_after_joint1 = Transform3D(frame_before(joint2), frame_after(joint1), SVector(0, 0., -l1))
attach!(fourbar, link1, link2, joint2,joint_pose = before_joint2_to_after_joint1)

joint3 = Joint("joint3", Revolute(axis))
inertia3 = SpatialInertia(frame_after(joint3),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
link3 = RigidBody(inertia3)
before_joint3_to_world = one(Transform3D,frame_before(joint3), default_frame(world))
attach!(fourbar, world, link3, joint3,joint_pose = before_joint3_to_world)

joint4 = Joint("joint4", Revolute(axis))
inertia4 = SpatialInertia(frame_after(joint4),com=SVector(0, 0, -l2/2),moment_about_com=I2*axis*transpose(axis),mass=m2)
link4 = RigidBody(inertia4)
before_joint4_to_after_joint3 = Transform3D(frame_before(joint4), frame_after(joint3), SVector(0, 0., -l1))
attach!(fourbar, link3, link4, joint4,joint_pose = before_joint4_to_after_joint3)

joint5 = Joint("joint5", Revolute(axis))
before_joint5_to_joint2 = Transform3D(frame_before(joint5), frame_after(joint2), SVector(0, 0., -l2))
joint4_to_after_joint5 = Transform3D(frame_after(joint4), frame_after(joint5), SVector(0, 0., l2))
attach!(fourbar, link2, link4, joint5,joint_pose = before_joint5_to_joint2, successor_pose = joint4_to_after_joint5)

state = MechanismState(fourbar)
result = DynamicsResult(fourbar);

set_configuration!(state, joint1, pi/2)
set_configuration!(state, joint2, -3*pi/4)
set_configuration!(state, joint3, 0.)
set_configuration!(state, joint4, 3*pi/4)

setdirty!(state)

ts, qs, vs = simulate(state, 10., Δt = 0.0001);


using MeshCatMechanisms

urdf = "paperexamples/twoTwoBarDiffLength.urdf"

mvis = MechanismVisualizer(fourbar, URDFVisuals(urdf));
open(mvis, Blink.Window())
MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 1.);
