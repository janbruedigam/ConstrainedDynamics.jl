using LinearAlgebra
using RigidBodyDynamics
using StaticArrays
using Blink

g = -9.81

l1 = 1.
m1 = l1
I1 = 0.0841667

axis = SVector(1., 0., 0.);

N = 3

world = RigidBody{Float64}("world")
fourbar = Mechanism(world; gravity = SVector(0., 0., g))

joint1 = Joint("joint1", Revolute(axis))
inertia1 = SpatialInertia(frame_after(joint1),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
link1 = RigidBody(inertia1)
before_joint1_to_world = one(Transform3D,frame_before(joint1), default_frame(world))
attach!(fourbar, world, link1, joint1,joint_pose = before_joint1_to_world)

joint2 = Joint("joint2", Revolute(axis))
inertia2 = SpatialInertia(frame_after(joint2),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
link2 = RigidBody(inertia2)
before_joint2_to_after_joint1 = Transform3D(frame_before(joint2), frame_after(joint1), SVector(0, 0., -l1))
attach!(fourbar, link1, link2, joint2,joint_pose = before_joint2_to_after_joint1)

joint3 = Joint("joint3", Revolute(axis))
inertia3 = SpatialInertia(frame_after(joint3),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
link3 = RigidBody(inertia3)
before_joint3_to_world = one(Transform3D,frame_before(joint3), default_frame(world))
attach!(fourbar, world, link3, joint3,joint_pose = before_joint3_to_world)

joint4 = Joint("joint4", Revolute(axis))
inertia4 = SpatialInertia(frame_after(joint4),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
link4 = RigidBody(inertia4)
before_joint4_to_after_joint3 = Transform3D(frame_before(joint4), frame_after(joint3), SVector(0, 0., -l1))
attach!(fourbar, link3, link4, joint4,joint_pose = before_joint4_to_after_joint3)


for i=2:N
    @eval begin
        $(Symbol("joint",(i-1)*4+1)) = Joint("joint"*string(($i-1)*4+1), Revolute(axis))
        $(Symbol("inertia",(i-1)*4+1)) = SpatialInertia(frame_after($(Symbol("joint",(i-1)*4+1))),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
        $(Symbol("link",(i-1)*4+1)) = RigidBody($(Symbol("inertia",(i-1)*4+1)))
        $(Symbol("before_joint",(i-1)*4+1,"_to_after_joint",(i-2)*4+2)) = Transform3D(frame_before($(Symbol("joint",(i-1)*4+1))), frame_after($(Symbol("joint",(i-2)*4+2))), SVector(0, 0., -l1))
        attach!(fourbar, $(Symbol("link",(i-2)*4+2)), $(Symbol("link",(i-1)*4+1)), $(Symbol("joint",(i-1)*4+1)),joint_pose = $(Symbol("before_joint",(i-1)*4+1,"_to_after_joint",(i-2)*4+2)))

        $(Symbol("joint",(i-1)*4+2)) = Joint("joint"*string(($i-1)*4+2), Revolute(axis))
        $(Symbol("inertia",(i-1)*4+2)) = SpatialInertia(frame_after($(Symbol("joint",(i-1)*4+2))),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
        $(Symbol("link",(i-1)*4+2)) = RigidBody($(Symbol("inertia",(i-1)*4+2)))
        $(Symbol("before_joint",(i-1)*4+2,"_to_after_joint",(i-1)*4+1)) = Transform3D(frame_before($(Symbol("joint",(i-1)*4+2))), frame_after($(Symbol("joint",(i-1)*4+1))), SVector(0, 0., -l1))
        attach!(fourbar, $(Symbol("link",(i-1)*4+1)), $(Symbol("link",(i-1)*4+2)), $(Symbol("joint",(i-1)*4+2)),joint_pose = $(Symbol("before_joint",(i-1)*4+2,"_to_after_joint",(i-1)*4+1)))

        $(Symbol("joint",(i-1)*4+3)) = Joint("joint"*string(($i-1)*4+3), Revolute(axis))
        $(Symbol("inertia",(i-1)*4+3)) = SpatialInertia(frame_after($(Symbol("joint",(i-1)*4+3))),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
        $(Symbol("link",(i-1)*4+3)) = RigidBody($(Symbol("inertia",(i-1)*4+3)))
        $(Symbol("before_joint",(i-1)*4+3,"_to_after_joint",(i-2)*4+2)) = Transform3D(frame_before($(Symbol("joint",(i-1)*4+3))), frame_after($(Symbol("joint",(i-2)*4+2))), SVector(0, 0., -l1))
        attach!(fourbar, $(Symbol("link",(i-2)*4+2)), $(Symbol("link",(i-1)*4+3)), $(Symbol("joint",(i-1)*4+3)),joint_pose = $(Symbol("before_joint",(i-1)*4+3,"_to_after_joint",(i-2)*4+2)))

        $(Symbol("joint",(i-1)*4+4)) = Joint("joint"*string(($i-1)*4+4), Revolute(axis))
        $(Symbol("inertia",(i-1)*4+4)) = SpatialInertia(frame_after($(Symbol("joint",(i-1)*4+4))),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
        $(Symbol("link",(i-1)*4+4)) = RigidBody($(Symbol("inertia",(i-1)*4+4)))
        $(Symbol("before_joint",(i-1)*4+4,"_to_after_joint",(i-1)*4+3)) = Transform3D(frame_before($(Symbol("joint",(i-1)*4+4))), frame_after($(Symbol("joint",(i-1)*4+3))), SVector(0, 0., -l1))
        attach!(fourbar, $(Symbol("link",(i-1)*4+3)), $(Symbol("link",(i-1)*4+4)), $(Symbol("joint",(i-1)*4+4)),joint_pose = $(Symbol("before_joint",(i-1)*4+4,"_to_after_joint",(i-1)*4+3)))
    end
end

for i=1:N
    @eval begin
        $(Symbol("joint",N*4+i)) = Joint("joint"*string(N*4+$i), Revolute(axis))
        $(Symbol("before_joint",N*4+i,"_to_joint",(i-1)*4+2)) = Transform3D(frame_before($(Symbol("joint",N*4+i))), frame_after($(Symbol("joint",(i-1)*4+2))), SVector(0, 0., -l1))
        $(Symbol("joint",(i-1)*4+4,"_to_after_joint",N*4+i)) = Transform3D(frame_after($(Symbol("joint",(i-1)*4+4))), frame_after($(Symbol("joint",N*4+i))), SVector(0, 0., l1))
        attach!(fourbar, $(Symbol("link",(i-1)*4+2)), $(Symbol("link",(i-1)*4+4)), $(Symbol("joint",N*4+i)),joint_pose = $(Symbol("before_joint",N*4+i,"_to_joint",(i-1)*4+2)), successor_pose = $(Symbol("joint",(i-1)*4+4,"_to_after_joint",N*4+i)))
    end
end


state = MechanismState(fourbar)
result = DynamicsResult(fourbar);

offset = pi/4

if N==1
    set_configuration!(state, joint1, pi/8+offset)
    set_configuration!(state, joint2, -2*pi/8)
    set_configuration!(state, joint3, -pi/8+offset)
    set_configuration!(state, joint4, 2*pi/8)
else
    set_configuration!(state, joint1, pi/8)
    set_configuration!(state, joint2, -2*pi/8)
    set_configuration!(state, joint3, -pi/8)
    set_configuration!(state, joint4, 2*pi/8)
end

for i=2:N-1
    @eval begin
        set_configuration!(state, $(Symbol("joint",(i-1)*4+1)), 2*pi/8)
        set_configuration!(state, $(Symbol("joint",(i-1)*4+2)), -2*pi/8)
        set_configuration!(state, $(Symbol("joint",(i-1)*4+3)), 0)
        set_configuration!(state, $(Symbol("joint",(i-1)*4+4)), 2*pi/8)
    end
end

if N>1
    for i=N
        @eval begin
            set_configuration!(state, $(Symbol("joint",(i-1)*4+1)), 2*pi/8+offset)
            set_configuration!(state, $(Symbol("joint",(i-1)*4+2)), -2*pi/8)
            set_configuration!(state, $(Symbol("joint",(i-1)*4+3)), 0+offset)
            set_configuration!(state, $(Symbol("joint",(i-1)*4+4)), 2*pi/8)
        end
    end
end

# set_configuration!(state, joint5, 0-offset)
# set_configuration!(state, joint6, 2*pi/8)
# set_configuration!(state, joint7, 2*pi/8-offset)
# set_configuration!(state, joint8, -2*pi/8)
#
# set_configuration!(state, joint5, 0)
# set_configuration!(state, joint6, 2*pi/8)
# set_configuration!(state, joint7, 2*pi/8)
# set_configuration!(state, joint8, -2*pi/8)
#
# set_configuration!(state, joint9, 0+offset)
# set_configuration!(state, joint10, -2*pi/8)
# set_configuration!(state, joint11, -2*pi/8+offset)
# set_configuration!(state, joint12, 2*pi/8)

setdirty!(state)

ts, qs, vs = simulate(state, 10., Δt = 0.01);


using MeshCatMechanisms

urdf = "paperexamples/3fourbars.urdf"

mvis = MechanismVisualizer(fourbar, URDFVisuals(urdf));
open(mvis, Blink.Window())
MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 0.1);

for i=1:length(qs)
    if isnan(qs[i][1])
        display(i)
        break;
    end
end
