using LinearAlgebra
using RigidBodyDynamics
using StaticArrays

g = -9.81

l = 1.
m = l
I = 0.0841667

joint_axis = SVector(1., 0., 0.);

N = 3

origin = RigidBody{Float64}("origin")
chain = Mechanism(origin; gravity = SVector(0., 0., g))

joint1 = Joint("joint1", Revolute(joint_axis))
inertia1 = SpatialInertia(frame_after(joint1),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
link1 = RigidBody(inertia1)
before_joint1_to_origin = one(Transform3D,frame_before(joint1), default_frame(origin))
attach!(chain, origin, link1, joint1,joint_pose = before_joint1_to_origin)

joint2 = Joint("joint2", Revolute(joint_axis))
inertia2 = SpatialInertia(frame_after(joint2),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
link2 = RigidBody(inertia2)
before_joint2_to_after_joint1 = Transform3D(frame_before(joint2), frame_after(joint1), SVector(0, 0., -l))
attach!(chain, link1, link2, joint2,joint_pose = before_joint2_to_after_joint1)

joint3 = Joint("joint3", Revolute(joint_axis))
inertia3 = SpatialInertia(frame_after(joint3),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
link3 = RigidBody(inertia3)
before_joint3_to_origin = one(Transform3D,frame_before(joint3), default_frame(origin))
attach!(chain, origin, link3, joint3,joint_pose = before_joint3_to_origin)

joint4 = Joint("joint4", Revolute(joint_axis))
inertia4 = SpatialInertia(frame_after(joint4),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
link4 = RigidBody(inertia4)
before_joint4_to_after_joint3 = Transform3D(frame_before(joint4), frame_after(joint3), SVector(0, 0., -l))
attach!(chain, link3, link4, joint4,joint_pose = before_joint4_to_after_joint3)


for i=2:N
    @eval begin
        $(Symbol("joint",(i-1)*4+1)) = Joint("joint"*string(($i-1)*4+1), Revolute(joint_axis))
        $(Symbol("inertia",(i-1)*4+1)) = SpatialInertia(frame_after($(Symbol("joint",(i-1)*4+1))),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
        $(Symbol("link",(i-1)*4+1)) = RigidBody($(Symbol("inertia",(i-1)*4+1)))
        $(Symbol("before_joint",(i-1)*4+1,"_to_after_joint",(i-2)*4+2)) = Transform3D(frame_before($(Symbol("joint",(i-1)*4+1))), frame_after($(Symbol("joint",(i-2)*4+2))), SVector(0, 0., -l))
        attach!(chain, $(Symbol("link",(i-2)*4+2)), $(Symbol("link",(i-1)*4+1)), $(Symbol("joint",(i-1)*4+1)),joint_pose = $(Symbol("before_joint",(i-1)*4+1,"_to_after_joint",(i-2)*4+2)))

        $(Symbol("joint",(i-1)*4+2)) = Joint("joint"*string(($i-1)*4+2), Revolute(joint_axis))
        $(Symbol("inertia",(i-1)*4+2)) = SpatialInertia(frame_after($(Symbol("joint",(i-1)*4+2))),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
        $(Symbol("link",(i-1)*4+2)) = RigidBody($(Symbol("inertia",(i-1)*4+2)))
        $(Symbol("before_joint",(i-1)*4+2,"_to_after_joint",(i-1)*4+1)) = Transform3D(frame_before($(Symbol("joint",(i-1)*4+2))), frame_after($(Symbol("joint",(i-1)*4+1))), SVector(0, 0., -l))
        attach!(chain, $(Symbol("link",(i-1)*4+1)), $(Symbol("link",(i-1)*4+2)), $(Symbol("joint",(i-1)*4+2)),joint_pose = $(Symbol("before_joint",(i-1)*4+2,"_to_after_joint",(i-1)*4+1)))

        $(Symbol("joint",(i-1)*4+3)) = Joint("joint"*string(($i-1)*4+3), Revolute(joint_axis))
        $(Symbol("inertia",(i-1)*4+3)) = SpatialInertia(frame_after($(Symbol("joint",(i-1)*4+3))),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
        $(Symbol("link",(i-1)*4+3)) = RigidBody($(Symbol("inertia",(i-1)*4+3)))
        $(Symbol("before_joint",(i-1)*4+3,"_to_after_joint",(i-2)*4+2)) = Transform3D(frame_before($(Symbol("joint",(i-1)*4+3))), frame_after($(Symbol("joint",(i-2)*4+2))), SVector(0, 0., -l))
        attach!(chain, $(Symbol("link",(i-2)*4+2)), $(Symbol("link",(i-1)*4+3)), $(Symbol("joint",(i-1)*4+3)),joint_pose = $(Symbol("before_joint",(i-1)*4+3,"_to_after_joint",(i-2)*4+2)))

        $(Symbol("joint",(i-1)*4+4)) = Joint("joint"*string(($i-1)*4+4), Revolute(joint_axis))
        $(Symbol("inertia",(i-1)*4+4)) = SpatialInertia(frame_after($(Symbol("joint",(i-1)*4+4))),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
        $(Symbol("link",(i-1)*4+4)) = RigidBody($(Symbol("inertia",(i-1)*4+4)))
        $(Symbol("before_joint",(i-1)*4+4,"_to_after_joint",(i-1)*4+3)) = Transform3D(frame_before($(Symbol("joint",(i-1)*4+4))), frame_after($(Symbol("joint",(i-1)*4+3))), SVector(0, 0., -l))
        attach!(chain, $(Symbol("link",(i-1)*4+3)), $(Symbol("link",(i-1)*4+4)), $(Symbol("joint",(i-1)*4+4)),joint_pose = $(Symbol("before_joint",(i-1)*4+4,"_to_after_joint",(i-1)*4+3)))
    end
end

for i=1:N
    @eval begin
        $(Symbol("joint",N*4+i)) = Joint("joint"*string(N*4+$i), Revolute(joint_axis))
        $(Symbol("before_joint",N*4+i,"_to_joint",(i-1)*4+2)) = Transform3D(frame_before($(Symbol("joint",N*4+i))), frame_after($(Symbol("joint",(i-1)*4+2))), SVector(0, 0., -l))
        $(Symbol("joint",(i-1)*4+4,"_to_after_joint",N*4+i)) = Transform3D(frame_after($(Symbol("joint",(i-1)*4+4))), frame_after($(Symbol("joint",N*4+i))), SVector(0, 0., l))
        attach!(chain, $(Symbol("link",(i-1)*4+2)), $(Symbol("link",(i-1)*4+4)), $(Symbol("joint",N*4+i)),joint_pose = $(Symbol("before_joint",N*4+i,"_to_joint",(i-1)*4+2)), successor_pose = $(Symbol("joint",(i-1)*4+4,"_to_after_joint",N*4+i)))
    end
end


state = MechanismState(chain)

offset = π/4

if N==1
    set_configuration!(state, joint1, π/8+offset)
    set_configuration!(state, joint2, -2*π/8)
    set_configuration!(state, joint3, -π/8+offset)
    set_configuration!(state, joint4, 2*π/8)
else
    set_configuration!(state, joint1, π/8)
    set_configuration!(state, joint2, -2*π/8)
    set_configuration!(state, joint3, -π/8)
    set_configuration!(state, joint4, 2*π/8)
end

for i=2:N-1
    @eval begin
        set_configuration!(state, $(Symbol("joint",(i-1)*4+1)), 2*π/8)
        set_configuration!(state, $(Symbol("joint",(i-1)*4+2)), -2*π/8)
        set_configuration!(state, $(Symbol("joint",(i-1)*4+3)), 0)
        set_configuration!(state, $(Symbol("joint",(i-1)*4+4)), 2*π/8)
    end
end

if N>1
    for i=N
        @eval begin
            set_configuration!(state, $(Symbol("joint",(i-1)*4+1)), 2*π/8+offset)
            set_configuration!(state, $(Symbol("joint",(i-1)*4+2)), -2*π/8)
            set_configuration!(state, $(Symbol("joint",(i-1)*4+3)), 0+offset)
            set_configuration!(state, $(Symbol("joint",(i-1)*4+4)), 2*π/8)
        end
    end
end

setdirty!(state)

ts, qs, vs = simulate(state, 10., Δt = 0.01);

for i=1:length(qs)
    if isnan(qs[i][1])
        display(string("Failure after ",i," time steps"))
        break;
    end
end
