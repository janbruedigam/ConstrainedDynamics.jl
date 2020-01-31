using LinearAlgebra
using RigidBodyDynamics
using StaticArrays
using Blink

g = -9.81

l1 = 1.
m1 = l1
I1 = 0.0841667

axis = SVector(1., 0., 0.);

function test(state,jointies,N)
    ang = 0.2
    set_configuration!(state, jointies[1], ang)
    for i = 2:floor(Int64,N/2)
        set_configuration!(state, jointies[i], 0.)
    end
    set_configuration!(state, jointies[ceil(Int64,N/2)], pi/2-ang)
    set_configuration!(state, jointies[ceil(Int64,N/2)+1], pi/2+ang)
    for i = ceil(Int64,N/2)+2:N-1
        set_configuration!(state, jointies[i], 0.)
    end
    set_configuration!(state, jointies[N], ang)

    setdirty!(state)
    simulate(state, 10., Δt = 0.1)
end

world = RigidBody{Float64}("world")
fourbar = Mechanism(world; gravity = SVector(0., 0., g))

N = 4

joint1 = Joint("joint1", Revolute(axis))
inertia1 = SpatialInertia(frame_after(joint1),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
link1 = RigidBody(inertia1)
before_joint1_to_world = one(Transform3D,frame_before(joint1), default_frame(world))
attach!(fourbar, world, link1, joint1,joint_pose = before_joint1_to_world)

links = [link1]
jointies = [joint1]

for i=2:N-1
    @eval begin
        $(Symbol("joint",i)) = Joint("joint"*string($i), Revolute(axis))
        $(Symbol("inertia",i)) = SpatialInertia(frame_after($(Symbol("joint",i))),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
        $(Symbol("link",i)) = RigidBody($(Symbol("inertia",i)))
        $(Symbol("before_joint",i,"_to_after_joint",i-1)) = Transform3D(frame_before($(Symbol("joint",i))), frame_after($(Symbol("joint",i-1))), SVector(0, 0., -l1))
        attach!(fourbar, $(Symbol("link",i-1)), $(Symbol("link",i)), $(Symbol("joint",i)),joint_pose = $(Symbol("before_joint",i,"_to_after_joint",i-1)))
        push!(links,$(Symbol("link",i)))
        push!(jointies,$(Symbol("joint",i)))

    end
end

@eval begin
    $(Symbol("joint",N)) = Joint("joint"*string(N), Revolute(axis))
    $(Symbol("inertia",N)) = SpatialInertia(frame_after($(Symbol("joint",N))),com=SVector(0, 0, -l1/2),moment_about_com=I1*axis*transpose(axis),mass=m1)
    $(Symbol("link",N)) = RigidBody($(Symbol("inertia",N)))
    $(Symbol("before_joint",N,"_to_world")) = one(Transform3D,frame_before($(Symbol("joint",N))), default_frame(world))
    attach!(fourbar, world, $(Symbol("link",N)), $(Symbol("joint",N)),joint_pose = $(Symbol("before_joint",N,"_to_world")))
    push!(links,$(Symbol("link",N)))
    push!(jointies,$(Symbol("joint",N)))

end

@eval begin
    $(Symbol("joint",N+1)) = Joint("joint"*string(N+1), Revolute(axis))
    $(Symbol("before_joint",N+1,"_to_joint",N-1)) = Transform3D(frame_before($(Symbol("joint",N+1))), frame_after(jointies[N-1]), SVector(0, 0., -l1))
    $(Symbol("joint",N,"_to_joint",N+1)) = Transform3D(frame_after(jointies[N]), frame_after($(Symbol("joint",N+1))), SVector(0, 0., l1))
    attach!(fourbar, links[N-1], links[N], $(Symbol("joint",N+1)),joint_pose = $(Symbol("before_joint",N+1,"_to_joint",N-1)), successor_pose = $(Symbol("joint",N,"_to_joint",N+1)))
    push!(jointies,$(Symbol("joint",N+1)))

end

state = MechanismState(fourbar)
result = DynamicsResult(fourbar);


phi = pi/4
offset = pi/4
set_configuration!(state, jointies[1], phi+offset)
set_configuration!(state, jointies[2], -2*phi)
set_configuration!(state, jointies[3], -pi+2*phi)
set_configuration!(state, jointies[4], -phi+offset)
# set_configuration!(state, jointies[1], phi+offset)
# set_configuration!(state, jointies[2], -phi)
# set_configuration!(state, jointies[3], 0.)
# set_configuration!(state, jointies[25], -phi)
# set_configuration!(state, jointies[26], -pi+2*phi)
# set_configuration!(state, jointies[27], -phi)
# set_configuration!(state, jointies[13], 0.)
# set_configuration!(state, jointies[50], -phi+offset)

setdirty!(state)

# ts, qs, vs = test(state,jointies,N);
ts, qs, vs = simulate(state, 10., Δt = 0.01);

drift = zeros(length(ts))
for (i,el) in enumerate(qs)
    for (j,el2) in enumerate(el)
        set_configuration!(state, jointies[j], el2)
    end
    p1 = transform_to_root(state, frame_after(jointies[3])).mat*[0;0;-l1;1]
    p2 = transform_to_root(state, frame_after(jointies[4])).mat*[0;0;-l1;1]
    drift[i] = norm(p2-p1)
end




using MeshCatMechanisms

urdf = "paperexamples/fourbarBar.urdf"
mvis = MechanismVisualizer(fourbar, URDFVisuals(urdf));
# open(mvis)
open(mvis, Blink.Window())
MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 1.);
